#!/usr/bin/env python3
"""
Accuracy check on real genotypes: the mtDNA of the 2,504 1000 Genomes phase 3
samples (mtDNA-Server calls from the haplocheck paper, MIT), compared with the
haplogroups HaploGrep 2 called from the same data on PhyloTree 17.

Two densities: the whole mitochondrial genome, and the 326 mtDNA positions of
the Illumina Omni 2.5 array (taken from HaploGrep 3's microarray example).
--errors adds random genotyping errors to the array calls; those runs are
reported without a pass mark.

Downloads about 1 MB from GitHub on first run (cached in tests/.cache).
"""

import argparse
import collections
import csv
import gzip
import importlib.util
import os
import random
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CACHE = os.path.join(HERE, ".cache")
HAPLOCHECK = "https://raw.githubusercontent.com/genepi/haplocheck/master/paper-data/1000g/"
CALLS = HAPLOCHECK + "results/mutserve.vcf.gz"
TRUTH = HAPLOCHECK + "1000g_haplogroups.csv"
OMNI = "https://raw.githubusercontent.com/genepi/haplogrep3/main/data/examples/example-microarray.vcf"

spec = importlib.util.spec_from_file_location("mhp", os.path.join(ROOT, "mt-haplogroup-predictor.py"))
mhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mhp)


def fetch(url):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, url.rsplit("/", 1)[1])
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    return path


def read_calls():
    """Each sample's non-reference calls, read with the predictor's own VCF rules."""
    samples = []
    with gzip.open(fetch(CALLS), "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = [mhp.Sample(name.split(".")[0]) for name in fields[9:]]
                continue
            pos = int(fields[1])
            ref = fields[3].upper()
            alts = fields[4].upper().split(",")
            fmt = fields[8].split(":")
            for s, value in zip(samples, fields[9:]):
                if value == "0" or value.startswith("0:"):
                    continue
                values = value.split(":")
                idx = mhp.vcf_alleles(values[0])
                if idx is None:
                    s.add(pos, None, "nocall")
                    continue
                pick = idx[0] if len(set(idx)) == 1 else mhp.major_allele(idx, fmt, values, len(alts), False)
                if pick is None:
                    s.add(pos, None, "het")
                elif pick:
                    mhp.add_vcf_allele(s, pos, ref, ([ref] + alts)[pick])
    return samples


def full_sequence(s, db):
    t = mhp.Sample(s.path)
    t.kind = "vcf"
    t.by_pos = dict(s.by_pos)
    t.dropped = set(s.dropped)
    t.blocked = set(s.blocked)
    t.ref_bases = {p: db.rcrs[p - 1] for p in s.by_pos}
    t.fill = True
    return t


def array(s, db, panel, error, rng):
    t = mhp.Sample(s.path)
    for p in panel:
        if p in s.dropped or p in s.blocked:
            continue
        base = s.by_pos.get(p, db.rcrs[p - 1])
        if error and rng.random() < error:
            base = rng.choice([b for b in "ACGT" if b != base])
        t.add(p, base, "call")
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--errors", default="", help="comma-separated genotyping error rates for extra array runs")
    ap.add_argument("--skip-full", action="store_true", help="only run the array comparison (faster)")
    args = ap.parse_args()

    db = mhp.TreeDB(os.path.join(ROOT, "mtdb.tsv.gz"))
    truth = {r["SampleID"]: r["HG mutserve"] for r in csv.DictReader(open(fetch(TRUTH)))}
    panel = sorted({int(line.split("\t")[1]) for line in open(fetch(OMNI)) if not line.startswith("#")})
    samples = read_calls()
    by_name = collections.defaultdict(set)
    for i, nid in enumerate(db.ids):
        by_name[nid.upper()].add(i)

    class Args:
        error_rate = 0.01
        recurrence_rate = None
        max_error_rate = 0.4
        min_support = 0.95
        coordinates = "auto"
        unlisted = "auto"
        complement = "auto"
        exclude = set()

    runs = [] if args.skip_full else [("whole mtDNA genome", None)]
    runs.append((f"Omni 2.5 array, {len(panel)} positions", 0.0))
    runs.extend((f"Omni 2.5 array, {float(e):.0%} genotype errors", float(e)) for e in args.errors.split(",") if e)
    ok = True
    for label, error in runs:
        rng = random.Random(20261010)
        tally = collections.Counter()
        misses = []
        for s in samples:
            sample = full_sequence(s, db) if error is None else array(s, db, panel, error, rng)
            r = mhp.analyse(sample, db, Args)
            expected = truth.get(s.path, "")
            nodes = by_name.get(re.split(r"[+\s(]", expected, maxsplit=1)[0].upper(), set())
            if r["status"] != "ok" or not nodes:
                tally["no call" if r["status"] != "ok" else "truth not in MTree"] += 1
                continue
            got = db.index_of(r["haplogroup"])
            if got in nodes:
                tally["same branch"] += 1
            elif any(db.contains(n, got) for n in nodes):
                tally["deeper"] += 1
            elif any(db.contains(got, n) for n in nodes):
                tally["shallower"] += 1
            else:
                tally["other lineage"] += 1
                misses.append(f"{s.path}: {r['haplogroup']} vs {expected}")
        compared = sum(tally[k] for k in ("same branch", "deeper", "shallower", "other lineage"))
        print(label)
        for key in ("same branch", "deeper", "shallower", "other lineage", "truth not in MTree", "no call"):
            print(f"  {key:<24}{tally[key]}")
        for m in misses:
            print("    " + m)
        summary = f"{compared - tally['other lineage']}/{compared} on the HaploGrep lineage"
        if error:
            print(f"      {summary}")
            continue
        passed = tally["other lineage"] <= 0.01 * compared
        ok = ok and passed
        print(f"{'PASS' if passed else 'FAIL'}  {summary}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
