#!/usr/bin/env python3
"""
Accuracy check on real genotypes: 1,233 1000 Genomes phase 3 males at the 2,366
chrY sites of the Illumina MEGA array (SNAPPY test data, GPL-3.0), compared with
the haplogroups yhaplo called from the full sequence data (Poznik et al. 2016).

Downloads about 1 MB from GitHub on first run (cached in tests/.cache).
"""

import collections
import gzip
import importlib.util
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CACHE = os.path.join(HERE, ".cache")
SNAPPY = "https://raw.githubusercontent.com/chrisgene/snappy/master/test_data/TGP_chrY_MEGA"
TRUTH = "https://raw.githubusercontent.com/23andMe/yhaplo/master/tests/fixtures/output/haplogroups.1000Y.all.txt"

spec = importlib.util.spec_from_file_location("yhp", os.path.join(ROOT, "y-haplogroup-predictor.py"))
yhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(yhp)


def fetch(url):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, url.rsplit("/", 1)[1])
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    return path


def main():
    bim = [line.split() for line in open(fetch(SNAPPY + ".bim"))]
    fam = [line.split()[1] for line in open(fetch(SNAPPY + ".fam"))]
    bed = open(fetch(SNAPPY + ".bed"), "rb").read()
    truth = {}
    for line in open(fetch(TRUTH)):
        p = line.split()
        truth[p[0]] = p[1].split("-", 1)[1]

    per = (len(fam) + 3) // 4
    samples = {sid: yhp.Sample(sid) for sid in fam}
    for j, row in enumerate(bim):
        block = bed[3 + j * per: 3 + (j + 1) * per]
        for i, sid in enumerate(fam):
            code = (block[i // 4] >> (2 * (i % 4))) & 3
            allele = row[4] if code == 0 else row[5] if code == 3 else None
            if allele in ("A", "C", "G", "T"):
                samples[sid].add("Y", int(row[3]), row[1], allele, "call")

    db_path = os.path.join(ROOT, "ydb.tsv.gz")
    names = collections.defaultdict(set)
    with gzip.open(db_path, "rt") as f:
        for line in f:
            if line[0] not in "#N":
                p = line.split("\t")
                for n in p[4].split("/"):
                    names[n.upper()].add(int(p[0]))
    db = yhp.TreeDB(db_path, {int(r[3]) for r in bim})

    class Args:
        error_rate = 0.01
        recurrent_error_rate = 0.1
        min_support = 0.95
        build = "auto"
        complement = "auto"
        exclude = set()

    tally = collections.Counter()
    misses = []
    for sid, s in samples.items():
        r = yhp.analyse(s, db, Args)
        snp = truth.get(sid, "").upper()
        nodes = names.get(snp) or names.get(snp.split(".")[0], set())
        if r["status"] != "ok" or not nodes:
            tally["no call" if r["status"] != "ok" else "truth SNP not in YFull"] += 1
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
            misses.append(f"{sid}: {r['haplogroup']} vs {truth[sid]}")
    compared = sum(tally[k] for k in ("same branch", "deeper", "shallower", "other lineage"))
    for key in ("same branch", "deeper", "shallower", "other lineage", "truth SNP not in YFull", "no call"):
        print(f"{key:<24}{tally[key]}")
    for m in misses:
        print("  " + m)
    ok = tally["other lineage"] <= 0.005 * compared
    print(f"{'PASS' if ok else 'FAIL'}  {compared - tally['other lineage']}/{compared} on the sequence-based lineage")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
