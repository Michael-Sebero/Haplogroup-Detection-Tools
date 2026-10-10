#!/usr/bin/env python3
"""
Build ydb.tsv.gz, the Y-SNP tree database used by y-haplogroup-predictor.py.

Tree:      YFull YTree (https://github.com/YFullTeam/YTree), CC BY 4.0
Positions: Yleaf position tables (https://github.com/genid/Yleaf), GPL-3.0,
           which carry GRCh37/GRCh38 coordinates and ancestral/derived alleles
           compiled from YBrowse, YFull, FTDNA and ISOGG
rsIDs:     Yleaf's chrY dbSNP table (hg19/snp_data.csv)

Needs network access once (or --ytree/--yleaf pointing at local copies).
The predictor itself never touches the network.
"""

import argparse
import collections
import gzip
import io
import json
import os
import sys
import urllib.request

YTREE_URL = "https://raw.githubusercontent.com/YFullTeam/YTree/master/"
YLEAF_URL = "https://raw.githubusercontent.com/genid/Yleaf/master/yleaf/data/"

POSITION_SOURCES = [
    ("yfull14", "{b}/new_positions.txt"),
    ("yfull10", "{b}/yfull_v10_positions_{b}.txt"),
    ("ftdna", "{b}/ftdna_positions_{b}.txt"),
    ("isogg", "{b}/isogg_positions_{b}.txt"),
]

PAR = {
    "hg19": ((10001, 2649520), (59034050, 59363566)),
    "hg38": ((10001, 2781479), (56887903, 57217415)),
}

BASES = set("ACGT")

# Checked against 1000 Genomes phase 3 chrY genotypes (GRCh37): 35 men at all 60,555 sites and
# 1,233 men at the 2,366 sites of the Illumina MEGA array.
# CTS2096: listed alleles disagree with the REF/ALT at the site.
# PF667, PF725, PF803, CTS9296, CTS5457: placed at CT, but every called CT man carries the listed
#   ancestral allele and the A and B men carry the listed derived one (reversed polarity).
# M11839: placed at A-L602, but all 29 called CT men carry the derived allele.
EXCLUDE = {"CTS2096", "PF667", "PF725", "PF803", "CTS9296", "CTS5457", "M11839"}

# Scored like recurrent markers. First two lines: contradict their branch (derived off it or ancestral
# on it) in 3 or more of the 1,233 array-density men, or in 2 or more of the 35 fully genotyped men.
# Rest: reported as unreliable on 23andMe arrays by 23andMe's yhaplo
# (github.com/23andMe/yhaplo, data/variants/isogg.omit.*.23andMe*.txt).
UNSTABLE = {
    "M8963", "L729", "Y4117", "F1449", "YSC0001295", "YP4158", "CTS3316", "L1009", "PAGE52",
    "YSC0000194", "Y512", "L957", "L970", "L977", "CTS3813", "Z3754",
    "M288", "P117", "P80", "P266", "CTS8443", "CTS9188", "CTS9338", "CTS10806", "M5072", "CTS689",
    "CTS3457", "Z4361", "M2706", "M2972", "M3038", "Z5881", "Z4489", "Z262", "S47", "PF7563",
}


def fetch(base, rel, local_dir, cache_dir):
    if local_dir:
        path = os.path.join(local_dir, rel)
        if os.path.exists(path):
            return path
    path = os.path.join(cache_dir, rel.replace("/", "_"))
    if os.path.exists(path):
        return path
    for attempt in range(1, 4):
        print(f"downloading {base}{rel}", file=sys.stderr)
        received = 0
        with urllib.request.urlopen(base + rel, timeout=120) as r, open(path + ".part", "wb") as f:
            expected = int(r.headers.get("Content-Length") or 0)
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
                received += len(chunk)
        if not expected or received == expected:
            os.replace(path + ".part", path)
            return path
        print(f"incomplete download ({received} of {expected} bytes), attempt {attempt} of 3", file=sys.stderr)
    raise SystemExit(f"could not download {base}{rel}")


def load_positions(path):
    table = collections.defaultdict(set)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            anc, der = p[5].upper(), p[6].upper()
            if anc in BASES and der in BASES and anc != der and p[3].isdigit():
                table[p[1]].add((int(p[3]), anc, der))
    return {k: next(iter(v)) for k, v in table.items() if len(v) == 1}


def load_dbsnp(path):
    table = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) >= 4 and p[0].isdigit() and p[1].isdigit():
                table[int(p[1])].append(("rs" + p[0], {p[2].upper(), p[3].upper()}))
    return table


def in_par(pos, build):
    return any(a <= pos <= b for a, b in PAR[build])


def resolve(names, tables):
    votes = collections.Counter()
    for name in names:
        for table in tables:
            hit = table.get(name)
            if hit:
                votes[hit] += 1
    if not votes:
        return None, "unplaced"
    ranked = votes.most_common()
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None, "conflict"
    return ranked[0][0], "ok"


def parse_age(value):
    return int(value) if isinstance(value, int) or (isinstance(value, str) and value.isdigit()) else -1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="ydb.tsv.gz")
    ap.add_argument("--ytree", help="local YTree checkout (current_tree.json, current_version.txt)")
    ap.add_argument("--yleaf", help="local Yleaf data directory (yleaf/data)")
    ap.add_argument("--cache", default=".ydb-cache")
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)

    tree_path = fetch(YTREE_URL, "current_tree.json", args.ytree, args.cache)
    version_path = fetch(YTREE_URL, "current_version.txt", args.ytree, args.cache)
    version = open(version_path).read().strip()

    tables = {}
    for build in ("hg19", "hg38"):
        tables[build] = []
        for _, pattern in POSITION_SOURCES:
            rel = pattern.format(b=build)
            tables[build].append(load_positions(fetch(YLEAF_URL, rel, args.yleaf, args.cache)))
    dbsnp = load_dbsnp(fetch(YLEAF_URL, "hg19/snp_data.csv", args.yleaf, args.cache))

    with open(tree_path) as f:
        root = json.load(f)

    nodes = []
    stack = [(root, -1)]
    while stack:
        node, parent = stack.pop()
        nid = node.get("id") or "ROOT"
        children = node.get("children") or []
        if nid.endswith("*") and not children and not node.get("snps"):
            continue
        idx = len(nodes)
        nodes.append({
            "id": nid,
            "parent": parent,
            "snps": node.get("snps") or "",
            "formed": parse_age(node.get("formed")),
            "tmrca": parse_age(node.get("tmrca")),
        })
        for child in reversed(children):
            stack.append((child, idx))

    stats = collections.Counter()
    markers = []
    for idx, node in enumerate(nodes):
        seen = set()
        for group in node["snps"].split(","):
            names = [n.strip() for n in group.split("/") if n.strip()]
            if not names:
                continue
            stats["groups"] += 1
            if EXCLUDE.intersection(names):
                stats["excluded"] += 1
                continue
            hit19, status = resolve(names, tables["hg19"])
            if hit19 is None:
                stats[status] += 1
                continue
            pos19, anc, der = hit19
            if in_par(pos19, "hg19"):
                stats["par"] += 1
                continue
            hit38, _ = resolve(names, tables["hg38"])
            pos38 = 0
            if hit38 is not None:
                if (hit38[1], hit38[2]) == (anc, der) and not in_par(hit38[0], "hg38"):
                    pos38 = hit38[0]
                else:
                    stats["hg38_mismatch"] += 1
            key = (pos19, der)
            if key in seen:
                stats["duplicate"] += 1
                continue
            seen.add(key)
            rsids = [rs for rs, alleles in dbsnp.get(pos19, ()) if alleles <= {anc, der}]
            bases = [n.split(".")[0] for n in names]
            unstable = bool(UNSTABLE.intersection(names)) or len(bases) != len(set(bases))
            stats["unstable"] += unstable
            markers.append([idx, pos19, pos38, anc, der, "/".join(names), ",".join(rsids), "u" if unstable else ""])
            stats["markers"] += 1

    depth = [0] * len(nodes)
    for i, node in enumerate(nodes):
        if node["parent"] >= 0:
            depth[i] = depth[node["parent"]] + 1

    def is_ancestor(a, b):
        while depth[b] > depth[a]:
            b = nodes[b]["parent"]
        return a == b

    by_pos = collections.defaultdict(list)
    for m in markers:
        by_pos[m[1]].append(m)
    drop = set()
    for pos, group in by_pos.items():
        alleles = {(m[3], m[4]) for m in group}
        ancs = {a for a, _ in alleles}
        if len(ancs) > 1:
            drop.update(id(m) for m in group)
            stats["dropped_polarity_conflict"] += len(group)
            continue
        if len(group) > 1:
            ids = [m[0] for m in group]
            nested = any(is_ancestor(a, b) or is_ancestor(b, a) for i, a in enumerate(ids) for b in ids[i + 1:])
            if nested:
                drop.update(id(m) for m in group)
                stats["dropped_nested_recurrence"] += len(group)
            else:
                stats["recurrent_markers"] += len(group)
    markers = [m for m in markers if id(m) not in drop]
    markers.sort(key=lambda m: (m[1], m[0]))

    buf = io.StringIO()
    buf.write("#ydb\t3\n")
    buf.write(f"#tree\tYFull YTree v{version}\thttps://github.com/YFullTeam/YTree\tCC BY 4.0\n")
    buf.write("#positions\tYleaf position tables (YBrowse, YFull, FTDNA, ISOGG)\thttps://github.com/genid/Yleaf\tGPL-3.0\n")
    buf.write(f"#counts\tnodes={len(nodes)}\tmarkers={len(markers)}\n")
    buf.write("#nodes\tN\tindex\tparent\tid\tformed_ybp\ttmrca_ybp\n")
    for i, node in enumerate(nodes):
        buf.write(f"N\t{i}\t{node['parent']}\t{node['id']}\t{node['formed']}\t{node['tmrca']}\n")
    buf.write("#markers\tnode\tgrch37\tgrch38_minus_grch37\tancestral_derived\tnames\trsids\tflags\n")
    for node_idx, pos19, pos38, anc, der, names, rsids, flags in markers:
        off = str(pos38 - pos19) if pos38 else ""
        buf.write(f"{node_idx}\t{pos19}\t{off}\t{anc}{der}\t{names}\t{rsids}\t{flags}\n")
    with open(args.out, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as f:
        f.write(buf.getvalue().encode("utf-8"))

    for key in sorted(stats):
        print(f"{key}\t{stats[key]}", file=sys.stderr)
    print(f"wrote {args.out}: {len(nodes)} nodes, {len(markers)} markers, YTree v{version}", file=sys.stderr)


if __name__ == "__main__":
    main()
