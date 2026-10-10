#!/usr/bin/env python3
"""
Build mtdb.tsv.gz, the mtDNA tree database used by mt-haplogroup-predictor.py.

Tree:       YFull MTree (https://github.com/YFullTeam/MTree), CC BY 4.0
Reference:  rCRS (NC_012920.1) and the UCSC hg19 chrM sequence (NC_001807.4), from NCBI

The MTree is rooted at the RSRS and lists every mutation in the forward
direction (C16278T on L0a, C16278T! for a back mutation further down). The
builder replays the tree from the root and stores each substitution as the
change between the allele above the branch and the allele below it, so the
predictor can work out the allele expected at any position on any branch.

Needs network access once (or --mtree/--fasta pointing at local copies).
The predictor itself never touches the network.
"""

import argparse
import collections
import gzip
import io
import json
import os
import re
import sys
import urllib.request

MTREE_URL = "https://raw.githubusercontent.com/YFullTeam/MTree/master/"
NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="

REFERENCES = [
    ("rcrs", "NC_012920.1", 16569),
    ("hg19_chrM", "NC_001807.4", 16571),
]

# Never used for placement, following HaploGrep's defaults: 3107 is the placeholder N of the rCRS;
# 16182C and 16183C are by-products of the 16189C length variant in the HVS-I poly-C tract;
# 16519 mutates so often that it says little about the lineage.
EXCLUDE_POSITIONS = {3107, 16182, 16183, 16519}

# Haplogroup of the rCRS itself; replaying the tree down to it must give back the rCRS.
RCRS_HAPLOGROUP = "H2a2a1"

SUB = re.compile(r"^([ACGT]?)(\d+)([ACGT])(!*)$")
DEL = re.compile(r"^([ACGT]+)(\d+)d$")


def fetch(url, name, local_dir, cache_dir):
    if local_dir:
        path = os.path.join(local_dir, name)
        if os.path.exists(path):
            return path
    path = os.path.join(cache_dir, name)
    if os.path.exists(path):
        return path
    for attempt in range(1, 4):
        print(f"downloading {url}", file=sys.stderr)
        received = 0
        with urllib.request.urlopen(url, timeout=120) as r, open(path + ".part", "wb") as f:
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
    raise SystemExit(f"could not download {url}")


def read_fasta(path, length):
    with open(path) as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">")).upper()
    if len(seq) != length:
        raise SystemExit(f"{path}: expected {length} bp, found {len(seq)}")
    return seq


def parse_age(value):
    return int(value) if isinstance(value, int) or (isinstance(value, str) and value.isdigit()) else -1


def tokens(node):
    return [t.strip() for t in node["snps"].split(",") if t.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="mtdb.tsv.gz")
    ap.add_argument("--mtree", help="local MTree checkout (mtree/current_mtree.json, current_version.txt)")
    ap.add_argument("--fasta", help="directory with NC_012920.1.fasta and NC_001807.4.fasta")
    ap.add_argument("--cache", default=".mtdb-cache")
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)

    tree_path = fetch(MTREE_URL + "mtree/current_mtree.json", "current_mtree.json",
                      args.mtree and os.path.join(args.mtree, "mtree"), args.cache)
    version_path = fetch(MTREE_URL + "current_version.txt", "current_version.txt", args.mtree, args.cache)
    version = open(version_path).read().strip()
    number, _, date = version.partition("_")
    seqs = {}
    for key, accession, length in REFERENCES:
        seqs[key] = read_fasta(fetch(NCBI_URL + accession, accession + ".fasta", args.fasta, args.cache), length)
    rcrs = seqs["rcrs"]

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
        # The root lists rCRS-to-RSRS differences flagged with a leading "!"; they are not mutations.
        nodes.append({
            "id": nid,
            "parent": parent,
            "snps": (node.get("snps") or "") if parent >= 0 else "",
            "formed": parse_age(node.get("formed")),
            "tmrca": parse_age(node.get("tmrca")),
        })
        for child in reversed(children):
            stack.append((child, idx))
    children = [[] for _ in nodes]
    for i, node in enumerate(nodes):
        if node["parent"] >= 0:
            children[node["parent"]].append(i)

    def walk(enter, leave):
        stack = [(0, False)]
        while stack:
            i, done = stack.pop()
            if done:
                leave(i)
                continue
            enter(i)
            stack.append((i, True))
            for c in reversed(children[i]):
                stack.append((c, False))

    # Root (RSRS) allele: the starting allele of the first mutation at a position on each lineage.
    first = collections.defaultdict(collections.Counter)
    on_path = collections.Counter()
    added = {}

    def enter_first(i):
        added[i] = []
        for t in tokens(nodes[i]):
            m = SUB.match(t)
            if m and m[1]:
                p = int(m[2])
                if not on_path[p]:
                    first[p][m[1]] += 1
                on_path[p] += 1
                added[i].append(p)

    def leave_first(i):
        for p in added.pop(i):
            on_path[p] -= 1

    walk(enter_first, leave_first)
    # Positions that never mutate on the way to the rCRS keep the rCRS allele at the root, whatever a
    # mislabelled starting allele elsewhere says.
    index = {node["id"]: i for i, node in enumerate(nodes)}
    if RCRS_HAPLOGROUP not in index:
        raise SystemExit(f"{RCRS_HAPLOGROUP} not found in the tree")
    rcrs_path = set()
    i = index[RCRS_HAPLOGROUP]
    while i >= 0:
        rcrs_path.update(int(m[2]) for m in map(SUB.match, tokens(nodes[i])) if m)
        i = nodes[i]["parent"]
    root_state = {}
    for p, votes in first.items():
        allele = votes.most_common(1)[0][0]
        if p in rcrs_path and allele != rcrs[p - 1]:
            root_state[p] = allele

    stats = collections.Counter()
    markers = []
    state = list(rcrs)
    for p, allele in root_state.items():
        state[p - 1] = allele
    undo = {}
    rcrs_check = []

    def enter(i):
        undo[i] = []

        def change(p, frm, to, back):
            markers.append((i, p, frm, to, back))
            undo[i].append((p, state[p - 1]))
            state[p - 1] = to

        for t in tokens(nodes[i]):
            m = SUB.match(t)
            if m:
                p = int(m[2])
                if p in EXCLUDE_POSITIONS:
                    stats["excluded_position"] += 1
                    continue
                frm, to = state[p - 1], m[3]
                if m[1] and frm not in (m[1], "?"):
                    stats["stated_allele_differs"] += 1
                if frm == to:
                    stats["no_change"] += 1
                    continue
                change(p, frm, to, len(m[4]))
                stats["substitutions"] += 1
                continue
            m = DEL.match(t)
            if m:
                start = int(m[2])
                for p in range(start, start + len(m[1])):
                    if p not in EXCLUDE_POSITIONS and state[p - 1] != "?":
                        change(p, state[p - 1], "?", 0)
                        stats["deleted_positions"] += 1
                continue
            stats["insertions_skipped"] += 1
        if nodes[i]["id"] == RCRS_HAPLOGROUP:
            rcrs_check.append("".join(state))

    def leave(i):
        for p, old in reversed(undo.pop(i)):
            state[p - 1] = old

    walk(enter, leave)
    if not rcrs_check:
        raise SystemExit(f"{RCRS_HAPLOGROUP} not found in the tree")
    diff = [p for p in range(1, len(rcrs) + 1)
            if p not in EXCLUDE_POSITIONS and rcrs_check[0][p - 1] not in (rcrs[p - 1], "?")]
    if diff:
        raise SystemExit(f"replaying the tree to {RCRS_HAPLOGROUP} does not give the rCRS at {diff[:20]}")

    markers.sort(key=lambda m: (m[1], m[0]))
    buf = io.StringIO()
    buf.write("#mtdb\t1\n")
    buf.write(f"#tree\tYFull MTree v{number}" + (f" ({date})" if date else "")
              + "\thttps://github.com/YFullTeam/MTree\tCC BY 4.0\n")
    for key, accession, _ in REFERENCES:
        buf.write(f"#{key}\t{seqs[key]}\t{accession}\n")
    buf.write(f"#excluded\t{','.join(map(str, sorted(EXCLUDE_POSITIONS)))}\n")
    buf.write(f"#counts\tnodes={len(nodes)}\tmarkers={len(markers)}\troot_differences={len(root_state)}\n")
    buf.write("#nodes\tN\tindex\tparent\tid\tformed_ybp\ttmrca_ybp\n")
    for i, node in enumerate(nodes):
        buf.write(f"N\t{i}\t{node['parent']}\t{node['id']}\t{node['formed']}\t{node['tmrca']}\n")
    buf.write("#root\tR\tposition\tRSRS allele where it differs from the rCRS\n")
    for p in sorted(root_state):
        buf.write(f"R\t{p}\t{root_state[p]}\n")
    buf.write("#markers\tnode\tposition\tbefore_after (? = deleted)\tback_mutation\n")
    for node, p, frm, to, back in markers:
        buf.write(f"{node}\t{p}\t{frm}{to}\t{'!' * back}\n")
    with open(args.out, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as f:
        f.write(buf.getvalue().encode("utf-8"))

    for key in sorted(stats):
        print(f"{key}\t{stats[key]}", file=sys.stderr)
    print(f"wrote {args.out}: {len(nodes)} nodes, {len(markers)} markers, {len(root_state)} RSRS differences, "
          f"MTree v{version}", file=sys.stderr)


if __name__ == "__main__":
    main()
