 <p align="center">
  <img
    src="https://upload.wikimedia.org/wikipedia/commons/c/ca/World_Map_of_Y-DNA_Haplogroups.png"
    alt="World Map of Y-DNA Haplogroups"
    width="800"
  />
</p>

 <p align="center">
  <img
    src="https://upload.wikimedia.org/wikipedia/commons/0/04/Human_migrations_and_mitochondrial_haplogroups.PNG?utm_source=commons.wikimedia.org&utm_campaign=index&utm_content=original"
    alt="World Map of Y-DNA Haplogroups"
    width="800"
  />
</p>

## How to Install

```
git clone https://github.com/Michael-Sebero/Haplogroup-Detection-Tools

python3 /path/to/Haplogroup-Detection-Tools/y-haplogroup-predictor.py

python3 /path/to/Haplogroup-Detection-Tools/mt-haplogroup-predictor.py
```

Requires Python 3.7 or newer. No other packages are needed. Keep `ydb.tsv.gz` and `mtdb.tsv.gz` in the same folder as the scripts.

## Description

Haplogroup-Detection-Tools reads the raw DNA file from a consumer test and places the sample on the YFull trees of human lineages:

- **Y-Haplogroup-Predictor** (`y-haplogroup-predictor.py`) places the paternal line on the YFull Y-chromosome tree, which has 43,758 branches covering haplogroups A00 through T.
- **MT-Haplogroup-Predictor** (`mt-haplogroup-predictor.py`) places the maternal line on the YFull mtDNA tree (MTree), which has 27,988 branches covering haplogroups L0 through Z. It also reads full mitochondrial sequences.

Both report the deepest branch your data supports, the mutations behind each step and any calls that disagree with the result.

The intention of this software is to provide a FOSS alternative to detecting haplogroups. Data breaches occur often and institutions shouldn't be trusted with processing or storing client genetic data over a long period of time. Haplogroup-Detection-Tools is completely private and runs offline: the tree databases ship with the scripts and your DNA never leaves your machine.

## How to use Haplogroup-Detection-Tools

- Acquire a copy of your raw DNA data from a genetic testing platform. This is a plain text file (`.txt`, `.csv` or `.vcf`) with one SNP per line.
- Run either script with the file path, or run it without arguments and enter the path when prompted. Y haplogroups need a male sample; everyone has an mtDNA haplogroup.

```
python3 y-haplogroup-predictor.py AncestryDNA.txt

python3 mt-haplogroup-predictor.py AncestryDNA.txt
```

Several files can be given at once:

```
python3 y-haplogroup-predictor.py genome_father.txt genome_son.txt
```

If the vendor delivered the file compressed (`.zip`, `.gz` or `.bz2`), it can be passed as it is without unpacking it first.

### Supported files

| Vendor | Y row | mtDNA row |
|---|---|---|
| 23andMe | `rs9786153  Y  22739367  C` | `<id>  MT  16311  C` |
| AncestryDNA | `rs9786153  24  22739367  C  C` | `<id>  26  16311  C  C` |
| FamilyTreeDNA, MyHeritage | `"rs9786153","Y","22739367","CC"` | `"<id>","MT","16311","CC"` |
| LivingDNA, tellmeGen, Genes for Good | same layout as 23andMe | same layout as 23andMe |
| Illumina GSA reports | `[Header]` / `[Data]` sections | `[Header]` / `[Data]` sections |
| VCF (Nebula, Dante Labs, other whole-genome services) | single- or multi-sample | single- or multi-sample, gVCF |

The mtDNA tool also reads:

| File | Example |
|---|---|
| FASTA (FamilyTreeDNA mtFull Sequence, GenBank) | complete or partial sequence |
| HaploGrep `.hsd` files and mutation lists | `73G 263G 2706G 16311C` or `A73G C2706T` |

Tab-, comma- and space-separated text files are all accepted, with or without a header line; the columns are recognised from the header or the row layout.

## Y-Haplogroup-Predictor

The reference build is detected from the data: GRCh37 or GRCh38 positions, with an rsID fallback for older build 36 files.

### Options

| Option | Effect |
|---|---|
| `--summary` | one tab-separated line per file, for batches |
| `--json` | full result as JSON |
| `-v` | list every branch on the path and every conflicting call |
| `--snp L21,M269,rs9786153` | where these SNPs sit on the tree and each file's call at them |
| `--node R-L21` | a branch's path, YFull age estimate and subclades |
| `--build GRCh37` | override build detection (`GRCh38`, `rsID`) |
| `--sample NAME` | sample to read from a multi-sample VCF |
| `--exclude L150,P25` | ignore these SNPs |
| `--error-rate 0.01` | per-call genotyping error rate assumed by the model |
| `--min-support 0.95` | support needed to report a branch |

### Example

1000 Genomes sample HG00096, limited to the 2,366 Y positions of the Illumina MEGA array and saved as a 23andMe file:

```
genome_HG00096.txt
  Source        23andMe, markers matched by GRCh37 positions
  Y calls       2,348 of 2,366 rows called, 2 heterozygous; 2,216 tree SNPs informative
  Haplogroup    R-S733  (major clade R1b)
  Support       99.00%
  YFull age     formed 2,700 ybp, TMRCA 1,750 ybp

  branch         derived ancestral   support  SNPs
  root..P-P226       176         0   >99.99%  20 branches above R
  R                   21         0   >99.99%  PF6013, Y472, PF6016, P229
  R1                  19         0   >99.99%  P242, P238, PF6116, F132
  R1b                  7         0   >99.99%  M343, PF6246, L822, CTS2134
  ...
  R-M269              48         0   >99.99%  S10, PF6400, PF6419, YSC0001293
  R-L151               3         0   >99.99%  L151, YSC0000191, L11
  R-P312               1         0   >99.99%  P312
  R-L21                1         0   >99.99%  L21
  R-DF13               1         0   >99.99%  CTS8221
  R-L1335              1         0    99.99%  L1335
  R-S733               1         0    99.00%  CTS6838

  Subclades of R-S733
    3 subclades without markers in this file
```

"Derived" means the sample carries the mutation that defines the branch; "ancestral" means it doesn't. YFull ages are in years before present.

### How it works

Each Y SNP in your file that matches a tree marker is scored as derived or ancestral. For every branch of the tree, the script checks how well your calls fit the assumption that you belong to that branch: markers on the path from the root to it should be derived and all others ancestral. Each call that doesn't fit counts against the branch, weighted by how likely a genotyping error is. The script then reports the deepest branch whose support is at least 95%. Support compares the best placement inside a branch with the best placement outside it.

A branch confirmed by a single SNP gets 99% support. Heterozygous Y calls, no-calls and indels are ignored. SNPs known to occur on several branches, or to misbehave on genotyping arrays, count for less. Files with almost no Y calls are reported as female samples.

### Accuracy

Checked against haplogroups called from full sequence data for 1000 Genomes men (Poznik et al. 2016):

- **Array density** (1,233 men, the 2,366 Y positions of the Illumina MEGA array): 1,213 of the 1,214 comparable men were placed on the correct lineage. The remaining one is on the matching O-M95 lineage, but its reference SNP only exists under haplogroup J in YFull.
- **Full sequence** (35 men, all 60,555 Y positions): all 34 comparable men were placed on the correct lineage.
- **Genotyping errors** (300 of the array-density samples with random errors added): 300/300 on the correct lineage at 1% errors, 299/300 at 3% and 300/300 at 5%.

`python3 tests/run_tests_y.py` runs the offline checks. `python3 tests/validate_1000g_y.py` repeats the array comparison; it downloads about 1 MB of 1000 Genomes data from GitHub.

### Data sources

`ydb.tsv.gz` is built by `build-ydb.py` from:

- **YFull YTree v14.06.0** (22 September 2026), https://github.com/YFullTeam/YTree, licensed CC BY 4.0. It is converted to tab-separated text with positions added; empty placeholder branches (`R1b*` and similar) are dropped.
- **Yleaf position tables** (Erasmus MC), https://github.com/genid/Yleaf, GPL-3.0. They provide the GRCh37 and GRCh38 positions and alleles, compiled from YBrowse, YFull, FTDNA and ISOGG, plus the rsIDs.
- **Quality control:**
  - Seven SNPs whose alleles or placement contradict 1000 Genomes data are removed.
  - 422 SNPs are down-weighted: multi-copy SNPs, SNPs that contradict their branch in 1000 Genomes men and SNPs that 23andMe's yhaplo reports as unreliable on its arrays.
  - The lists are at the top of `build-ydb.py`.

To update the tree after a new YFull release (a one-time download of 122 MB):

```
python3 build-ydb.py
```

## MT-Haplogroup-Predictor

mtDNA rows are found under `MT`, `M`, `chrM` or `26`.

- **Positions** are read on the rCRS, the numbering every consumer test uses. VCFs aligned to UCSC hg19, whose `chrM` is a different sequence (NC_001807, 16,571 bp), are recognised from the data and converted.
- **VCFs and mutation lists** that only list differences from the reference are filled in with the reference allele everywhere else.
- **Lists written against the RSRS**, as FamilyTreeDNA can show them, are recognised from their alleles.
- **FASTA files** are aligned to the rCRS, so insertions and deletions in the sequence are handled.

### Options

| Option | Effect |
|---|---|
| `--summary` | one tab-separated line per file, for batches |
| `--json` | full result as JSON |
| `-v` | list every branch on the path and every conflicting call |
| `--snp 16519,T152C` | which branches mutate at these positions and each file's call there |
| `--node H1c3` | a haplogroup's path, YFull age, mutations, differences from the rCRS and subclades |
| `--coordinates hg19` | override detection of the position numbering (`rCRS`) |
| `--unlisted unknown` | do not fill positions a VCF or mutation list leaves out (`reference`) |
| `--sample NAME` | sample to read from a multi-sample VCF, FASTA or `.hsd` file |
| `--exclude 16189,152` | ignore these positions |
| `--error-rate 0.01` | per-call error rate at positions that mutate once in the tree |
| `--recurrence-rate 0.008` | added to the error rate for each further time a position mutates in the tree |
| `--min-support 0.95` | support needed to report a branch |

### Example

1000 Genomes sample NA12878, limited to the 326 mtDNA positions of the Illumina Omni 2.5 array and saved as a 23andMe file:

```
genome_NA12878.txt
  Source        23andMe, positions on the rCRS
  MT calls      326 of 326 rows called; 326 tree positions informative
  Haplogroup    H13a  (major clade H)
  Support       99.00%
  YFull age     formed 9,700 ybp, TMRCA 7,000 ybp

  branch         found  missing   support  mutations
  L..HV             10        0   >99.99%  12 branches above H
  H                  1        0    99.96%  T7028C
  H13                1        0    99.96%  C14872T
  H13a               1        0    99.00%  C2259T

  Subclades of H13a
    H13a1       not in the file; downstream mutations missing: T3396C, A15758G, G14323A, C16167T, T11025C, A12397G
    H13a2       not in the file; downstream mutations missing: A15218G, G1438A, A15662G, T10463C
```

The same sample's whole mtDNA genome, as a VCF listing its 19 differences from the rCRS, is placed four branches further down at H13a1a1a.

"Found" means the sample carries the mutation that defines the branch; "missing" means it still has the allele from before that mutation. Mutations are written the YFull way, from the allele above the branch to the allele below it; a trailing `!` marks a back mutation. YFull ages are in years before present.

### How it works

The MTree is rooted at the RSRS, the reconstructed sequence of the most recent common maternal ancestor, and each branch lists the mutations that separate it from the branch above, back mutations included. Replaying them from the root gives the allele expected at every position on every branch.

Each mtDNA position in your file that the tree uses is compared with that expected allele. For every branch of the tree, the script checks how well your calls fit the assumption that you belong to that branch, and each call that doesn't fit counts against it, weighted by how likely a mismatch is at that position. The script then reports the deepest branch whose support is at least 95%. Support compares the best placement inside a branch with the best placement outside it.

Many mtDNA positions mutate again and again in different lineages: 16311 does so 401 times in the MTree, 152 does 233 times. A call at such a position says little about the lineage, so the allowed error rate rises with each further time a position mutates in the tree. The increase is 0.002 per occurrence for whole sequences (files with calls at 90% or more of the positions the tree uses), where every branch on the way is tested, and 0.008 for arrays and partial sequences, where most of the tree goes unseen and one parallel mutation could otherwise carry a placement on its own. A branch confirmed by a single mutation at a position that occurs only once in the tree gets 99% support.

The following are left out:

- Heteroplasmic calls (two alleles at one position), no-calls, insertions and deletions.
- Positions 3107, 16182, 16183 and 16519, following HaploGrep's defaults.

### Accuracy

Checked against the haplogroups HaploGrep 2 called on PhyloTree 17 for the 2,504 1000 Genomes phase 3 samples, from the same mtDNA calls (haplocheck study, Weissensteiner et al. 2021). 19 samples carry a PhyloTree name that does not exist in the MTree and are left out.

- **Whole mtDNA genome** (2,485 samples): 2,469 placed on the HaploGrep lineage. 1,025 are on the same branch, 1,212 are deeper because the MTree has five times as many branches as PhyloTree 17, and 232 stop above it. The 16 others sit where the two trees disagree. 13 of them fit the MTree branch they were given with no missing mutation and the other 3 miss one, but PhyloTree 17 places the lineage elsewhere, for example R11b1 against the MTree's R11b2a1.
- **Array density** (2,485 samples, the 326 mtDNA positions of the Illumina Omni 2.5 array): 2,484 placed on the HaploGrep lineage. Most stop a few branches above it, because the array has no positions for them.
- **Genotyping errors** (the array-density samples with random errors added): 2,476 of 2,485 on the HaploGrep lineage at 1% errors, 2,469 at 3% and 2,448 at 5%.

`python3 tests/run_tests_mt.py` runs the offline checks. `python3 tests/validate_1000g_mt.py` repeats the 1000 Genomes comparison (`--errors 0.01,0.03,0.05` adds the error runs); it downloads about 1 MB of data from GitHub.

### Data sources

`mtdb.tsv.gz` is built by `build-mtdb.py` from:

- **YFull MTree v1.02.24043** (22 September 2026), https://github.com/YFullTeam/MTree, licensed CC BY 4.0. It is converted to tab-separated text, with each mutation stored as the change from the allele above the branch to the allele below it; empty placeholder branches (`H1*` and similar) are dropped.
- **rCRS** (NC_012920.1) and the **UCSC hg19 chrM** sequence (NC_001807.4), from NCBI GenBank.
- **Quality control:** the RSRS root is reconstructed from the starting alleles of the first mutation on each lineage, and the build stops unless replaying the tree down to H2a2a1, the rCRS's own haplogroup, reproduces the rCRS exactly.

To update the tree after a new YFull release (a one-time download of 5 MB):

```
python3 build-mtdb.py
```

## Limitations

- Consumer arrays carry a few thousand Y SNPs (a 2011 23andMe file has 3,488 Y rows) and a few hundred to a few thousand mtDNA positions. Results usually stop several branches above what sequencing (Big Y, mtFull Sequence, YSEQ, whole genome) would reach.
- Y branch names follow YFull (R-L21, I-Z59). ISOGG long-form names such as R1b1a1b1a1a2c1 are not produced.
- mtDNA branch names follow the YFull MTree, which keeps PhyloTree names for established haplogroups and adds its own for newer branches (H16-a, L2a1-a2h5). HaploGrep, which uses PhyloTree 17, gives different names to some branches.
- mtDNA insertions and deletions are not used, including the 9-bp deletion that defines haplogroup B; B is found from its other mutations.
- Y-STRs are not analysed.

*This repository has been repurposed; its existing stars are from PrivateGPT4Linux-CLI*
