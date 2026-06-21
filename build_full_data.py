#!/usr/bin/env python3
"""Build the full 3-species sequence store for the alignment tool.

One longest isoform per gene symbol, per species:
  - X. laevis : Xenbase v10.1 proteome FASTA (gene = 2nd '|' field; .L/.S kept
                separate so homeologs are distinct).
  - human/mouse: UniProt reviewed (Swiss-Prot) FASTA (gene = GN= field).

Output: sequences.tsv  (species, gene, accession, length, sequence)
matching the schema the app already reads.

Inputs are expected at:
  /tmp/XENLA_10.1.pep.fa     (Xenbase: download.xenbase.org/.../XENLA_10.1_Xenbase.pep.fa.gz)
  /tmp/hs_proteome.fasta     (UniProt: organism_id:9606  AND reviewed:true, fasta)
  /tmp/mm_proteome.fasta     (UniProt: organism_id:10090 AND reviewed:true, fasta)
"""
import os, re, csv
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences.tsv")

XL = "/tmp/XENLA_10.1.pep.fa"
HS = "/tmp/hs_proteome.fasta"
MM = "/tmp/mm_proteome.fasta"


def parse_fasta(path):
    recs, hdr, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(buf)))
                hdr, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if hdr is not None:
            recs.append((hdr, "".join(buf)))
    return recs


def collect(path, species, gene_fn, acc_fn):
    """Keep the longest sequence per gene symbol."""
    best = {}
    for hdr, seq in parse_fasta(path):
        gene = gene_fn(hdr)
        if not gene:
            continue
        acc = acc_fn(hdr)
        if gene not in best or len(seq) > len(best[gene][1]):
            best[gene] = (acc, seq)
    return [(species, g, acc, len(s), s) for g, (acc, s) in best.items()]


# X. laevis: header 'RefSeq:ACC|gene|XBmRNA...|...'
def xl_gene(h):
    p = h.split("|")
    return p[1] if len(p) >= 2 else None
def xl_acc(h):
    return h.split("|", 1)[0].replace("RefSeq:", "")

# UniProt: 'sp|ACC|NAME ... GN=GENE ...'
def up_gene(h):
    m = re.search(r"\bGN=(\S+)", h)
    return m.group(1) if m else None
def up_acc(h):
    return h.split("|")[1] if "|" in h else h.split()[0]


rows = []
rows += collect(XL, "X. laevis", xl_gene, xl_acc)
rows += collect(HS, "human", up_gene, up_acc)
rows += collect(MM, "mouse", up_gene, up_acc)

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["species", "gene", "accession", "length", "sequence"])
    for sp, g, acc, n, s in rows:
        w.writerow([sp, g, acc, n, s])

by_sp = defaultdict(int)
for sp, *_ in rows:
    by_sp[sp] += 1
print(f"wrote {OUT}: {len(rows):,} sequences")
for sp in ("X. laevis", "human", "mouse"):
    print(f"  {sp}: {by_sp[sp]:,}")
print(f"file size: {os.path.getsize(OUT)/1e6:.1f} MB")
