#!/usr/bin/env python3
"""Recover X. laevis .L/.S homeologs that have a transcript model but no protein
in the Xenbase v10.1 proteome, by translating their longest ORF.

Most transcript-only homeologs are RefSeq-predicted ncRNA (XR_/NR_) whose ORFs
are spurious. We keep only the ones that are clearly real protein-coding genes
mis-classified as ncRNA: the translated ORF must be within 0.6-1.4x the length
of the partner homeolog's annotated protein (e.g. dsp.S 2865aa vs dsp.L 2866aa).

Inputs:
  /tmp/XENLA_10.1.pep.fa     (proteome; symbol -> existing protein length)
  /tmp/xl_transcripts.fa.gz  (Xenbase v10.1 transcripts, gzipped)

Outputs:
  homeolog_orfs.tsv          (the validated translated homeologs)
  + appends them to sequences.tsv (idempotent: drops prior "(ORF)" rows first)
"""
import os, re, gzip, csv

HERE = os.path.dirname(os.path.abspath(__file__))
PEP = "/tmp/XENLA_10.1.pep.fa"
TX  = "/tmp/xl_transcripts.fa.gz"
OUT = os.path.join(HERE, "homeolog_orfs.tsv")

CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}


def translate(nt):
    return ''.join(CODON.get(nt[i:i+3], 'X') for i in range(0, len(nt) - 2, 3))


def longest_orf(nt):
    """Longest M..stop peptide across the 3 forward frames."""
    nt = nt.upper().replace('U', 'T')
    best = ''
    for f in range(3):
        aa = translate(nt[f:])
        for seg in aa.split('*'):
            m = seg.find('M')
            if m >= 0 and len(seg) - m > len(best):
                best = seg[m:]
    return best


def pep_lengths(path):
    """symbol -> longest protein length in the proteome."""
    out, sym, ln = {}, None, 0
    def flush():
        if sym:
            out[sym] = max(out.get(sym, 0), ln)
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
                p = line[1:].split('|'); sym = p[1] if len(p) >= 2 and p[1] else None; ln = 0
            else:
                ln += len(line.strip())
    flush()
    return out


def iter_fasta_gz(path):
    hdr, buf = None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:
                    yield hdr, ''.join(buf)
                hdr, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if hdr is not None:
            yield hdr, ''.join(buf)


prot_len = pep_lengths(PEP)
have_protein = set(prot_len)
print(f"proteome symbols: {len(have_protein):,}")

# Best (longest-ORF) protein per transcript-only .L/.S homeolog symbol.
best = {}   # symbol -> (accession, protein)
scanned = 0
for hdr, seq in iter_fasta_gz(TX):
    p = hdr.split('|')
    if len(p) < 2 or not p[1]:
        continue
    sym = p[1]
    if not re.search(r'\.[LS]$', sym) or sym in have_protein:
        continue
    scanned += 1
    prot = longest_orf(seq)
    if not prot:
        continue
    acc = p[0].replace('RefSeq:', '') + " (ORF)"
    if sym not in best or len(prot) > len(best[sym][1]):
        best[sym] = (acc, prot)

# Keep only those validated against the partner homeolog's protein length.
def partner(sym):
    base = re.sub(r'\.[LS]$', '', sym)
    return base + ('.L' if sym.endswith('.S') else '.S')

validated = {}
for sym, (acc, prot) in best.items():
    pl = prot_len.get(partner(sym))
    if pl and 0.6 <= len(prot) / pl <= 1.4:
        validated[sym] = (acc, prot)

with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(["species", "gene", "accession", "length", "sequence"])
    for sym in sorted(validated, key=lambda s: -len(validated[s][1])):
        acc, prot = validated[sym]
        w.writerow(["X. laevis", sym, acc, len(prot), prot])

print(f"transcript-only homeolog symbols scanned: {scanned}")
print(f"candidate ORFs: {len(best)} | validated against partner protein: {len(validated)}")
print(f"-> {OUT}")

# Append to sequences.tsv (idempotent: drop any prior "(ORF)" rows first).
SEQ = os.path.join(HERE, "sequences.tsv")
if os.path.exists(SEQ):
    with open(SEQ) as fh:
        kept = [ln for ln in fh if "(ORF)" not in ln]
    with open(SEQ, "w") as fh:
        fh.writelines(kept)
        for sym in sorted(validated, key=lambda s: -len(validated[s][1])):
            acc, prot = validated[sym]
            fh.write(f"X. laevis\t{sym}\t{acc}\t{len(prot)}\t{prot}\n")
    print(f"appended {len(validated)} validated homeolog proteins to {SEQ}")
