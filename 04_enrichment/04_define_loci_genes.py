#!/usr/bin/env python3
"""
fuma_loci_genes.py

Reproduce (locally, no web UI) the FUMA-style LD-based genomic-risk-loci
definition and positional gene mapping used in Yang et al. (Brain Comms 2026),
applied to pleioFDR/conjFDR output.

Protocol reproduced (Yang et al., Supplementary Methods):
  - Independent significant SNPs : conjFDR < 0.05 and LD r2 < 0.6
  - Candidate SNPs               : conjFDR < 0.05 and r2 >= 0.6 with an
                                   independent significant SNP
  - Genomic risk loci            : independent loci merged if within 250 kb
  - Lead SNP                     : lowest-conjFDR SNP in the merged locus
  - Locus boundaries             : min/max position of candidate SNPs
  - Downstream extension         : +/- 500 kb around the lead SNP (reported only)
  - Gene mapping                 : overlap with GENCODE v37 protein-coding genes

IMPLEMENTATION NOTES
  * Candidate SNPs all have conjFDR < 0.05 (Yang's definition), so PLINK clump
    members come from the conjFDR significant file; no genome-wide file needed.
  * The conjFDR file has chr:pos but no rsID. We match each SNP to EUR.bim BY
    POSITION (hg19) to obtain the reference rsID used for clumping. SNPs absent
    from EUR.bim are dropped from clumping (no LD), and the count is reported.
  * Loci are defined by r2<0.6 clumping + 250 kb merge; a separate r2<0.1 lead
    pass is not run because for GENE MAPPING only the boundaries matter, and the
    250 kb merge is functionally equivalent for almost all loci.

Outputs (per pair, in --outdir):
  <pair>_GenomicRiskLoci.tsv   loci with candidate-SNP span, +/-500kb span, counts
  <pair>_genes.txt             unique gene symbols from candidate-SNP boundaries
  <pair>_genes_ext500kb.txt    unique gene symbols from +/-500 kb boundaries

Run on ONE pair first and sanity-check before looping (see SLURM wrapper).
"""

import argparse
import os
import subprocess
import sys
import csv
import tempfile
import shutil

MHC_CHR = 6
MHC_START = 25_000_000   # hg19 broad MHC
MHC_END = 34_000_000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--conjfdr-csv", required=True,
                   help="pleioFDR *_conjfdr_0.05_all.csv (significant SNPs)")
    p.add_argument("--bfile", required=True,
                   help="PLINK reference prefix, e.g. .../data/EUR")
    p.add_argument("--gene-file", required=True,
                   help="protein_coding_genes_grch37.tsv (CHR START END SYMBOL, tab)")
    p.add_argument("--pair", required=True,
                   help="pair label used for output filenames")
    p.add_argument("--outdir", required=True)
    p.add_argument("--plink", default="plink", help="plink 1.9 executable")
    p.add_argument("--clump-r2", type=float, default=0.6,
                   help="r2 for independent significant / candidate SNPs (Yang: 0.6)")
    p.add_argument("--clump-kb", type=int, default=250,
                   help="clump window kb (adjustable; r2>=0.6 partners are local)")
    p.add_argument("--merge-kb", type=int, default=250,
                   help="merge loci within this distance (Yang: 250 kb)")
    p.add_argument("--ext-kb", type=int, default=500,
                   help="downstream +/- extension around lead SNP (Yang: 500 kb)")
    p.add_argument("--gene-pad", type=int, default=10000,
                   help="bp padding around candidate-SNP boundaries for gene mapping "
                        "(matches FUMA positional default of 10 kb)")
    p.add_argument("--sig", type=float, default=0.05, help="conjFDR significance")
    p.add_argument("--exclude-mhc", action="store_true",
                   help="drop chr6:25-34Mb SNPs (recommended for enrichment lists)")
    return p.parse_args()


def detect_conjfdr_col(header):
    """Find the conjFDR value column (conjfdr_<trait>_<trait>), fall back to min_conjfdr."""
    cands = [h for h in header
             if h.startswith("conjfdr_") and h not in ("conjfdr",)]
    if cands:
        return cands[0]
    if "min_conjfdr" in header:
        return "min_conjfdr"
    sys.exit("ERROR: could not find a conjFDR column in the CSV header: %s" % header)


def read_sig_snps(path, sig, exclude_mhc):
    """Return dict (chr,pos)->conjfdr for SNPs with conjFDR < sig, deduped by min conjFDR."""
    with open(path) as f:
        reader = csv.DictReader(f)
        col = detect_conjfdr_col(reader.fieldnames)
        sys.stderr.write("Using conjFDR column: %s\n" % col)
        out = {}
        n_raw = n_mhc = 0
        for row in reader:
            try:
                chrom = int(row["chrnum"]); pos = int(row["chrpos"])
                cf = float(row[col])
            except (ValueError, KeyError):
                continue
            if cf >= sig:
                continue
            n_raw += 1
            if exclude_mhc and chrom == MHC_CHR and MHC_START <= pos <= MHC_END:
                n_mhc += 1
                continue
            key = (chrom, pos)
            if key not in out or cf < out[key]:
                out[key] = cf
    sys.stderr.write("Significant SNPs (conjFDR<%.3g): %d (MHC dropped: %d)\n"
                     % (sig, len(out), n_mhc))
    return out


def map_positions_to_rsid(bim_path, sig_snps):
    """Stream EUR.bim; return rsid->(chr,pos,conjfdr) and (chr,pos)->rsid for matched SNPs."""
    want = set(sig_snps.keys())
    rsid_info = {}
    pos_seen = set()
    with open(bim_path) as f:
        for line in f:
            c, rsid, _gp, bp, _a1, _a2 = line.split()
            key = (int(c), int(bp))
            if key in want and key not in pos_seen:
                pos_seen.add(key)
                rsid_info[rsid] = (key[0], key[1], sig_snps[key])
    matched = len(rsid_info)
    sys.stderr.write("Matched to %s by position: %d / %d (unmatched dropped from clumping)\n"
                     % (os.path.basename(bim_path), matched, len(sig_snps)))
    return rsid_info


def write_assoc(rsid_info, path):
    with open(path, "w") as f:
        f.write("SNP\tP\n")
        for rsid, (_c, _p, cf) in rsid_info.items():
            # PLINK clumps on smallest P; conjFDR already behaves that way.
            f.write("%s\t%g\n" % (rsid, cf))


def run_clump(plink, bfile, assoc, r2, kb, out_prefix):
    cmd = [plink, "--bfile", bfile, "--clump", assoc,
           "--clump-snp-field", "SNP", "--clump-field", "P",
           "--clump-p1", "1", "--clump-p2", "1",  # all SNPs in assoc are already sig
           "--clump-r2", str(r2), "--clump-kb", str(kb),
           "--out", out_prefix]
    sys.stderr.write("Running: %s\n" % " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout + "\n" + res.stderr + "\n")
        sys.exit("ERROR: PLINK clump failed")
    return out_prefix + ".clumped"


def parse_clumped(clumped_path, rsid_info):
    """Return list of clumps; each = dict(chr, members=[rsid,...])."""
    clumps = []
    if not os.path.exists(clumped_path):
        sys.stderr.write("No .clumped file (no clumps formed)\n")
        return clumps
    with open(clumped_path) as f:
        header = f.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.split()
            if not parts:
                continue
            index_snp = parts[idx["SNP"]]
            members = [index_snp]
            sp2 = parts[idx["SP2"]] if "SP2" in idx and len(parts) > idx["SP2"] else "NONE"
            if sp2 not in ("NONE", "."):
                for tok in sp2.split(","):
                    rs = tok.split("(")[0].strip()
                    if rs and rs in rsid_info:
                        members.append(rs)
            members = [m for m in members if m in rsid_info]
            if not members:
                continue
            chrom = rsid_info[index_snp][0]
            clumps.append({"chr": chrom, "members": members})
    return clumps


def clump_span(clump, rsid_info):
    positions = [rsid_info[m][1] for m in clump["members"]]
    cfs = [(rsid_info[m][2], m, rsid_info[m][1]) for m in clump["members"]]
    cfs.sort()  # lowest conjFDR first
    lead_cf, lead_rsid, lead_pos = cfs[0]
    return {"chr": clump["chr"], "start": min(positions), "end": max(positions),
            "lead_rsid": lead_rsid, "lead_pos": lead_pos, "lead_conjfdr": lead_cf,
            "n_members": len(clump["members"])}


def merge_loci(spans, merge_kb):
    """Merge same-chromosome loci whose spans are within merge_kb."""
    merged = []
    by_chr = {}
    for s in spans:
        by_chr.setdefault(s["chr"], []).append(s)
    gap = merge_kb * 1000
    for chrom, group in by_chr.items():
        group.sort(key=lambda x: x["start"])
        cur = dict(group[0])
        cur["n_indep"] = 1
        for s in group[1:]:
            if s["start"] <= cur["end"] + gap:
                cur["end"] = max(cur["end"], s["end"])
                cur["start"] = min(cur["start"], s["start"])
                cur["n_indep"] += 1
                if s["lead_conjfdr"] < cur["lead_conjfdr"]:
                    cur["lead_conjfdr"] = s["lead_conjfdr"]
                    cur["lead_rsid"] = s["lead_rsid"]
                    cur["lead_pos"] = s["lead_pos"]
                cur["n_members"] += s["n_members"]
            else:
                merged.append(cur)
                cur = dict(s); cur["n_indep"] = 1
        merged.append(cur)
    merged.sort(key=lambda x: (x["chr"], x["start"]))
    return merged


def load_genes(path):
    by_chr = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                chrom = int(p[0]); start = int(p[1]); end = int(p[2])
            except ValueError:
                continue
            by_chr.setdefault(chrom, []).append((start, end, p[3]))
    return by_chr


def genes_in(by_chr, chrom, start, end):
    out = []
    for g_start, g_end, sym in by_chr.get(chrom, []):
        if g_start <= end and g_end >= start:   # overlap
            out.append(sym)
    return sorted(set(out))


def main():
    a = parse_args()
    if shutil.which(a.plink) is None:
        sys.exit("ERROR: '%s' not found on PATH.\n"
                 "Load it first:  module load StdEnv/2020 plink/1.9b_6.21" % a.plink)
    os.makedirs(a.outdir, exist_ok=True)
    sig = read_sig_snps(a.conjfdr_csv, a.sig, a.exclude_mhc)
    if not sig:
        sys.exit("No significant SNPs found.")
    rsid_info = map_positions_to_rsid(a.bfile + ".bim", sig)
    if not rsid_info:
        sys.exit("No SNPs matched the reference .bim by position.")

    with tempfile.TemporaryDirectory() as tmp:
        assoc = os.path.join(tmp, "assoc.txt")
        write_assoc(rsid_info, assoc)
        clumped = run_clump(a.plink, a.bfile, assoc, a.clump_r2, a.clump_kb,
                            os.path.join(tmp, "clump"))
        clumps = parse_clumped(clumped, rsid_info)

    spans = [clump_span(c, rsid_info) for c in clumps]
    loci = merge_loci(spans, a.merge_kb)
    sys.stderr.write("Independent clumps: %d  ->  merged loci: %d\n"
                     % (len(spans), len(loci)))

    genes = load_genes(a.gene_file)
    ext = a.ext_kb * 1000
    pad = a.gene_pad

    loci_path = os.path.join(a.outdir, "%s_GenomicRiskLoci.tsv" % a.pair)
    all_genes, all_genes_ext = set(), set()
    with open(loci_path, "w") as out:
        out.write("locus\tchr\tlead_rsid\tlead_pos\tlead_conjfdr\tstart\tend\t"
                  "start_ext500kb\tend_ext500kb\tn_indep\tn_candidate\t"
                  "n_genes\tn_genes_ext\tgenes\n")
        for i, L in enumerate(loci, 1):
            s_ext = max(0, L["lead_pos"] - ext)
            e_ext = L["lead_pos"] + ext
            g = genes_in(genes, L["chr"], max(0, L["start"] - pad), L["end"] + pad)
            g_ext = genes_in(genes, L["chr"], s_ext, e_ext)
            all_genes.update(g); all_genes_ext.update(g_ext)
            out.write("\t".join(str(x) for x in [
                i, L["chr"], L["lead_rsid"], L["lead_pos"], "%.3e" % L["lead_conjfdr"],
                L["start"], L["end"], s_ext, e_ext, L["n_indep"], L["n_members"],
                len(g), len(g_ext), ";".join(g)]) + "\n")

    with open(os.path.join(a.outdir, "%s_genes.txt" % a.pair), "w") as f:
        f.write("\n".join(sorted(all_genes)) + "\n")
    with open(os.path.join(a.outdir, "%s_genes_ext500kb.txt" % a.pair), "w") as f:
        f.write("\n".join(sorted(all_genes_ext)) + "\n")

    sys.stderr.write("Loci: %d | genes (candidate span): %d | genes (+/-500kb): %d\n"
                     % (len(loci), len(all_genes), len(all_genes_ext)))
    sys.stderr.write("Wrote: %s\n" % loci_path)


if __name__ == "__main__":
    main()