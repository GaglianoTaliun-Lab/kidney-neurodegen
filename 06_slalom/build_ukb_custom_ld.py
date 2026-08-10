#!/usr/bin/env python3
"""
build_ukb_custom_ld.py -- build a *custom UK Biobank LD panel* for SLALOM's
`--ld-reference custom`, with NO Hail and NO Google Cloud.

For each per-locus .snp file (from build_slalom_snp.R), it:
  1. extracts UKB LD for that locus's variants via PLINK (`--r square`, EUR-unrelated);
  2. aligns the LD matrix to the .snp EFFECT allele (allele2) orientation;
  3. writes a Hail BlockMatrix (pure Python, via hail_bm_io) + an ldcov variant-index
     Parquet (contig, position, ref=allele1, alt=allele2, idx) matched to the .snp;
  4. writes a filtered <base>.snp containing only variants present in the panel (so the
     lead is always in-panel and z<->r orientation is consistent by construction).

SLALOM then runs with EXACT allele matching (no --align-alleles needed):
  slalom --snp <base>.snp --ld-reference custom \
         --custom-ld-path <base>.bm --custom-ld-variant-index-path <base>.vi.parquet \
         --custom-ld-label ukb --dentist-s --abf --summary [--case-control]

Because the variant index stores ref=allele1/alt=allele2 exactly as the .snp, and R is
aligned so its orientation is the effect allele (allele2, the one z=beta/se is expressed
for), SLALOM reads r(lead, j) already sign-consistent with z -> DENTIST-S is correct.

Usage:
  python3 build_ukb_custom_ld.py <manifest.tsv> <ukb_bfile_pattern> <eur_keep> <out_dir>
    <manifest.tsv>       from build_slalom_snp.R (cols incl snp_file, locus_idx, gene,
                         trait, trait_type)
    <ukb_bfile_pattern>  PLINK bfile with %d for chromosome, e.g.
                         /lustre06/.../plink/ukb_chr%d_v2
    <eur_keep>           FID IID list of EUR-unrelated UKB samples (MANDATORY)
    <out_dir>            output dir for <base>.snp/.bm/.vi.parquet + manifest_custom.tsv

Requires: PLINK 1.9 on PATH; numpy, pandas, pyarrow, lz4; hail_bm_io.py alongside this file.
"""
from __future__ import annotations
import os, sys, subprocess, tempfile, shutil
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hail_bm_io import write_hail_blockmatrix

PLINK = os.environ.get("PLINK", "plink")
LD_MAF = float(os.environ.get("LD_MAF", "0.001"))
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _ambiguous(a1: str, a2: str) -> bool:
    return _COMP.get(a1) == a2


def read_snp(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    for c in ("chromosome", "position", "allele1", "allele2", "beta", "se"):
        if c not in df.columns:
            raise ValueError(f"{path}: missing required column {c}")
    df["allele1"] = df["allele1"].astype(str).str.upper()
    df["allele2"] = df["allele2"].astype(str).str.upper()
    df["chromosome"] = df["chromosome"].astype(str)
    return df.reset_index(drop=True)


def extract_plink_ld(chrom, start, end, bfile_pat, eur_keep, tag, workdir):
    """Run PLINK --r square for the window; return (R matrix, bim DataFrame) or (None, None)."""
    bfile = bfile_pat % int(chrom) if "%d" in bfile_pat else bfile_pat
    pref = os.path.join(workdir, f"ld_{tag}")
    cmd = [PLINK, "--bfile", bfile, "--keep", eur_keep, "--chr", str(int(chrom)),
           "--from-bp", str(int(start)), "--to-bp", str(int(end)), "--maf", str(LD_MAF),
           "--keep-allele-order", "--r", "square", "gz", "--make-just-bim", "--out", pref]
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ldf, bmf = pref + ".ld.gz", pref + ".bim"
    if r.returncode != 0 or not (os.path.exists(ldf) and os.path.exists(bmf)):
        return None, None
    bim = pd.read_csv(bmf, sep=r"\s+", header=None,
                      names=["chr", "id", "cm", "pos", "a1", "a2"], dtype={"chr": str})
    bim["a1"] = bim["a1"].astype(str).str.upper()
    bim["a2"] = bim["a2"].astype(str).str.upper()
    R = np.loadtxt(ldf)
    if R.ndim == 1:           # 1-variant window edge case
        R = R.reshape(1, 1)
    return R, bim


def align_r_to_snp(snp: pd.DataFrame, bim: pd.DataFrame, R: np.ndarray):
    """Match .snp variants to the PLINK panel by position; align R to the EFFECT allele.

    Returns (kept_snp, R_eff, keep_mask):
      kept_snp : .snp rows present in the panel (biallelic, unambiguous, allele-set match),
                 in panel order, with an added integer 'idx' (0..m-1)
      R_eff    : m x m LD matrix in the allele2(effect) orientation, aligned to kept_snp
      s_k = +1 if allele2==a1 else -1  ->  R_eff = (s x s^T) * R[panel_rows][:,panel_rows]
    """
    pos_to_bim = {}
    for bi, p in enumerate(bim["pos"].to_numpy()):
        pos_to_bim.setdefault(int(p), bi)          # first occurrence per position
    bim_a1 = bim["a1"].to_numpy(); bim_a2 = bim["a2"].to_numpy()

    keep_snp_rows, bim_rows, signs = [], [], []
    for si in range(len(snp)):
        p = int(snp["position"].iloc[si]); a1 = snp["allele1"].iloc[si]; a2 = snp["allele2"].iloc[si]
        bi = pos_to_bim.get(p)
        if bi is None:
            continue
        A1, A2 = bim_a1[bi], bim_a2[bi]
        if _ambiguous(a1, a2):
            continue
        if {a1, a2} != {A1, A2}:                    # allele set must match
            continue
        # allele2 is the effect allele; align R (plink r is for a1) to allele2 dosage.
        s = 1.0 if a2 == A1 else -1.0
        keep_snp_rows.append(si); bim_rows.append(bi); signs.append(s)

    if len(keep_snp_rows) < 2:
        return snp.iloc[[]].copy(), np.zeros((0, 0)), keep_snp_rows

    br = np.array(bim_rows); s = np.array(signs, dtype=np.float64)
    R_sub = np.asarray(R, dtype=np.float64)[np.ix_(br, br)]
    R_eff = (s[:, None] * s[None, :]) * R_sub
    np.fill_diagonal(R_eff, 1.0)
    kept = snp.iloc[keep_snp_rows].copy().reset_index(drop=True)
    kept["idx"] = np.arange(len(kept), dtype=np.int64)
    return kept, R_eff, keep_snp_rows


def write_variant_index(kept: pd.DataFrame, path: str) -> None:
    """ldcov VariantIndex Parquet: contig, position, ref(=allele1), alt(=allele2), idx."""
    vi = pd.DataFrame({
        "contig": kept["chromosome"].astype(str).to_numpy(),
        "position": kept["position"].astype(np.int64).to_numpy(),
        "ref": kept["allele1"].astype(str).to_numpy(),
        "alt": kept["allele2"].astype(str).to_numpy(),
        "idx": kept["idx"].astype(np.int64).to_numpy(),
    })
    vi.to_parquet(path, index=False)


def build_one(snp_path, chrom, bfile_pat, eur_keep, out_dir, workdir):
    base = os.path.basename(snp_path)[:-4] if snp_path.endswith(".snp") else os.path.basename(snp_path)
    snp = read_snp(snp_path)
    start, end = int(snp["position"].min()), int(snp["position"].max())
    R, bim = extract_plink_ld(chrom, start, end, bfile_pat, eur_keep, base, workdir)
    if R is None:
        return None, "PLINK LD extraction failed"
    kept, R_eff, _ = align_r_to_snp(snp, bim, R)
    if len(kept) < 2:
        return None, f"only {len(kept)} variants matched the UKB panel"
    snp_out = os.path.join(out_dir, base + ".snp")
    bm_out = os.path.join(out_dir, base + ".bm")
    vi_out = os.path.join(out_dir, base + ".vi.parquet")
    kept.drop(columns=["idx"]).to_csv(snp_out, sep=" ", index=False)
    write_hail_blockmatrix(R_eff, bm_out)
    write_variant_index(kept, vi_out)
    return (snp_out, bm_out, vi_out, len(kept)), None


def main():
    if len(sys.argv) < 5:
        sys.exit("usage: build_ukb_custom_ld.py <manifest.tsv> <ukb_bfile_pattern> <eur_keep> <out_dir>")
    manifest, bfile_pat, eur_keep, out_dir = sys.argv[1:5]
    if not os.path.exists(eur_keep):
        sys.exit(f"EUR keep-list not found: {eur_keep} (EUR-unrelated FID IID list is mandatory)")
    os.makedirs(out_dir, exist_ok=True)
    man = pd.read_csv(manifest, sep="\t", dtype=str)
    workdir = tempfile.mkdtemp()
    rows = []
    try:
        for _, m in man.iterrows():
            res, err = build_one(m["snp_file"], m["chr"], bfile_pat, eur_keep, out_dir, workdir)
            tag = f"locus {m.get('locus_idx','?')} {m.get('gene','')} [{m.get('trait','')}]"
            if err:
                print(f"  [skip] {tag}: {err}"); continue
            snp_out, bm_out, vi_out, m_kept = res
            rows.append({"snp_file": os.path.abspath(snp_out), "bm_path": os.path.abspath(bm_out),
                         "vi_path": os.path.abspath(vi_out), "locus_idx": m.get("locus_idx", ""),
                         "chr": m["chr"], "pos": m.get("pos", ""), "gene": m.get("gene", ""),
                         "trait": m.get("trait", ""), "trait_type": m.get("trait_type", "quant"),
                         "n_panel": m_kept})
            print(f"  built {os.path.basename(bm_out)}  ({m_kept} variants)  {tag}")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    if not rows:
        sys.exit("no panels built (check EUR keep-list, bfile pattern, PLINK).")
    out_man = os.path.join(out_dir, "manifest_custom.tsv")
    pd.DataFrame(rows).to_csv(out_man, sep="\t", index=False)
    print(f"\n{len(rows)} UKB LD panels built -> {out_dir}\nmanifest: {out_man}")
    print("Next: run_slalom_custom.sh", out_man, "results_slalom_ukb")


if __name__ == "__main__":
    main()
