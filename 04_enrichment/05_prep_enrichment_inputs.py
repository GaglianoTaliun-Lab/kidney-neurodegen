#!/usr/bin/env python3

import argparse, glob, os, sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci-dir", default="results/fuma_loci")
    p.add_argument("--gene-file", required=True,
                   help="protein_coding_genes_grch37.tsv (symbol in column 4) for background")
    p.add_argument("--outdir", default="results/enrichment")
    p.add_argument("--min-genes", type=int, default=10,
                   help="pairs below this are flagged too_small (not worth submitting)")
    return p.parse_args()


def read_gene_list(path):
    with open(path) as f:
        genes = {line.strip() for line in f if line.strip()}
    return sorted(genes)


def build_background(gene_file):
    syms = set()
    with open(gene_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[3].strip():
                syms.add(parts[3].strip())
    return sorted(syms)


def main():
    a = parse_args()
    inputs_dir = os.path.join(a.outdir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)

    # only candidate-span lists: <pair>_genes.txt, not <pair>_genes_ext500kb.txt
    files = sorted(f for f in glob.glob(os.path.join(a.loci_dir, "*_genes.txt"))
                   if not f.endswith("_ext500kb.txt"))
    if not files:
        sys.exit("No *_genes.txt files found in %s" % a.loci_dir)

    bg = build_background(a.gene_file)
    bg_path = os.path.join(a.outdir, "background_protein_coding.txt")
    with open(bg_path, "w") as out:
        out.write("\n".join(bg) + "\n")

    manifest = os.path.join(a.outdir, "MANIFEST.tsv")
    rows = []
    for f in files:
        pair = os.path.basename(f)[:-len("_genes.txt")]
        genes = read_gene_list(f)
        n = len(genes)
        if n == 0:
            status = "empty"
        elif n < a.min_genes:
            status = "too_small"
        else:
            status = "ready"
        if n > 0:
            with open(os.path.join(inputs_dir, "%s.txt" % pair), "w") as out:
                out.write("\n".join(genes) + "\n")
        rows.append((pair, n, status))

    rows.sort(key=lambda r: (-r[1], r[0]))
    with open(manifest, "w") as out:
        out.write("pair\tn_genes\tstatus\n")
        for pair, n, status in rows:
            out.write("%s\t%d\t%s\n" % (pair, n, status))

    print("Background (protein-coding universe): %d genes -> %s" % (len(bg), bg_path))
    print("\npair\tn_genes\tstatus")
    for pair, n, status in rows:
        print("%s\t%d\t%s" % (pair, n, status))
    ready = [r for r in rows if r[2] == "ready"]
    print("\nReady to submit (>=%d genes): %d pairs -> %s/"
          % (a.min_genes, len(ready), inputs_dir))


if __name__ == "__main__":
    main()
