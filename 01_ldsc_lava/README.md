# Stage 01 — Global and local genetic correlation (LDSC + LAVA)

This stage estimates genome-wide genetic correlation and SNP heritability (LDSC) and
local genetic correlation across LD blocks (LAVA) for each PD-kidney trait pair.

**These analyses were run with an external Nextflow pipeline, not with code in this
repository.** This directory contains only the project-specific configuration and the
exact command used, so the run can be reproduced. The pipeline source is *not* copied
here — it lives in (and is versioned by) its own repository.

## Pipeline used

- **nf-genetic-correlations** — a Nextflow pipeline for LDSC and LAVA developed by
  Dr. Frida Lona Durazo (lab member / co-author).
  Repository: https://github.com/ape4fld/nf-genetic-correlations

