#!/usr/bin/env bash
# Download 1000 Genomes Phase 3 chr22 genotypes and the sample panel.
# Public data, no registration. chr22 only (~200 MB) keeps this laptop-sized.
set -euo pipefail
cd "$(dirname "$0")/../data"
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
curl -fLO -C - "$BASE/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
curl -fLO -C - "$BASE/integrated_call_samples_v3.20130502.ALL.panel"
ls -lh
