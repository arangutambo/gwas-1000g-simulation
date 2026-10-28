#!/usr/bin/env bash
# Stage 1: variant/sample QC, LD pruning, PCA and GRM from the 1000 Genomes chr22 VCF.
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p "$QC"
cd "$PROJECT_ROOT"

# VCF -> pgen. Biallelic SNPs only. MAF >= 1% and <= 2% missing.
# No HWE filter on purpose: the sample mixes 26 populations, so Hardy-Weinberg
# departures (the Wahlund effect) would flag true variants, not genotyping errors.
plink2 --vcf "$VCF" \
  --max-alleles 2 --snps-only just-acgt --rm-dup exclude-mismatch \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 20 truncate \
  --maf 0.01 --geno 0.02 \
  --make-pgen --out "$QC/chr22_qc"

# Common variants (MAF >= 5%) used for the association scans.
plink2 --pfile "$QC/chr22_qc" --maf 0.05 --make-pgen --freq --out "$QC/chr22_common"

# LD-pruned set for PCA: window 200 SNPs, step 50, r^2 < 0.2.
plink2 --pfile "$QC/chr22_common" --indep-pairwise 200 50 0.2 --out "$QC/ld"
plink2 --pfile "$QC/chr22_common" --extract "$QC/ld.prune.in" --pca 10 --out "$QC/pcs"

# GRM from all QC'd variants (square text matrix + ID file).
plink2 --pfile "$QC/chr22_qc" --make-rel square --out "$QC/grm"

# Sample table (IID, population, super-population) from the 1000G panel. POSIX awk only.
awk 'BEGIN{OFS="\t"; print "IID","pop","super_pop"} NR>1{print $1,$2,$3}' "$PANEL" > "$QC/samples.tsv"

# Cache the GRM eigendecomposition so each task can reuse it.
Rscript "$PROJECT_ROOT/R/prep_grm.R"
echo "QC done: $(wc -l < "$QC/chr22_common.pvar") common variant lines, $(( $(wc -l < "$QC/samples.tsv") - 1 )) samples"
