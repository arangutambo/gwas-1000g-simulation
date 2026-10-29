# LD scores for the chromosome-22 common set, computed from the same 2,504 samples.
#   l_j = sum_k r2~_jk over SNPs k within +/- 1 Mb, with the standard unbiased
#   adjustment r2~ = r2 - (1 - r2)/(n - 2).
# Usage: Rscript R/ld_scores.R <plink1 bfile prefix> <output tsv>
# Takes about 45 minutes on a laptop. Note that LDSC normally uses an
# ancestry-matched reference panel; here the panel is the same 26-population
# sample being analysed, which inflates r2 between ancestry-informative SNPs.

# LD scores for the chr22 common set, computed from the same 2,504 samples.
# l_j = sum_k r2~_jk over SNPs k within +/- WINDOW bp, with the standard
# unbiased adjustment r2~ = r2 - (1 - r2)/(n - 2).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); bfile <- args[1]; out <- args[2]
WINDOW <- 1e6; STEP <- 1e6

bim <- fread(paste0(bfile, ".bim"), col.names = c("CHR","ID","CM","POS","A1","A2"))
fam <- fread(paste0(bfile, ".fam"))
n <- nrow(fam); m <- nrow(bim); bytes <- ceiling(n / 4)
cat(sprintf("samples %d  variants %d\n", n, m))

con <- file(paste0(bfile, ".bed"), "rb")
hdr <- readBin(con, "raw", 3)
stopifnot(hdr[1] == as.raw(0x6c), hdr[2] == as.raw(0x1b), hdr[3] == as.raw(0x01))  # SNP-major
raw_all <- readBin(con, "raw", bytes * m); close(con)
dim(raw_all) <- c(bytes, m)

# plink .bed 2-bit codes -> dosage; 01 is missing
LUT <- c(2, NA, 1, 0)
unpack <- function(idx) {            # whole window in one rawToBits call
  b <- as.integer(rawToBits(as.vector(raw_all[, idx, drop = FALSE])))
  code <- b[c(TRUE, FALSE)] + 2L * b[c(FALSE, TRUE)]
  dim(code) <- c(bytes * 4L, length(idx))
  g <- LUT[code[seq_len(n), , drop = FALSE] + 1L]
  dim(g) <- c(n, length(idx))
  g
}
standardise <- function(g) {            # mean-impute, then centre and scale
  mu <- colMeans(g, na.rm = TRUE)
  k <- which(is.na(g), arr.ind = TRUE); if (nrow(k)) g[k] <- mu[k[, 2]]
  s <- sqrt(colSums(sweep(g, 2, mu)^2))
  s[s == 0] <- NA                       # monomorphic -> drop from the sum
  sweep(sweep(g, 2, mu), 2, s, "/")
}

l2 <- rep(NA_real_, m)
starts <- seq(min(bim$POS), max(bim$POS), by = STEP)
for (w in seq_along(starts)) {
  s0 <- starts[w]; s1 <- s0 + STEP
  centre <- which(bim$POS >= s0 & bim$POS < s1)
  if (!length(centre)) next
  nb <- which(bim$POS >= s0 - WINDOW & bim$POS < s1 + WINDOW)
  Z <- standardise(unpack(nb))
  ok <- !is.na(Z[1, ])
  ci <- match(centre, nb)
  R <- crossprod(Z[, ci[ok[ci]], drop = FALSE], Z[, ok, drop = FALSE])
  r2 <- R^2
  r2 <- r2 - (1 - r2) / (n - 2)         # unbiased
  l2[centre[ok[ci]]] <- rowSums(r2)
  if (w %% 10 == 0) cat(sprintf("  window %d/%d  pos %.1f Mb\n", w, length(starts), s0/1e6))
}
res <- data.table(ID = bim$ID, POS = bim$POS, L2 = l2)
fwrite(res[!is.na(L2)], out, sep = "\t")
cat(sprintf("wrote %s: %d SNPs with LD scores; mean %.2f  median %.2f  range %.2f-%.2f\n",
            out, sum(!is.na(l2)), mean(l2, na.rm=TRUE), median(l2, na.rm=TRUE),
            min(l2, na.rm=TRUE), max(l2, na.rm=TRUE)))
