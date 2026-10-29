# LD score regression on the simulated scans (Bulik-Sullivan et al. 2015).
# Usage: Rscript R/ldsc.R <work dir> <results dir>
#   E[chi2_j] = 1 + Na + (N h2 / M) * l_j
# The intercept estimates 1 + Na, the confounding term. The simulation knows the
# truth: scenario A has no environmental confounder, scenario B has one.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(TRUE)
WORK <- args[1]          # scratch dir holding ldscores.tsv and the re-run scans
OUT  <- args[2]          # repo results directory
ld <- fread(file.path(WORK, "ldscores.tsv"))
N <- 2504; NBLOCK <- 200

# weighted least squares from sufficient statistics, so the block jackknife is O(1) per block
wls <- function(Sw, Sx, Sy, Sxx, Sxy) {
  den <- Sw * Sxx - Sx^2
  slope <- (Sw * Sxy - Sx * Sy) / den
  c(intercept = (Sy - slope * Sx) / Sw, slope = slope)
}
fit_blocks <- function(x, y, w, blk) {
  s <- function(v) as.numeric(tapply(v, blk, sum))
  Sw <- s(w); Sx <- s(w * x); Sy <- s(w * y); Sxx <- s(w * x * x); Sxy <- s(w * x * y)
  tot <- wls(sum(Sw), sum(Sx), sum(Sy), sum(Sxx), sum(Sxy))
  jk <- vapply(seq_along(Sw), function(b)
    wls(sum(Sw)-Sw[b], sum(Sx)-Sx[b], sum(Sy)-Sy[b], sum(Sxx)-Sxx[b], sum(Sxy)-Sxy[b])[1],
    numeric(1))
  list(fit = tot, jk = jk)
}
ldsc <- function(chi2, l2) {
  keep <- chi2 <= max(80, 0.001 * N)             # LDSC's own outlier rule
  c2 <- chi2[keep]; L <- l2[keep]; M <- length(L)
  blk <- cut(seq_len(M), NBLOCK, labels = FALSE)
  b0 <- fit_blocks(L, c2, rep(1, M), blk)$fit    # step 1: unweighted
  h2_0 <- max(0, (b0[2] * M) / N)
  w <- 1 / (pmax(L, 1) * (1 + N * h2_0 * pmax(L, 1) / M)^2)   # LDSC weights
  r <- fit_blocks(L, c2, w, blk)
  se <- sqrt((NBLOCK - 1) / NBLOCK * sum((r$jk - mean(r$jk))^2))
  c(intercept = unname(r$fit[1]), se = se, h2 = unname((r$fit[2] * M) / N),
    n_snps = M, intercept_unweighted = unname(b0[1]))
}
lam <- function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)

rows <- list()
for (s in c("A", "B")) for (i in 1:10) for (mdl in c("naive", "pc")) {
  d <- fread(file.path(WORK, sprintf("scans/%s_%d_%s.y.glm.linear", s, i, mdl)),
             select = c("ID", "P"))
  d <- merge(d, ld, by = "ID")
  d <- d[is.finite(P) & P > 0]
  setorder(d, POS)                                # jackknife blocks must be contiguous on the genome
  r <- ldsc(qchisq(d$P, 1, lower.tail = FALSE), d$L2)
  rows[[length(rows) + 1]] <- data.table(scenario = s, seed = i, model = mdl,
    lambda = lam(d$P), intercept = r[["intercept"]], se = r[["se"]],
    h2_ldsc = r[["h2"]], n_snps = r[["n_snps"]], intercept_unw = r[["intercept_unweighted"]])
}
res <- rbindlist(rows)
fwrite(res, file.path(OUT, "tables/ldsc.csv"))

summ <- res[, .(lambda_median = round(median(lambda), 2),
                intercept_median = round(median(intercept), 3),
                intercept_q1 = round(quantile(intercept, .25), 3),
                intercept_q3 = round(quantile(intercept, .75), 3),
                intercept_min = round(min(intercept), 3),
                intercept_max = round(max(intercept), 3),
                mean_jackknife_se = round(mean(se), 3),
                n_sig_above_1 = sum(intercept - 2 * se > 1),
                h2_ldsc_median = round(median(h2_ldsc), 3)), by = .(scenario, model)]
write_json(list(per_scan = res, summary = summ, n_replicates = 10L,
  n_snps_median = as.integer(median(res$n_snps)),
  ld_score_mean = round(mean(ld$L2), 2), ld_score_median = round(median(ld$L2), 2),
  ld_score_min = round(min(ld$L2), 2), ld_score_max = round(max(ld$L2), 2),
  n_snps_with_ldscore = nrow(ld)),
  file.path(OUT, "ldsc_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)

cat("\n=== LD score regression, median over 10 replicates ===\n"); print(summ)
cat("\nTruth: A has no environmental confounder (intercept should be ~1);",
    "\n       B has an ancestry-linked shift (intercept should exceed 1).\n")
cat("\n=== unadjusted scans, per replicate ===\n")
print(res[model == "naive", .(scenario, seed, lambda = round(lambda, 1),
                              intercept = round(intercept, 3), se = round(se, 3))])
