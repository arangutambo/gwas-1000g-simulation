source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "R", "common.R"))
a <- commandArgs(TRUE); seed <- as.integer(a[1]); m <- as.integer(a[2]); out <- a[3]
pv <- fread(cmd = paste("grep -v '^##'", shQuote(file.path(qc_dir, "chr22_common.pvar"))))
setnames(pv, sub("^#", "", names(pv)))
af <- fread(file.path(qc_dir, "chr22_common.afreq")); setnames(af, sub("^#", "", names(af)))
d <- merge(pv[, .(ID, POS)], af[, .(ID, F = ALT_FREQS)], by = "ID")
d <- d[pmin(F, 1 - F) >= 0.05]
set.seed(seed)
d <- d[sample(.N)]
keep <- integer(0); pos <- numeric(0)
for (i in seq_len(nrow(d))) {                     # greedy: at least 500 kb between causal variants
  if (all(abs(d$POS[i] - pos) >= 5e5)) { keep <- c(keep, i); pos <- c(pos, d$POS[i]) }
  if (length(keep) == m) break
}
if (length(keep) < m) stop(sprintf("only found %d causal variants with 500 kb spacing", length(keep)))
ch <- d[keep][order(POS)]
ch[, beta := rnorm(.N)]                             # effect on standardised genotype
writeLines(ch$ID, file.path(out, "causal_ids.txt"))
fwrite(ch[, .(ID, POS, beta)], file.path(out, "causal.tsv"), sep = "\t")
