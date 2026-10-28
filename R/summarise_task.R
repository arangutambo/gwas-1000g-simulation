source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "R", "common.R"))
source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "R", "reml.R"))
a <- commandArgs(TRUE); scen <- a[1]; seed <- as.integer(a[2]); out <- a[3]
ch <- fread(file.path(out, "causal.tsv")); truth <- fread(file.path(out, "truth.csv"))
ALPHA <- 5e-8; WIN <- 5e5

score <- function(path, label) {
  d <- fread(path); setnames(d, sub("^#", "", names(d)))
  d <- d[TEST == "ADD" & !is.na(P)]
  sig <- d[P < ALPHA]
  near <- function(pos) any(abs(pos - ch$POS) <= WIN)
  detected <- sum(vapply(ch$POS, function(p) any(abs(sig$POS - p) <= WIN), logical(1)))
  far <- sig[!vapply(POS, near, logical(1))]
  # thinned scan for plots: every P < 0.01 plus 5% of the rest
  set.seed(1); keep <- d$P < 0.01 | runif(nrow(d)) < 0.05
  saveRDS(d[keep, .(POS, P)], file.path(out, paste0("scan_", label, ".rds")))
  setNames(c(lambda_gc(d$P), nrow(sig), detected, length(unique(floor(far$POS / 1e6)))),
           paste0(c("lambda_", "n_sig_", "detected_", "false_loci_"), label))
}
r_naive <- score(file.path(out, "naive.y.glm.linear"), "naive")
r_pc    <- score(file.path(out, "pcadj.y.glm.linear"), "pc")

eig <- readRDS(file.path(qc_dir, "grm_eigen.rds"))
ph <- fread(file.path(out, "pheno.tsv")); y <- ph$y[match(eig$iid, ph[["#IID"]])]
pcs <- fread(file.path(qc_dir, "pcs.eigenvec")); X <- as.matrix(pcs[match(eig$iid, pcs[["#IID"]]), paste0("PC", 1:10), with = FALSE])
h_naive <- reml_h2(y, eig); h_pc <- reml_h2(y, eig, X)

res <- data.table(scenario = scen, seed = seed, n_causal = nrow(ch), true_h2 = truth$true_h2,
                  t(r_naive), t(r_pc), h2_naive = h_naive[["h2"]], h2_pc = h_pc[["h2"]],
                  reml_ok = h_naive[["converged"]] == 1 && h_pc[["converged"]] == 1)
fwrite(res, file.path(out, "summary.csv"))
unlink(file.path(out, c("naive.y.glm.linear", "pcadj.y.glm.linear")))   # large; the thinned scans are kept
print(res)
