# Self-test: REML should recover a known heritability from a simulated GRM.
root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..")))
source(file.path(root, "R", "reml.R"))
set.seed(42)
n <- 400; m <- 2000; h2 <- 0.5
est <- replicate(15, {
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  Z <- scale(G); A <- tcrossprod(Z) / m
  y <- Z %*% rnorm(m, sd = sqrt(h2 / m)) + rnorm(n, sd = sqrt(1 - h2))
  e <- eigen(A, symmetric = TRUE)
  reml_h2(as.numeric(y), list(U = e$vectors, lam = e$values))[["h2"]]
})
cat(sprintf("mean REML h2 = %.3f (truth %.2f), sd %.3f\n", mean(est), h2, sd(est)))
stopifnot(abs(mean(est) - h2) < 0.08)
cat("OK\n")
