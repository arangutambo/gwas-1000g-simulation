source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "R", "common.R"))
a <- commandArgs(TRUE); scen <- a[1]; seed <- as.integer(a[2]); out <- a[3]
ch <- fread(file.path(out, "causal.tsv"))
raw <- fread(file.path(out, "causal.raw"))
iid <- as.character(raw[[grep("^#?IID$", names(raw))[1]]])
gcols <- setdiff(names(raw), c("#FID", "FID", "IID", "#IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G <- as.matrix(raw[, ..gcols])
colnames(G) <- sub("_[^_]+$", "", gcols)
G <- G[, ch$ID, drop = FALSE]
G <- apply(G, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); (x - mean(x)) / sd(x) })
g <- as.numeric(G %*% ch$beta)
g <- (g - mean(g)) / sd(g) * sqrt(VAR_G)
set.seed(100000L * (scen == "B") + seed)
sm <- fread(file.path(qc_dir, "samples.tsv")); sm <- sm[match(iid, IID)]
if (scen == "A") {
  s <- rep(0, length(g)); e <- rnorm(length(g), sd = sqrt(1 - VAR_G))
} else {
  s <- unname(super_pop_shift[sm$super_pop]); s <- (s - mean(s)) / sd(s) * sqrt(VAR_S)
  e <- rnorm(length(g), sd = sqrt(1 - VAR_G - VAR_S))
}
y <- g + s + e
fwrite(data.table(`#IID` = iid, y = y), file.path(out, "pheno.tsv"), sep = "\t")
fwrite(data.table(scenario = scen, seed = seed, n = length(y), var_g = var(g), var_s = var(s),
                  var_y = var(y), true_h2 = var(g) / var(y)), file.path(out, "truth.csv"))
