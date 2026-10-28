source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "R", "common.R"))
A <- as.matrix(fread(file.path(qc_dir, "grm.rel"), header = FALSE))
ids <- fread(file.path(qc_dir, "grm.rel.id"), header = FALSE)
iid <- as.character(ids[[ncol(ids)]]); if (iid[1] == "IID" || startsWith(iid[1], "#")) iid <- iid[-1]
stopifnot(nrow(A) == ncol(A), nrow(A) == length(iid))
e <- eigen(A, symmetric = TRUE)
saveRDS(list(U = e$vectors, lam = e$values, iid = iid), file.path(qc_dir, "grm_eigen.rds"))
cat(sprintf("GRM %d x %d; mean diagonal %.3f\n", nrow(A), ncol(A), mean(diag(A))))
