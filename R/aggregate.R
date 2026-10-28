source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "R", "common.R"))
f <- list.files(file.path(root, "results", "tasks"), pattern = "^summary.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(f) > 0)
all <- rbindlist(lapply(f, fread))[order(scenario, seed)]
fwrite(all, file.path(root, "results", "tables", "summary.csv"))
num <- setdiff(names(all)[sapply(all, is.numeric)], c("seed", "n_causal"))
agg <- all[, c(list(reps = .N), lapply(.SD, mean)), by = scenario, .SDcols = num]
fwrite(agg, file.path(root, "results", "tables", "summary_by_scenario.csv"))
print(agg)
