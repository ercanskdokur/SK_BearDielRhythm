# Parse-check every analysis script (no execution) — quick pre-submit sanity check.
fs <- list.files(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"),
                           "BrownBearDielAct", "scripts"),
                 pattern = "[.]R$", full.names = TRUE)
fs <- fs[!grepl("_syntax_check", fs)]
bad <- 0
for (f in fs) {
  r <- tryCatch({ parse(f); "OK" }, error = function(e) conditionMessage(e))
  if (identical(r, "OK")) cat(sprintf("  OK     %s\n", basename(f)))
  else { bad <- bad + 1; cat(sprintf("  ERROR  %s\n        -> %s\n", basename(f), r)) }
}
cat(sprintf("\n%d files, %d parse error(s)\n", length(fs), bad))
