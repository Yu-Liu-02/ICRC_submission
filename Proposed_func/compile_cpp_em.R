library(Rcpp)
library(RcppArmadillo)

.find_em_cpp <- function() {
  candidates <- c(
    "EM_proposed.cpp",
    "Proposed_func/EM_proposed.cpp",
    "../Proposed_func/EM_proposed.cpp",
    "../../Proposed_func/EM_proposed.cpp",
    "../../../Proposed_func/EM_proposed.cpp"
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  NULL
}

cpp_file <- .find_em_cpp()
if (is.null(cpp_file)) stop("EM_proposed.cpp not found.")

# Use a per-process cache directory to avoid race conditions when multiple
# SLURM array jobs compile simultaneously and corrupt the shared cache index.
cache_dir <- file.path(tempdir(), paste0("rcpp_cache_", Sys.getpid()))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("Compiling EM Rcpp code...\n")
sourceCpp(cpp_file, cacheDir = cache_dir)
cat("Compilation successful and cached!\n")
