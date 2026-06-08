library(Rcpp)

# Locate Gibbs.cpp by trying candidate paths relative to working directory.
# This handles being called from any simulation subfolder in the submission.
.find_gibbs_cpp <- function() {
  candidates <- c(
    "Gibbs.cpp",
    "Goggins_func/Gibbs.cpp",
    "../Goggins_func/Gibbs.cpp",
    "../../Goggins_func/Gibbs.cpp",
    "../../../Goggins_func/Gibbs.cpp"
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  NULL
}

cpp_file <- .find_gibbs_cpp()
if (is.null(cpp_file)) stop("Gibbs.cpp not found. Please run from within the submission folder.")

# Use a per-process cache directory to avoid race conditions when multiple
# SLURM array jobs compile simultaneously and corrupt the shared cache index.
cache_dir <- file.path(tempdir(), paste0("rcpp_cache_", Sys.getpid()))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("Compiling Gibbs sampler Rcpp code...\n")
sourceCpp(cpp_file, cacheDir = cache_dir)
cat("Compilation successful!\n")
