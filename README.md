# Code Overview

Joint model for interval-censored change point (IC) and recurrent event (RC) outcomes, estimated via EM algorithm with profile likelihood variance estimation.

---

## Repository Structure

```
ICRC_submission/
├── data_generation/
│   ├── data_generate.R        # Data-generating functions
│   └── g_function.R           # g function (indicator / ReLU)
├── Proposed_func/
│   ├── EM_proposed.cpp        # C++ E-step and M-step kernels
│   ├── compile_cpp_em.R       # Compiles EM_proposed.cpp via Rcpp
│   ├── EM_proposed.R          # E-step, M-step, EM loop
│   └── variance_estimation.R  # Profile likelihood variance estimation
├── Ahn_func/
│   └── EM_Ahn.R               # Ahn competitor method
├── Goggins_func/
│   ├── Gibbs.cpp              # C++ Gibbs sampler
│   ├── compile_cpp_mcem.R     # Compiles Gibbs.cpp
│   └── MCEM_Goggins.R         # Goggins competitor method
└── simulation/
    └── Table1/
        ├── Proposed.R         # Simulation: proposed method
        ├── Ahn.R              # Simulation: Ahn method
        ├── Goggins.R          # Simulation: Goggins method
        └── Midpoint.R         # Simulation: midpoint imputation
```

---

## Proposed Method

### 1. Setup

Source all required files (done automatically in `Proposed.R`):

```r
library(Rcpp)
library(RcppArmadillo)
library(bayess)

source("data_generation/data_generate.R")
source("data_generation/g_function.R")
source("Proposed_func/compile_cpp_em.R")   # compiles C++ on first run
source("Proposed_func/EM_proposed.R")
source("Proposed_func/variance_estimation.R")
```

### 2. Choose the g function

The model supports two forms of the threshold effect:

```r
set_g_type("indicator")   # g(s, t) = 1(s >= t)  [default]
set_g_type("relu")        # g(s, t) = max(s - t, 0)
```

This must be called before running the EM. Both the R code and the C++ kernels read `g_type_int()` internally, so a single call propagates everywhere.

### 3. Prepare data

The expected data structure is a `data.frame` with columns `Y`, `delta`, `L`, `R`, and a matrix `X` attached as `dat$X`:

```r
dat        <- data.frame(Y = Y, delta = delta, L = LR$L, R = LR$R)
dat$X      <- X                          # n x p covariate matrix
k.start    <- findInterval(dat$L, t) + 1 # index of first t >= L
k.end      <- findInterval(dat$R, t)     # index of last  t <= R
```

where `t` is the IC event time grid (interior points of the unique examination times) and `s` is the vector of unique observed RC event times among uncensored subjects.

### 4. Initialize parameters

```r
m1 <- length(t); m2 <- length(s)
params <- list(
  IC = list(gamma     = rep(0, p),
            lambda    = rep(1/m1, m1),
            cumlambda = c(0, cumsum(rep(1/m1, m1)))),
  RC = list(alpha = rep(0, p),
            beta  = 0,
            h     = rep(1/m2, m2))
)
```

### 5. Run the EM

```r
em_fit <- EM_proposed(dat, t, s, params, k.start, k.end,
                      M = 8000, epsilon = 5e-4)
```

Returns a list:
- `em_fit$params` — converged parameter estimates (`$IC$gamma`, `$IC$lambda`, `$RC$alpha`, `$RC$beta`, `$RC$h`)
- `em_fit$E` — final E-step quantities (`$p`, `$q`, `$w`)

### 6. Variance estimation

Profile likelihood SEs for all regression parameters (γ, α, β):

```r
var_fit <- variance_est(dat, t, s, em_fit, k.start, k.end, n,
                        M = 8000, epsilon = 5e-4)
```

Returns:
- `var_fit$se` — SEs in order: `gamma[1:p]`, `alpha[1:p]`, `beta`

---

## Simulation Scripts

All four simulation scripts in `simulation/Table1/` follow the same structure:

1. Generate data (`draw_X`, `rT_PH_bump`, `draw_Z`, `draw_examination`, `get_LR`)
2. Fit the method
3. Store point estimates in `pe` / `coef` and SEs in `se` / `htsis`
4. After the loop, print a summary table of **Bias, SE, SEE, CP** for α and β

### True parameter values (all scripts)

| Parameter | Value |
|-----------|-------|
| γ | (1, 1.5, −1.5) |
| α | (0.45, 0.5, −0.25) |
| β | 1 |

### Alternative methods

| Script | Method | Key function | SE source |
|--------|--------|--------------|-----------|
| `Proposed.R` | Proposed EM | `EM_proposed()` + `variance_est()` | Profile likelihood |
| `Ahn.R` | Ahn et al. (2018) | Newton–Raphson loop in `EM_Ahn.R` | `covProbf()` (model-based) |
| `Goggins.R` | Goggins MCEM | `fit_mcem_interval()` in `MCEM_Goggins.R` | `out$se` from MCEM |
| `Midpoint.R` | Midpoint imputation | `coxph()` with `tt()` | `sqrt(diag(fit$var))` |

**Ahn:** estimates α and β jointly (no γ; Weibull IC model fit separately via `survreg`). Output is `p+1` columns.

**Goggins:** MCEM with a Gibbs sampler (C++ backend in `Gibbs.cpp`). Starting values are obtained from a Cox model with midpoint imputation. Output is `p+1` columns (α, β).

**Midpoint:** Cox model with a time-varying covariate `g(t, midpoint)`. Uses the global `g` function, so `set_g_type()` controls whether indicator or ReLU is applied. Output is `p+1` columns (α, β).

### Switching g type in simulations

Only `Proposed.R` and `Midpoint.R` use `set_g_type()`. Change the call near the top of each script:

```r
set_g_type("relu")       # ReLU threshold effect
set_g_type("indicator")  # Hard threshold effect (default)
```

`Ahn.R` and `Goggins.R` do not model the threshold effect via g; they assume a simpler parametric structure for the IC event.

---

## Dependencies

- R packages: `Rcpp`, `RcppArmadillo`, `bayess`, `survival`
- A C++ compiler accessible to R (e.g. via Xcode CLT on macOS, or `gcc`/`gfortran` on Linux)
- On macOS with Homebrew GCC, ensure `~/.R/Makevars` points to the correct gfortran library path
