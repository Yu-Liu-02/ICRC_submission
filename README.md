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
│   ├── EM_proposed.cpp        # C++ E-step and M-step 
│   ├── compile_cpp_em.R       # Compiles EM_proposed.cpp 
│   ├── EM_proposed.R          # E-step, M-step, EM loop
│   └── variance_estimation.R  # Profile likelihood variance estimation
├── Ahn_func/
│   └── EM_Ahn.R               # Functions for Ahn method
├── Goggins_func/
│   ├── Gibbs.cpp              # C++ Gibbs sampler
│   ├── compile_cpp_mcem.R     # Compiles Gibbs.cpp
│   └── MCEM_Goggins.R         # Functions for Goggins method
└── simulation/
    ├── Table1/                # Indicator g, beta = 1
    │   ├── Proposed.R         # Proposed method
    │   ├── Ahn.R              # Ahn method
    │   ├── Goggins.R          # Goggins method
    │   └── Midpoint.R         # Midpoint imputation
    └── Table2/                # ReLU g, beta = 0.1
        ├── Proposed.R         # Proposed method
        └── Midpoint.R         # Midpoint imputation
```

---

## Proposed Method

### 1. Setup

Source all required files and load all required packages (done automatically in `Proposed.R`):

```r
library(Rcpp)
library(RcppArmadillo)
library(bayess)

source("data_generation/data_generate.R")
source("data_generation/g_function.R")
source("Proposed_func/compile_cpp_em.R")   
source("Proposed_func/EM_proposed.R")
source("Proposed_func/variance_estimation.R")
```

### 2. Choose the g function

The model supports two forms of the pre-clinical effect:

```r
set_g_type("indicator")   # g(s, t) = 1(s >= t)  [default]
set_g_type("relu")        # g(s, t) = max(s - t, 0)
```

This must be called before running the EM. Both the R code and the C++ kernels read `g_type_int()` internally, so a single call propagates everywhere.

### 3. Prepare data

| Variable | Description |
|----------|-------------|
| `Y` | Observed outcome event time: `min(Z, C)`, where `Z` is the true outcome event time and `C` is the censoring time |
| `delta` | Outcome event indicator: 1 if `Z <= C` (event observed), 0 if censored |
| `L` | Left endpoint of the interval containing the preclinical event `TT` (last examination before `TT`; 0 if none) |
| `R` | Right endpoint of the interval containing the preclinical event `TT` (first examination after `TT`; `Inf` if never observed) |
| `X` | n × p covariate matrix |

The expected data structure is a `data.frame` with these columns, and `X` attached as `dat$X`:

```r
dat        <- data.frame(Y = Y, delta = delta, L = LR$L, R = LR$R)
dat$X      <- X                          # n x p covariate matrix
k.start    <- findInterval(dat$L, t) + 1 # index of first t >= L
k.end      <- findInterval(dat$R, t)     # index of last  t <= R
```

where `t` is the pre-clinical event time grid (interior points of the unique examination times) and `s` is the vector of unique observed outcome event times among uncensored subjects.

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


## Simulation Data Generation

All DGP functions are defined in `data_generation/data_generate.R`. A full simulation replicate is generated as follows:

```r
X  <- draw_X(n, p)                                      # covariates
TT <- rT_PH_bump(n, X, gamma, c, a, omega, b, tau, sigma)  # Pre-clinical event time
U  <- draw_examination(n, count = 3, X, psi)            # examination times
LR <- get_LR(TT, U)                                     # interval (L, R) for pre-clinical event time
Z  <- draw_Z(n, X, TT, alpha, beta, g_type)             # outcome event time
C  <- rexp(n, 1/6)                                      # censoring time
Y  <- ifelse(Z <= C, Z, C)                              # observed outcome time
delta <- ifelse(Z <= C, 1, 0)                           # outcome event indicator
```

### Function details

**`draw_X(n, p)`**
Draws an n × p covariate matrix: `X1 ~ Bernoulli(0.6)`, `X2, X3 ~ TruncNormal(0, 1, -2, 2)`.

**`rT_PH_bump(n, X, gamma, c, a, omega, b, tau, sigma)`**
Draws pre-clinical event times `TT` from a proportional hazards model with a bump baseline hazard.

**`draw_examination(n, count, X, psi)`**
Draws `count` covariate-dependent examination times per subject. 

**`get_LR(TT, U)`**
Converts examination times `U` and true pre-clinical event time `TT` into interval-censored brackets `(L, R)`.

**`draw_Z(n, X, TT, alpha, beta, g_type)`**
Draws outcome event times from a proportional hazards model with preclinical effect `g(Z, TT)`. Supports `g_type = "indicator"` and `g_type = "relu"`.

---



## Simulation Scripts

All simulation scripts follow the same structure:

1. Generate data (`draw_X`, `rT_PH_bump`, `draw_Z`, `draw_examination`, `get_LR`)
2. Fit the method
3. Store point estimates in `pe` / `coef` and SEs in `se` / `htsis`
4. After the loop, print a summary table of **Bias, SE, SEE, CP** for α and β

### Table 1 — Indicator g, β = 1

| Script | Method | Key function | SE source |
|--------|--------|--------------|-----------|
| `Proposed.R` | Proposed EM | `EM_proposed()` + `variance_est()` | Profile likelihood |
| `Ahn.R` | Ahn et al. (2018) | Newton–Raphson in `EM_Ahn.R` | `covProbf()` (model-based) |
| `Goggins.R` | Goggins MCEM | `fit_mcem_interval()` in `MCEM_Goggins.R` | `out$se` from MCEM |
| `Midpoint.R` | Midpoint imputation | `coxph()` with `tt()` | `sqrt(diag(fit$var))` |

True parameter values: γ = (1, 1.5, −1.5), α = (0.45, 0.5, −0.25), **β = 1**, `g_type = "indicator"`.

**Ahn:** estimates α and β only (Parametric Weibull model for preclinical event fit separately via `survreg`). Output is `p+1` columns.

**Goggins:** MCEM with a Gibbs sampler (C++ backend in `Gibbs.cpp`). Starting values from a Cox model with midpoint imputation. Output is `p+1` columns (α, β).

**Midpoint:** Cox model with time-varying covariate `g(t, midpoint)`. Uses the global `g` function set by `set_g_type()`. Output is `p+1` columns (α, β).

### Table 2 — ReLU g, β = 0.1

| Script | Method | Key function | SE source |
|--------|--------|--------------|-----------|
| `Proposed.R` | Proposed EM | `EM_proposed()` + `variance_est()` | Profile likelihood |
| `Midpoint.R` | Midpoint imputation | `coxph()` with `tt()` | `sqrt(diag(fit$var))` |

True parameter values: γ = (1, 1.5, −1.5), α = (0.45, 0.5, −0.25), **β = 0.1**, `g_type = "relu"`.

Ahn and Goggins are not included in Table 2 as their method do not support the ReLU g function.

### Switching g type in simulations

`Proposed.R` and `Midpoint.R` call `set_g_type()` near the top of the script. Change it to switch scenarios:

```r
set_g_type("indicator")   # Table 1 setting
set_g_type("relu")        # Table 2 setting
```

---

## Dependencies

- R packages: `Rcpp`, `RcppArmadillo`, `bayess`, `survival`
- A C++ compiler accessible to R (e.g. via Xcode CLT on macOS, or `gcc`/`gfortran` on Linux)
- On macOS with Homebrew GCC, ensure `~/.R/Makevars` points to the correct gfortran library path
