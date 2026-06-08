suppressPackageStartupMessages({
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
})

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

.is_null <- function(x) is.null(x) || length(x) == 0
`%||%` <- function(a, b) if (.is_null(a)) b else a


.has_cpp_full_gibbs <- function() exists("gibbs_sample_k_with_w_cpp", mode = "function")

make_T_endpoints <- function(L, R, time, status, include_zero = FALSE) {
  stopifnot(length(time) == length(status), length(L) == length(R), length(L) == length(time))
  Tcore <- c(L, R[is.finite(R)], time[status == 1L])
  if (!include_zero) Tcore <- Tcore[Tcore > 0]
  Tcore <- sort(unique(Tcore))
  if (length(Tcore) == 0L) stop("Empty time grid: check L/R/time inputs.")
  c(Tcore, Inf) 
}


feasible_k_T <- function(Li, Ri, T) {
  J <- length(T)
  if (is.infinite(Ri)) {
    which(Li < T)  
  } else {
    which(Li < T & T <= Ri)
  }
}


Z_at_time_T <- function(t, k_i, T) {
  if (k_i >= length(T)) return(0L)      
  as.integer(t >= T[k_i])
}

.make_event_info <- function(time, status) {
  evt_idx <- which(status == 1L)
  if (length(evt_idx) == 0L) stop("No events (status==1) present.")
  o <- order(time[evt_idx])
  evt_idx <- evt_idx[o]
  list(evt_idx = evt_idx, evt_time = time[evt_idx])
}
.riskset_indices <- function(time, t) which(time >= t)

# ---------------------------------------
# Cox partial loglik given k (R version) 
# ---------------------------------------
ploglik_ab <- function(alpha, beta, time, status, X, k, T) {
  X <- as.matrix(X)
  n <- length(time)
  stopifnot(length(status) == n, nrow(X) == n, length(k) == n)
  ev <- .make_event_info(time, status)
  ll <- 0
  for (e in seq_along(ev$evt_idx)) {
    i <- ev$evt_idx[e]
    t <- ev$evt_time[e]
    Rset <- .riskset_indices(time, t)
    
    Zi <- Z_at_time_T(t, k[i], T)
    eta_i <- sum(X[i, ] * alpha) + beta * Zi  
    
    Z_R <- vapply(Rset, function(j) Z_at_time_T(t, k[j], T), integer(1))
    eta_R <- as.vector(X[Rset, , drop = FALSE] %*% alpha) + beta * Z_R 
    
    ll <- ll + eta_i - log(sum(exp(eta_R)))
  }
  ll
}

# ---------------------------------------
# Multinomial MLE for w from MC draws
# ---------------------------------------
estimate_w_mle <- function(draws, J) {
  B <- nrow(draws); n <- ncol(draws)
  w_acc <- rep(0, J)
  for (b in seq_len(B)) {
    tab <- tabulate(draws[b, ], nbins = J)
    w_acc <- w_acc + tab / n
  }
  w_hat <- w_acc / B
  w_hat <- pmax(w_hat, 1e-12)
  w_hat / sum(w_hat)
}

# ---------------------------------------
# Score/Hessian
# ---------------------------------------
score_hess_ab <- function(alpha, beta, time, status, X, k, T) {
  X <- as.matrix(X)
  p <- length(alpha)
  ev <- .make_event_info(time, status)
  
  U <- numeric(p + 1)
  H <- matrix(0, p + 1, p + 1)
  
  for (e in seq_along(ev$evt_idx)) {
    i <- ev$evt_idx[e]
    t <- ev$evt_time[e]
    Rset <- .riskset_indices(time, t)
    
    Z_R <- vapply(Rset, function(j) Z_at_time_T(t, k[j], T), integer(1))
    
    C_R <- cbind(X[Rset, , drop = FALSE], Z_R)
    
    eta_R <- as.vector(X[Rset, , drop = FALSE] %*% alpha) + beta * Z_R
    wgt <- exp(eta_R)
    S0 <- sum(wgt)
    
    S1 <- colSums(wgt * C_R)
    S2 <- crossprod(C_R, wgt * C_R)
    
    Zi <- Z_at_time_T(t, k[i], T)
    C_i <- c(X[i, ], Zi)
    
    E_C <- S1 / S0
    U <- U + (C_i - E_C)
    
    V <- (S2 / S0) - tcrossprod(E_C)
    H <- H - V
  }
  
  list(score = U, hess = H)
}

variance_louis_ab_only <- function(alpha, beta, time, status, X, draws, T) {
  X <- as.matrix(X)
  p <- length(alpha)
  B <- nrow(draws)
  
  S <- matrix(0, nrow = B, ncol = p + 1)
  H <- array(0, dim = c(p + 1, p + 1, B))
  
  for (b in seq_len(B)) {
    k_b <- draws[b, ]
    sh <- score_hess_ab(alpha, beta, time, status, X, k_b, T)
    S[b, ] <- sh$score
    H[,,b] <- sh$hess
  }
  
  E_H <- apply(H, c(1,2), mean)
  S_center <- sweep(S, 2, colMeans(S), "-")
  V_S <- crossprod(S_center) / (B - 1)
  
  Iobs <- -E_H - V_S
  V <- solve(Iobs)
  
  param_names <- c(paste0("alpha", 1:p), "beta")
  colnames(V) <- rownames(V) <- param_names
  U_mean <- colMeans(S)
  names(U_mean) <- param_names

  
  list(vcov = V, Iobs = Iobs, score_mc_mean = U_mean)
}

# ---------------------------------------
# Rcpp loader
# ---------------------------------------
.load_rcpp_backend <- function() {
  if (.has_cpp_full_gibbs()) return(TRUE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
  cpp <- .find_gibbs_cpp()
  if (is.null(cpp)) return(FALSE)
  Rcpp::sourceCpp(cpp, verbose = FALSE)
  TRUE
}

# ---------------------------------------
# E-step Gibbs sampler
# ---------------------------------------
mcem_gibbs_sample_k_with_w <- function(alpha, beta, w, time, status, X, k_init, T, L, R,
                                       burn = 250, thin = 1, keep = 250, seed = NULL) {
  X <- as.matrix(X)
  n <- length(time)
  J <- length(T)
  stopifnot(length(w) == J)
  
  feasible_list <- vector("list", n)
  for (i in seq_len(n)) {
    ks <- feasible_k_T(L[i], R[i], T)
    if (length(ks) == 0L) stop("Empty feasible set for i=", i)
    feasible_list[[i]] <- ks
  }
  .load_rcpp_backend()
  if (.has_cpp_full_gibbs()) {
    out <- gibbs_sample_k_with_w_cpp(alpha, beta, w,
                                     time, as.integer(status), X,
                                     as.integer(k_init),
                                     T, feasible_list,
                                     burn = as.integer(burn),
                                     thin = as.integer(thin),
                                     keep = as.integer(keep),
                                     seed = as.integer(if (is.null(seed)) 1L else seed))
    return(list(draws = out$draws, k_last = out$k_last, feasible = feasible_list))
  }
  
  # R Fallback
  k <- k_init
  draws <- matrix(NA_integer_, nrow = keep, ncol = n)
  saved <- 0L
  total_iters <- burn + keep * thin
  
  for (iter in seq_len(total_iters)) {
    for (i in seq_len(n)) {
      cand <- feasible_list[[i]]
      k_old <- k[i]
      lp <- numeric(length(cand))
      for (cc in seq_along(cand)) {
        k[i] <- cand[cc]
        lp[cc] <- ploglik_ab(alpha, beta, time, status, X, k, T) + log(w[k[i]])
      }
      k[i] <- k_old
      lp <- lp - max(lp)
      p <- exp(lp); p <- p / sum(p)
      k[i] <- cand[sample.int(length(cand), 1L, prob = p)]
    }
    if (iter > burn && ((iter - burn) %% thin == 0L)) {
      saved <- saved + 1L
      draws[saved, ] <- k
      if (saved >= keep) break
    }
  }
  list(draws = draws, k_last = k, feasible = feasible_list)
}

# ---------------------------------------
# M-step for (alpha,beta)
# ---------------------------------------
mstep_optimize <- function(time, status, X, draws, TT, start) {
  p <- ncol(X)
  f <- function(par) {
    a_vec <- par[1:p]
    b_val <- par[p + 1]
    mean_ll <- mean_ploglik_ab_cpp(a_vec, b_val, time, as.integer(status), X, draws, TT)
    -mean_ll 
  }
  opt <- optim(par = start, fn = f, method = "BFGS")
  list(alpha = opt$par[1:p], beta = opt$par[p + 1], value = opt$value, converged = (opt$convergence == 0L))
}

# ---------------------------------------
# Main fit
# ---------------------------------------
fit_mcem_interval <- function(time, status, X, L, R,
                                     em_max = 150, tol = 1e-3,
                                     burn = 250, thin = 1, keep = 250,
                                     seed = 1,
                                     start = NULL,
                                     verbose = TRUE) {
  
  X <- as.matrix(X)
  p <- ncol(X)
  n <- length(time)
  stopifnot(length(status) == n, nrow(X) == n, length(L) == n, length(R) == n)
  
  TT <- make_T_endpoints(L, R, time, status)
  J <- length(TT)
  
  k_init <- integer(n)
  for (i in seq_len(n)) {
    ks <- feasible_k_T(L[i], R[i], TT)
    mid <- if (is.infinite(R[i])) L[i] else 0.5 * (L[i] + R[i])
    d <- rep(Inf, length(ks))
    for (jj in seq_along(ks)) {
      if (is.finite(TT[ks[jj]])) d[jj] <- abs(TT[ks[jj]] - mid)
      else if (is.infinite(R[i])) d[jj] <- 0 
    }
    k_init[i] <- ks[which.min(d)]
  }
  
  if (is.null(start)) start <- rep(0, p + 1)
  alpha <- start[1:p]
  beta <- start[p + 1]
  k_last <- k_init
  w <- rep(1 / J, J)
  
  trace <- data.frame()
  
  for (m in seq_len(em_max)) {
    if (verbose) {
      message(sprintf("EM %d | beta=%.4g | w[1:3]=%s",
                      m, beta, paste(signif(w[1:min(3,J)], 3), collapse = ",")))
    }
    
    gib <- mcem_gibbs_sample_k_with_w(alpha, beta, w, time, status, X, k_last, TT, L, R,
                                      burn = burn, thin = thin, keep = keep, seed = seed + m)
    draws <- gib$draws
    k_last <- gib$k_last
    
    w <- estimate_w_mle(draws, J)
    
    opt <- mstep_optimize(time, status, X, draws, TT, c(alpha, beta))
    alpha_new <- opt$alpha; beta_new <- opt$beta
    
    # Save trace
    trace_row <- data.frame(iter = m)
    for(d in 1:p) trace_row[[paste0("alpha", d)]] <- alpha_new[d]
    trace_row$beta <- beta_new
    trace <- rbind(trace, trace_row)
    
    if (max(abs(c(alpha_new - alpha, beta_new - beta))) < tol) {
      alpha <- alpha_new; beta <- beta_new
      if (verbose) message("Converged.")
      break
    }
    alpha <- alpha_new; beta <- beta_new
  }
  print(c(alpha, beta))
  
  gib_fin <- mcem_gibbs_sample_k_with_w(alpha, beta, w, time, status, X, k_last, TT, L, R,
                                        burn = burn, thin = thin, keep = keep, seed = seed + 999)
  draws_fin <- gib_fin$draws
  w <- estimate_w_mle(draws_fin, J)
  
  louis_ab <- variance_louis_ab_only(alpha, beta, time, status, X, draws_fin, TT)
  Vab <- louis_ab$vcov
  se_ab <- sqrt(diag(Vab))
  
  param_vals <- c(alpha, beta)
  names(param_vals) <- c(paste0("alpha", 1:p), "beta")
  
  z_ab <- param_vals / se_ab
  p_ab <- 2 * stats::pnorm(-abs(z_ab))
  
  list(
    coef = param_vals,
    vcov_ab = Vab,
    se = se_ab,
    z = z_ab,
    p = p_ab,
    TT = TT,          
    w = w,          
    trace = trace,
    draws_final = draws_fin,
    louis_ab = louis_ab
  )
}