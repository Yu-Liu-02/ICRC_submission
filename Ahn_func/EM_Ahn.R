#######################################################
#  Faster sub-functions for Ahn method — multivariate z
#######################################################

.fast_subject_integrals <- function(lower, upper, zind, beta1, beta2, Lambda0) {
  # Returns cumulative integrals F(j) = int_G(lower, tt[j]) for all tt[j] in
  # [lower, upper], together with total int_G(lower, upper).
  # This replaces many repeated intG.EM() calls for the same subject.

  uconst  <- exp(sum(zind * beta1)) * (exp(beta2) - 1)
  lambind <- exp(h0) * exp(sum(zind * gamma))

  sind <- sum(tt < lower) + 1L
  eind <- sum(tt < upper)

  # If there is no event time between lower and upper, closed form is trivial.
  if (eind < sind) {
    total <- exp((if (sind == 1L) 0 else Lambda0[sind - 1L]) * uconst) *
      (exp(-lambind * lower^alpha) - exp(-lambind * upper^alpha))
    return(list(sind = sind, eind = eind, F = numeric(0), total = total))
  }

  idx <- sind:eind
  m   <- length(idx)

  starts <- c(lower, tt[idx[-m]])
  ends   <- tt[idx]
  aa_vec <- c(if (sind == 1L) 0 else Lambda0[sind - 1L], Lambda0[idx[-m]])

  contrib <- exp(aa_vec * uconst) *
    (exp(-lambind * starts^alpha) - exp(-lambind * ends^alpha))
  F <- cumsum(contrib)

  aa_last <- Lambda0[eind]
  tail_contrib <- exp(aa_last * uconst) *
    (exp(-lambind * tt[eind]^alpha) - exp(-lambind * upper^alpha))

  list(sind = sind, eind = eind, F = F, total = F[m] + tail_contrib)
}


intG.EM <- function(beta, lower, upper, zind, Lambda0) {
  # Compatibility wrapper.
  p     <- ncol(z)
  beta1 <- beta[1:p]
  beta2 <- beta[p + 1L]
  .fast_subject_integrals(lower, upper, zind, beta1, beta2, Lambda0)$total
}


lambda0.EM <- function(beta, Dout) {
  p     <- ncol(z)
  beta1 <- beta[1:p]
  beta2 <- beta[p + 1L]

  zbeta <- as.vector(z %*% beta1)
  ez    <- exp(zbeta)

  eD <- Dout[1:N, , drop = FALSE]
  D0 <- Dout[(N + 1L):(2L * N), , drop = FALSE]

  # Vectorized replacement of the inner j-loop.
  obs_term <- (Y * R) * exp(zbeta + D * beta2)
  unk_term <- (Y * (1 - R)) * (ez * D0)

  denom <- colSums(obs_term + unk_term)
  lambda0 <- delta / denom
  cumsum(lambda0)
}


expD.EM <- function(beta, cumlambda0) {
  p     <- ncol(z)
  beta1 <- beta[1:p]
  beta2 <- beta[p + 1L]
  eb2   <- exp(beta2)

  Lambda0 <- cumlambda0
  eD      <- matrix(0, N, N)

  for (i in seq_len(N)) {
    zind    <- z[i, ]
    zbeta_i <- sum(zind * beta1)
    lambind <- exp(h0) * exp(sum(zind * gamma))

    ints <- .fast_subject_integrals(x.left[i], x.right[i], zind, beta1, beta2, Lambda0)
    sind <- ints$sind
    eind <- ints$eind
    F    <- ints$F
    Gtot <- ints$total

    # Regions outside [L_i, R_i]
    if (sind > 1L) eD[i, seq_len(sind - 1L)] <- 0
    if (eind < N)  eD[i, (eind + 1L):N] <- 1

    # Nothing to do inside interval if there are no internal tt points.
    if (length(F) == 0L) next

    idx <- sind:eind

    if (tt[i] > x.right[i]) {
      # Figure 3 in the original code.
      eD[i, idx] <- F / Gtot
    } else if (tt[i] < x.right[i]) {
      # Figure 4 in the original code.
      num_base <- exp(delta[i] * beta2) * exp(-Lambda0[i] * exp(zbeta_i + beta2))
      denom <- num_base * intG.EM(beta, x.left[i], tt[i], zind, cumlambda0) +
        exp(-Lambda0[i] * exp(zbeta_i)) *
        (exp(-lambind * tt[i]^alpha) - exp(-lambind * x.right[i]^alpha))

      # Only columns j with tt[j] <= tt[i] can be nonzero.
      upto <- min(eind, i)
      if (upto >= sind) {
        loc <- seq_len(upto - sind + 1L)
        eD[i, sind:upto] <- num_base * F[loc] / denom
      }
    }
  }

  D0 <- eD * eb2 + 1 - eD
  D1 <- eD * eb2

  Dout <- matrix(0, 3L * N, N)
  Dout[1:N, ]                <- eD
  Dout[(N + 1L):(2L * N), ]  <- D0
  Dout[(2L * N + 1L):(3L * N), ] <- D1
  Dout
}


Sobs.EM <- function(beta, xt.data, Dout) {
  p     <- ncol(z)
  beta1 <- beta[1:p]
  beta2 <- beta[p + 1L]

  zbeta <- as.vector(z %*% beta1)
  ez    <- exp(zbeta)

  eD <- Dout[1:N, , drop = FALSE]
  D0 <- Dout[(N + 1L):(2L * N), , drop = FALSE]
  D1 <- Dout[(2L * N + 1L):(3L * N), , drop = FALSE]

  eeta_obs  <- exp(zbeta + D * beta2)
  eeta_base <- ez

  w_obs <- Y * R
  w_unk <- Y * (1 - R)

  Aobs <- w_obs * eeta_obs
  Aunk0 <- w_unk * (eeta_base * D0)
  Aunk1 <- w_unk * (eeta_base * D1)

  S0 <- (colSums(Aobs) + colSums(Aunk0)) / N

  S1 <- matrix(0, N, p + 1L)
  Zweighted_obs_base <- Aobs
  for (k in seq_len(p)) {
    zk <- z[, k]
    S1[, k] <- (colSums(Zweighted_obs_base * zk) + colSums(Aunk0 * zk)) / N
  }
  S1[, p + 1L] <- (colSums(Aobs * D) + colSums(Aunk1)) / N

  S2 <- array(0, c(N, p + 1L, p + 1L))
  for (k in seq_len(p)) {
    zk <- z[, k]
    for (l in k:p) {
      zl <- z[, l]
      v <- (colSums(Aobs * (zk * zl)) + colSums(Aunk0 * (zk * zl))) / N
      S2[, k, l] <- v
      S2[, l, k] <- v
    }
    v <- (colSums(Aobs * (zk * D)) + colSums(Aunk1 * zk)) / N
    S2[, k, p + 1L] <- v
    S2[, p + 1L, k] <- v
  }
  S2[, p + 1L, p + 1L] <- S1[, p + 1L]

  list(S0 = S0, S1 = S1, S2 = S2)
}


score.EM <- function(beta, Sout, Dout, cumlambda0) {
  p  <- ncol(z)
  S0 <- Sout$S0
  S1 <- Sout$S1

  eD <- Dout[1:N, , drop = FALSE]
  diag_eD <- diag(eD)
  dR <- diag(R)
  dD <- diag(D)

  A <- numeric(p + 1L)
  A[1:p] <- colSums(delta * (z - S1[, 1:p, drop = FALSE] / S0)) / N
  A[p + 1L] <- sum(delta * (dR * dD + (1 - dR) * diag_eD - S1[, p + 1L] / S0)) / N

  matrix(A, p + 1L, 1L)
}


hessian.EM <- function(beta, Sout, Dout, cumlambda0) {
  p  <- ncol(z)
  S0 <- Sout$S0
  S1 <- Sout$S1
  S2 <- Sout$S2

  H <- matrix(0, p + 1L, p + 1L)
  invS0sq <- 1 / (S0^2)
  for (k in seq_len(p + 1L)) {
    for (l in k:(p + 1L)) {
      val <- sum(delta * (S2[, k, l] * S0 - S1[, k] * S1[, l]) * invS0sq) / N
      H[k, l] <- val
      H[l, k] <- val
    }
  }
  H
}


covProbf <- function() {
  p            <- ncol(z)
  beta_em      <- as.vector(beta.EM)
  beta_new_vec <- as.vector(beta.new)
  beta2        <- beta_em[p + 1L]

  d.cumlambda0.EM <- c(cumlambda0.hat.EM[1L], diff(cumlambda0.hat.EM))
  z_eta <- as.vector(z %*% beta_em[1:p])

  k_lo   <- integer(N)
  is_obs <- logical(N)
  prLi   <- vector("list", N)

  for (i in seq_len(N)) {
    lambind <- exp(h0) * exp(sum(z[i, ] * gamma))

    if (sum(R[i, ] == 0L) == 0L) {
      ones      <- which(D[i, ] == 1L)
      k_lo[i]   <- if (length(ones) > 0L) min(ones) else N + 1L
      is_obs[i] <- TRUE
      prLi[[i]] <- 1
    } else {
      uncertain <- which(R[i, ] == 0L)
      k_hi_i    <- max(uncertain)
      k_lo[i]   <- min(uncertain)
      is_obs[i] <- FALSE
      n.paths   <- k_hi_i - k_lo[i] + 2L

      CL_i     <- cumlambda0.hat.EM[i]
      e_zeta_i <- exp(z_eta[i])
      e_beta2  <- exp(beta2)

      prWi_vec <- numeric(n.paths)
      for (row.k in seq_len(n.paths)) {
        k <- k_lo[i] + row.k - 1L
        if (k > i) {
          cum_haz <- e_zeta_i * CL_i
        } else {
          CL_km1  <- if (k == 1L) 0 else cumlambda0.hat.EM[k - 1L]
          cum_haz <- e_zeta_i * (CL_km1 + e_beta2 * (CL_i - CL_km1))
        }

        eta_ii <- z_eta[i] + beta2 * (i >= k)
        prWi   <- (d.cumlambda0.EM[i] * exp(eta_ii))^delta[i] * exp(-cum_haz)

        if (k == 1L) {
          p_T <- 1 - exp(-lambind * tt[1L]^alpha)
        } else if (k == k_hi_i + 1L) {
          p_T <- exp(-lambind * tt[k - 1L]^alpha)
        } else {
          p_T <- exp(-lambind * tt[k - 1L]^alpha) - exp(-lambind * tt[k]^alpha)
        }
        prWi_vec[row.k] <- prWi * p_T
      }

      s         <- sum(prWi_vec)
      prLi[[i]] <- if (s > 0) prWi_vec / s else rep(1 / n.paths, n.paths)
    }
  }

  idx_nonzero  <- which(d.cumlambda0.EM != 0)
  q            <- length(idx_nonzero)
  p_weib       <- dim(fit$var)[1]
  dim.EM       <- (p + 1L) + q + p_weib
  base_idx     <- (p + 1L) + seq_len(q)

  neg_solve_fitvar <- if (p_weib > 0L) -solve(fit$var) else NULL

  part  <- mat.EM3 <- matrix(0, dim.EM, 1L)
  part2 <- mat.EM1 <- mat.EM2 <- matrix(0, dim.EM, dim.EM)

  MC.B <- 500L
  for (MC.ite in seq_len(MC.B)) {
    k_vec <- k_lo
    for (i in seq_len(N)) {
      if (!is_obs[i]) {
        idx <- sample.int(length(prLi[[i]]), size = 1L, prob = prLi[[i]])
        k_vec[i] <- k_lo[i] + idx - 1L
      }
    }

    part1.sub <- numeric(p)
    part2.sub <- 0
    part3.sub <- numeric(N)
    part4.sub <- matrix(0, p, p)
    part5.sub <- 0
    part6.sub <- numeric(p)
    part7.sub <- matrix(0, N, p)
    part8.sub <- numeric(N)

    for (l in seq_len(N)) {
      d_l   <- as.integer(k_vec <= l)
      eta_l <- z_eta + d_l * beta2
      e_eta <- exp(eta_l)
      w_l   <- Y[, l] * d.cumlambda0.EM[l] * e_eta

      part1.sub <- part1.sub + as.vector(crossprod(z, w_l))
      part2.sub <- part2.sub + sum(w_l * d_l)
      part4.sub <- part4.sub + crossprod(z, w_l * z)
      part5.sub <- part5.sub + sum(w_l * d_l)
      part6.sub <- part6.sub + as.vector(crossprod(z, w_l * d_l))
      part3.sub[l]   <- sum(Y[, l] * e_eta)
      part7.sub[l, ] <- as.vector(crossprod(z, Y[, l] * e_eta))
      part8.sub[l]   <- sum(Y[, l] * d_l * e_eta)
    }

    d_diag <- as.integer(k_vec <= seq_len(N))
    part[1:p, 1L] <- colSums(delta * z) - part1.sub
    part[p + 1L, 1L] <- sum(delta * d_diag) - part2.sub
    if (q > 0L) {
      part[base_idx, 1L] <- (delta[idx_nonzero] / d.cumlambda0.EM[idx_nonzero]) -
        part3.sub[idx_nonzero]
    }

    part2[, ] <- 0
    part2[1:p, 1:p] <- -part4.sub
    part2[p + 1L, p + 1L] <- -part5.sub
    part2[1:p, p + 1L] <- -part6.sub
    part2[p + 1L, 1:p] <- -part6.sub
    if (q > 0L) {
      part2[1:p, base_idx] <- -t(part7.sub[idx_nonzero, , drop = FALSE])
      part2[p + 1L, base_idx] <- -part8.sub[idx_nonzero]
      part2[base_idx, 1:p] <- t(part2[1:p, base_idx, drop = FALSE])
      part2[base_idx, p + 1L] <- part2[p + 1L, base_idx]
      diag(part2)[base_idx] <- -(delta[idx_nonzero] / d.cumlambda0.EM[idx_nonzero])^2
    }
    if (p_weib > 0L) {
      w.idx <- (dim.EM - p_weib + 1L):dim.EM
      part2[w.idx, w.idx] <- neg_solve_fitvar
    }

    mat.EM1 <- mat.EM1 + part2
    mat.EM2 <- mat.EM2 + part %*% t(part)
    mat.EM3 <- mat.EM3 + part
  }

  Iobs_hat <- -mat.EM1 / MC.B - mat.EM2 / MC.B + (mat.EM3 / MC.B) %*% t(mat.EM3 / MC.B)
  EM.var   <- solve(Iobs_hat)[1:(p + 1L), 1:(p + 1L), drop = FALSE]

  pvalue <- numeric(p + 1L)
  for (k in seq_len(p + 1L)) {
    pvalue[k] <- 2 * (1 - pnorm(abs(beta_em[k] - beta_new_vec[k]), 0, sqrt(EM.var[k, k])))
  }

  list(p = pvalue, var = diag(EM.var))
}
