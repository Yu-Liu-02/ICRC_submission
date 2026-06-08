draw_X = function(n,p){
  
  x1 <- rbinom(n,1,0.6)
  x2 <- truncnorm(n, 0, 1, -2, 2)
  x3 <- truncnorm(n, 0, 1, -2, 2)
  cbind(x1, x2, x3)
}



draw_TT = function(n, X, gamma){
  
  lp <- as.vector(X %*% gamma)
  E  <- rexp(n)
  3 * (E / exp(lp))^(1 / 2)
  
}



draw_Z = function(n, X, TT, alpha, beta, g_type = "indicator"){
  
  E   <- rexp(n)
  eAX <- as.vector(exp(X %*% alpha))
  HT  <- TT^2 * eAX / 80   # cumulative hazard at TT
  
  if (g_type == "indicator") {
    return(ifelse(E < HT,
                  sqrt(80 * E / eAX),
                  sqrt(TT^2 + 80 * (E - HT) / (eAX * exp(beta)))))
  }
  
  # ReLU: g(z,TT) = max(z-TT, 0)
  H_relu <- function(z, TT_i, eAX_i) {
    w <- z - TT_i
    eAX_i * (TT_i^2 / 80 +
               (1/40) * (exp(beta * w) * (z / beta - 1/beta^2) + 1/beta^2 - TT_i/beta))
  }
  
  Z <- numeric(n)
  for (i in 1:n) {
    if (E[i] < HT[i]) {
      Z[i] <- sqrt(80 * E[i] / eAX[i])
    } else {
      f  <- function(z) H_relu(z, TT[i], eAX[i]) - E[i]
      hi <- TT[i] + 1
      while (f(hi) < 0) hi <- hi * 2
      Z[i] <- uniroot(f, lower = TT[i], upper = hi, tol = 1e-10)$root
    }
  }
  Z
}




draw_examination = function(n, count, X, psi) {
  lp  <- as.vector(X %*% psi)          
  mean_gap <- 8 / exp(lp)                 
  
  E <- matrix(0, n, count - 1)
  for (i in 1:n) {
    E[i, ] <- rexp(count - 1, rate = 1 / mean_gap[i])
  }
  
  U      <- matrix(0, n, count)
  U[, 1] <- runif(n, 0, 4)
  for (j in 2:count) {
    U[, j] <- U[, j-1] + 0.5 + E[, j-1]
  }
  return(U)
}



get_LR = function(TT, U){
  
  is.prev = U < TT
  t.prev = U
  t.prev[!is.prev] = -Inf
  L = do.call(pmax, as.data.frame(t.prev))
  L[is.infinite(L)] = 0
  
  is.after = U >= TT
  t.after = U
  t.after[!is.after] = Inf
  R = do.call(pmin, as.data.frame(t.after))
  
  return(data.frame(L = L, R = R))
  
}




# draw preclinical event from a bump hazard
Lambda0_bump <- function(t, c, a, omega, b, tau, sigma) {
  base <- c * (t - (a / omega) * sin(omega * t))
  bump <- b * sqrt(2*pi) * sigma *
    (pnorm((t - tau) / sigma) - pnorm((0 - tau) / sigma))
  base + bump
}

rT_PH_bump <- function(n, X, gamma,
                       c=0.10, a=0.98, omega=pi/3,
                       b=0.80, tau=1.5, sigma=0.05,
                       t_upper=1e4) {
  stopifnot(nrow(X) == n, ncol(X) == length(gamma))
  stopifnot(a > 0 && a < 1, c > 0, omega > 0, b >= 0, sigma > 0)
  
  eta <- as.vector(X %*% gamma)
  U <- runif(n)
  target <- -log(U) / exp(eta)
  
  vapply(seq_len(n), function(i) {
    f <- function(t) Lambda0_bump(t,c,a,omega,b,tau,sigma) - target[i]
    hi <- t_upper
    while (f(hi) < 0) {
      hi <- hi * 2
      if (hi > 1e8) stop("Failed to bracket root.")
    }
    uniroot(f, lower=0, upper=hi, tol=1e-10)$root
  }, numeric(1))
}


