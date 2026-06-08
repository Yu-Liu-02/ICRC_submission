
# E-step


calculate_p <- function(dat, t, s, params, k.start, k.end) {
    calculate_p_cpp(
      X         = dat$X,
      Y         = dat$Y,
      delta     = as.integer(dat$delta),
      R         = dat$R,
      gamma     = params$IC$gamma,
      alpha     = params$RC$alpha,
      beta      = as.double(params$RC$beta),
      cumlambda = params$IC$cumlambda,
      h         = params$RC$h,
      t         = t,
      s         = s,
      k_start   = as.integer(k.start),
      k_end     = as.integer(k.end),
      g_type    = g_type_int()
    )
}

calculate_w <- function(dat, t, params, q, k.start, k.end) {
    calculate_w_cpp(
      X       = dat$X,
      gamma   = params$IC$gamma,
      lambda  = params$IC$lambda,
      q       = q,
      k_start = as.integer(k.start),
      k_end   = as.integer(k.end)
    )
}





calculate_q = function(p){
  
  return(sweep(p,1,rowSums(p),"/"))
  
}







E_step = function(dat, t, s, params, k.start, k.end){
  
  
  p = calculate_p(dat,t,s,params,k.start,k.end)
  
  q = calculate_q(p)
  
  w = calculate_w(dat,t,params,q,k.start,k.end)
  
  return(list(p=p, q = q, w = w))
  
  
}

obs_loglik <- function(p){
  rs = rowSums(p)
  if (any(rs <= 0 | !is.finite(rs))) return(-Inf)  # guard
  sum(log(rs))
}



# M-step

update_gamma_R = function(dat, t, w, gamma){
  
  
  p = ncol(dat$X); m1 = length(t)
  
  score = rep(0,p); hessian = matrix(0,p,p)
  
  for(k in 1:m1){
    
    at.risk = which(dat$R >= t[k])
    
    w.risk = w[at.risk,k]
    X.risk = dat$X[at.risk,,drop = F]
    
    e2.gammaX = as.numeric(exp(X.risk%*%gamma))
    
    
    s0 = sum(e2.gammaX); s1 = colSums(X.risk*e2.gammaX)
    temp = s1/s0
    centered_X = sweep(X.risk, 2, temp, "-") 
    score = score + colSums(centered_X * w.risk)
    
    s2 = t(X.risk)%*%(X.risk*e2.gammaX)
    V = s2/s0 - tcrossprod(temp)
    hessian = hessian - sum(w.risk)*V
    
  }
  
  
  return(as.vector(gamma-solve(hessian)%*%score))
  
}

update_gamma = function(dat, t, w, gamma){
    update_gamma_cpp(X = dat$X, R = dat$R, t = t, w = w, gamma = gamma)
}




update_lambda = function(dat, t, gamma, w){
  
  n = nrow(dat)
  
  e2.gammaX = as.numeric(exp(dat$X%*%gamma))
  risk.set = data.frame(id = 1:n, R = dat$R, e2.gammaX = e2.gammaX)
  risk.set = risk.set[order(risk.set$R),]
  
  
  at.risk = rep(T,n)
  m1 = length(t)
  denom.now = sum(e2.gammaX)
  
  lambda = rep(0, m1)
  pt = 1
  for(k in 1:m1){
    
    tk = t[k]
    while(pt <= n && risk.set$R[pt] < tk){
      drop = risk.set$id[pt]
      
      denom.now = denom.now - risk.set$e2.gammaX[pt]
      
      at.risk[drop] = F
      
      pt = pt+1
    }
    numer = sum(w[at.risk,k])
    lambda[k] = numer/denom.now
  }
  
  return(lambda)
}



update_h <- function(dat, t, s, q, alpha, beta) {
    update_h_cpp(
      X      = dat$X,
      Y      = dat$Y,
      delta  = as.integer(dat$delta),
      s      = s,
      t      = t,
      q      = q,
      alpha  = alpha,
      beta   = as.double(beta),
      g_type = g_type_int()
    )
}


update_alphabeta <- function(dat, t, s, q, alpha, beta) {
    update_alphabeta_cpp(
      X      = dat$X,
      Y      = dat$Y,
      delta  = as.integer(dat$delta),
      t      = t,
      s      = s,
      q      = q,
      alpha  = alpha,
      beta   = as.double(beta),
      g_type = g_type_int()
    )
}




M_step = function(dat, t, s, params, E){
  
  gamma = update_gamma(dat, t, E$w, params$IC$gamma)
  
  lambda = update_lambda(dat, t, gamma, E$w);cumlambda = c(0,cumsum(lambda))
  
  theta = update_alphabeta(dat, t, s, E$q, params$RC$alpha, params$RC$beta)
  
  h = update_h(dat,t,s,E$q,theta$alpha, theta$beta)
  
  params.new = list();
  params.new$IC = list(); params.new$RC = list()
  params.new$IC$gamma = gamma; params.new$IC$lambda = lambda; params.new$IC$cumlambda = cumlambda
  params.new$RC$alpha = theta$alpha; params.new$RC$beta = theta$beta; params.new$RC$h = h
  
  return(params.new)
}










# EM

# Convergence metric: relative difference for large values, absolute near zero.

.Delta <- function(v, w, delta = 0.01, eps = 1e-6) {
  use_abs <- abs(v) <= delta | abs(w) <= delta | pmin(abs(v), abs(w)) <= eps
  diffs   <- ifelse(use_abs, abs(v - w), abs(v - w) / pmin(abs(v), abs(w)))
  max(diffs)
}




EM_proposed = function(dat,t,s,params,k.start,k.end,
                       M,epsilon){
  
  diff = Inf
  i = 0
  
  
  params.old = params
  while(diff > epsilon && i < M){
    
    E = E_step(dat,t,s,params.old,k.start,k.end)
    
    params.new = M_step(dat,t,s,params.old,E)
    
    diff = .Delta(c(params.new$IC$gamma, params.new$RC$alpha, params.new$RC$beta),
                  c(params.old$IC$gamma, params.old$RC$alpha, params.old$RC$beta)) +
      .Delta(c(params.new$IC$cumlambda, cumsum(params.new$RC$h)),
             c(params.old$IC$cumlambda, cumsum(params.old$RC$h)))
    
    i = i+1
    params.old = params.new
  }
  
  if(diff > epsilon){print("nonconvergence!")}
  
  E = E_step(dat, t, s, params.new, k.start, k.end)
  
  out = list();
  out$E = E; out$params = params.new
  return(out)
  
  
}
