
M_step_var = function(dat, t, s, params, E){
  
  
  lambda = update_lambda(dat, t, params$IC$gamma, E$w);cumlambda = c(0,cumsum(lambda))
  
  h = update_h(dat,t,s,E$q,params$RC$alpha, params$RC$beta)
  
  params.new = list();
  params.new$IC = list(); params.new$RC = list()
  params.new$IC$gamma = params$IC$gamma; params.new$IC$lambda = lambda; params.new$IC$cumlambda = cumlambda
  params.new$RC$alpha = params$RC$alpha; params.new$RC$beta = params$RC$beta; params.new$RC$h = h
  
  return(params.new)
}

EM_var = function(dat,t,s,params,k.start,k.end,
                      M,epsilon){

  diff = Inf
  i = 0


  params.old = params
  while(abs(diff) > epsilon && i < M){

    E = E_step(dat,t,s,params.old,k.start,k.end)

    params.new = M_step_var(dat,t,s,params.old,E)

    diff = .Delta(c(params.new$IC$cumlambda, cumsum(params.new$RC$h)),
                  c(params.old$IC$cumlambda, cumsum(params.old$RC$h)))

    i = i+1
    params.old = params.new
  }
  
  if(diff > epsilon){print("nonconvergence!")}
  
  out = list();
  out$E = E; out$params = params.new
  return(out)


}


# Profile likelihood variance for gamma, alpha, and beta.
variance_est = function(dat, t, s, em_fit, k.start, k.end, n,
                        M = 8000, epsilon = 5e-4) {
  hn          <- 2 / sqrt(n)
  p           <- ncol(dat$X)
  total_params <- 2 * p + 1

  log_lik_raw <- log(rowSums(em_fit$E$p))
  score_mat   <- matrix(0, nrow = total_params, ncol = nrow(dat))

  for (k in seq_len(total_params)) {
    params_k <- em_fit$params
    if (k <= p) {
      params_k$IC$gamma[k]    <- params_k$IC$gamma[k] + hn
    } else if (k <= 2 * p) {
      params_k$RC$alpha[k - p] <- params_k$RC$alpha[k - p] + hn
    } else {
      params_k$RC$beta         <- params_k$RC$beta + hn
    }
    var_out       <- EM_var(dat, t, s, params_k, k.start, k.end, M, epsilon)
    E_k           <- E_step(dat, t, s, var_out$params, k.start, k.end)
    score_mat[k, ] <- (log(rowSums(E_k$p)) - log_lik_raw) / hn
  }

  cov_mat <- solve(score_mat %*% t(score_mat))
  return(sqrt(diag(cov_mat)))
}