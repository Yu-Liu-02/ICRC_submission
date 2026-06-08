##-------------------------------------------##
library(Rcpp)
library(RcppArmadillo)
library(bayess)
source("../../data_generation/data_generate.R")
source("../../data_generation/g_function.R")
source("../../Proposed_func/compile_cpp_em.R")
source("../../Proposed_func/EM_proposed.R")
source("../../Proposed_func/variance_estimation.R")
##-------------------------------------------##


##-------------------------------------------##
n = 1000 # sample size
Num_INSTANCES = 1000
INSTANCES_PER_JOB = 10
set_g_type("indicator")
c     <- 0.03
a     <- 0.99
omega <- 2*pi/2
b     <- 7
tau   <- 2.5
sigma <- 0.12
gamma <- c(1,1.5,-1.5)
alpha <- c(0.9,1,-0.5)/2
psi <- c(0.4,0.4,-0.2)
beta  <- 1
p = 3
##-------------------------------------------##





##-------------------------------------------##
# EM algorithm ###
pe=matrix(0,INSTANCES_PER_JOB*Num_INSTANCES,2*p+1)
se=matrix(0,INSTANCES_PER_JOB*Num_INSTANCES,2*p+1)
for(iter in 1:Num_INSTANCES){
  set.seed(iter)
  for(iter_in in 1:INSTANCES_PER_JOB){
    # Point estimate
    X = draw_X(n,p)
    TT = rT_PH_bump(n, X, gamma, c, a, omega, b, tau, sigma)
    U = draw_examination(n,3,X,psi)
    LR = get_LR(TT, U)
    Z = draw_Z(n, X, TT, alpha, beta, g_type = "indicator")
    C = rexp(n, 1/6)
    Y = ifelse(Z<=C, Z, C)
    delta = ifelse(Z<=C, 1, 0)
    
    
    t = sort(unique(unlist(LR)))
    t = t[c(-1,-length(t))]
    s = sort(unique(Y[delta==1]))
    
    dat = data.frame(Y = Y, delta = delta, L = LR$L, R = LR$R)
    dat$X = X
    k.start = findInterval(dat$L, t)+1
    k.end = findInterval(dat$R, t)
    
    m1 = length(t); m2 = length(s)
    lambda = rep(1/m1,m1); h = rep(1/m2,m2)
    cumlambda = c(0,cumsum(lambda))
    params = list()
    params$IC = list(); params$RC = list()
    params$IC$gamma = rep(0,p); params$IC$lambda = lambda; params$IC$cumlambda = cumlambda
    params$RC$alpha = rep(0,p); params$RC$beta = 0; params$RC$h = h
    
    EM.fit = EM_proposed(dat,t,s,params,k.start,k.end,8000,5e-4)
    pe[(iter-1)*INSTANCES_PER_JOB+iter_in,] = c(EM.fit$params$IC$gamma, EM.fit$params$RC$alpha, EM.fit$params$RC$beta)
    # Variance estimation
    se[(iter-1)*INSTANCES_PER_JOB+iter_in,] = variance_est(dat, t, s, EM.fit, k.start, k.end, n,
                                                      8000,5e-4)
  }
}

true_ab     <- c(alpha, beta)
pe_ab       <- pe[, (p+1):(2*p+1), drop = FALSE]
se_ab       <- se[, (p+1):(2*p+1), drop = FALSE]
param_names <- c(paste0("alpha[", 1:p, "]"), "beta")

result <- rbind(
  Bias = colMeans(pe_ab) - true_ab,
  SE   = apply(pe_ab, 2, sd),
  SEE  = colMeans(se_ab),
  CP   = colMeans(abs(sweep(pe_ab, 2, true_ab)) <= 1.96 * se_ab)
)
colnames(result) <- param_names
print(round(result, 4))
##-------------------------------------------##

