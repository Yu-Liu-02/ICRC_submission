##-------------------------------------------##
library(survival)
library(bayess)
source("../../data_generation/data_generate.R")
source("../../Goggins_func/compile_cpp_mcem.R")
source("../../Goggins_func/MCEM_Goggins.R")
##-------------------------------------------##

##-------------------------------------------##
library(survival)
library(bayess)
n = 1000
Num_INSTANCES = 1000
INSTANCES_PER_JOB = 10
c0     <- 0.03
a0     <- 0.99
omega0 <- pi
b0    <- 7
tau0   <- 2.5
sigma0 <- 0.12
gamma <- c(1,1.5,-1.5)
alpha <- c(0.45,5,-0.25)
beta  <- 1
psi <- c(0.4,0.4,-0.2)
p=3
##-------------------------------------------##







pe=matrix(0,INSTANCES_PER_JOB*Num_INSTANCES,p+1)
se=matrix(0,INSTANCES_PER_JOB*Num_INSTANCES,p+1)
for(iter in 1:Num_INSTANCES){
  set.seed(iter)
  for(iter_in in 1:INSTANCES_PER_JOB){
    X = draw_X(n,p)
    TT <- rT_PH_bump(n, X, gamma,c0,a0,omega0,b0,tau0,sigma0)
    U = draw_examination(n,3,X,psi)
    LR = get_LR(TT, U)
    Z = draw_Z(n, X, TT, alpha, beta)
    C = rexp(n, 1/6)
    Y = ifelse(Z<=C, Z, C)
    delta = ifelse(Z<=C, 1, 0)
  
    time = Y
    status = delta
    L = LR$L; R = LR$R
    midpoints <- ifelse(is.infinite(R), R, (L + R) / 2)
    df_base <- data.frame(id = 1:n, time = time, status = status, X1 = X[,1], X2 = X[,2], X3 = X[,3], mid = midpoints)
    df_tv <- tmerge(data1 = df_base, data2 = df_base, id = id,
                    event = event(time, status),
                    Z_tv = tdc(mid))
    fit_mid <- coxph(Surv(tstart, tstop, event) ~ X1+ X2+ X3 + I(Z_tv), data = df_tv)
    start_alpha <- unname(coef(fit_mid)[1:p])
    start_beta  <- unname(coef(fit_mid)[p+1])
    out = fit_mcem_interval(time, status, X, L, R, start = c(alpha = start_alpha, beta = start_beta))
    pe[(iter-1)*INSTANCES_PER_JOB+iter_in,] = out$coef
    se[(iter-1)*INSTANCES_PER_JOB+iter_in,] = out$se
  }
}

true_ab     <- c(alpha, beta)
param_names <- c(paste0("alpha[", 1:p, "]"), "beta")

result <- rbind(
  Bias = colMeans(pe) - true_ab,
  SE   = apply(pe, 2, sd),
  SEE  = colMeans(se),
  CP   = colMeans(abs(sweep(pe, 2, true_ab)) <= 1.96 * se)
)
colnames(result) <- param_names
print(round(result, 4))
