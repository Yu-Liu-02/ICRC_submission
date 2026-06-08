##-------------------------------------------##
library(survival)
library(bayess)
source("../../data_generation/data_generate.R")
source("../../data_generation/g_function.R")
##-------------------------------------------##



##-------------------------------------------##
n = 1000
Num_INSTANCES = 1000
INSTANCES_PER_JOB = 10
p=3
set_g_type("relu")
gamma <- c(1, 1.5, -1.5)
alpha <- c(0.45,0.5,-0.25)
psi = c(0.4,0.4,-0.2)
beta  <- 0.1
c0     <- 0.03
a0     <- 0.99
omega0 <- 2*pi/2
b0     <- 7
tau0   <- 2.5
sigma0 <- 0.12
##-------------------------------------------##




##-------------------------------------------##
# Midpoint imputation ##
coef = matrix(0,nrow =INSTANCES_PER_JOB*Num_INSTANCES, ncol = p+1)
htsis = matrix(0,nrow =INSTANCES_PER_JOB*Num_INSTANCES, ncol = 2*p+2)
for(i in 1:Num_INSTANCES){
  set.seed(i)
  for(j in 1:INSTANCES_PER_JOB){
    
    
    X = draw_X(n,p)
    TT <- rT_PH_bump(n, X, gamma,c0,a0,omega0,b0,tau0,sigma0)
    Z = draw_Z(n, X, TT, alpha, beta, "relu")
    C = rexp(n, 1/6)
    Y = ifelse(Z<=C, Z, C)
    delta = ifelse(Z<=C, 1, 0)
    U = draw_examination(n, 3, X, psi)
    LR = get_LR(TT, U)
    
    time = Y
    status = delta
    L = LR$L; R = LR$R
    midpoints <- ifelse(is.infinite(R), R, (L + R) / 2)
    df_base <- data.frame(id = 1:n, time = time, status = status, X1 = X[,1], X2 = X[,2], X3 = X[,3], mid = midpoints)
    fit <- coxph(
      Surv(time, status) ~ X1 + X2 + X3 + tt(mid),
      data = df_base,
      tt = function(x, t, ...) g(t, x)
    )
    coef[(i-1)*INSTANCES_PER_JOB+j,] = fit$coefficients
    htsis[(i-1)*INSTANCES_PER_JOB+j,1:4] = c(2*(1-pnorm(abs(fit$coefficients-c(alpha,beta))[1],0,fit$var[1,1]^0.5)),2*(1-pnorm(abs(fit$coefficients-c(alpha,beta))[2],0,fit$var[2,2]^0.5)),
                                             2*(1-pnorm(abs(fit$coefficients-c(alpha,beta))[3],0,fit$var[3,3]^0.5)),2*(1-pnorm(abs(fit$coefficients-c(alpha,beta))[4],0,fit$var[4,4]^0.5)))
    htsis[(i-1)*INSTANCES_PER_JOB+j,5:8] = sqrt(diag(fit$var))
  }
}

true_ab     <- c(alpha, beta)
param_names <- c(paste0("alpha[", 1:p, "]"), "beta")

result <- rbind(
  Bias = colMeans(coef) - true_ab,
  SE   = apply(coef, 2, sd),
  SEE  = colMeans(htsis[, (p+2):(2*p+2)]),
  CP   = colMeans(htsis[, 1:(p+1)] > 0.05)
)
colnames(result) <- param_names
print(round(result, 4))
##-------------------------------------------##


