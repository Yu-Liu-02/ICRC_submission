##-------------------------------------------##
library(survival)
library(bayess)
n = 1000 # sample size
source("../../data_generation/data_generate.R")
source("../../Ahn_func/EM_Ahn.R")
##-------------------------------------------##

Num_INSTANCES = 1000
INSTANCES_PER_JOB = 10
p = 3
true_alpha  <- c(0.45, 0.5, -0.25)
true_beta   <- 1
pointestimate=matrix(0,INSTANCES_PER_JOB*Num_INSTANCES,p+1)
htsis = matrix(0,INSTANCES_PER_JOB*Num_INSTANCES,2*(p+1))
for(iter in 1:Num_INSTANCES){
  set.seed(iter)
  for(iter_in in 1:INSTANCES_PER_JOB){
    gamma <- c(1,1.5,-1.5)
    alpha <- c(0.45,0.5,-0.25)
    psi <- c(0.4,0.4,-0.2)
    beta  <- 1
    c     <- 0.03
    a     <- 0.99
    omega <- 2*pi/2
    b     <- 7
    tau   <- 2.5
    sigma <- 0.12
    
    X = draw_X(n,p)
    TT = rT_PH_bump(n, X, gamma, c, a, omega, b, tau, sigma)
    U = draw_examination(n,3,X,psi)
    LR = get_LR(TT, U)
    Z = draw_Z(n, X, TT, alpha, beta)
    C = rexp(n, 1/6)
    Y = ifelse(Z<=C, Z, C)
    delta = ifelse(Z<=C, 1, 0)
    t = sort(unique(unlist(LR)))
    t = t[c(-1,-length(t))]
    s = sort(unique(Y[delta==1]))
    
    
    dat = data.frame(Y = Y, delta = delta, L = LR$L, R = LR$R)
    dat = cbind(dat, as.data.frame(X))   
    
    xt.data = dat
    xt.data$id = 1:n
    names(xt.data) = c("T","delta","u","v", paste0("z",1:p), "id")
    xt.data = xt.data[, c(p+5, 1, 2, 3, 4, 5:(p+4))]   
    xt.data$u[xt.data$u==0] = 0.00001
    xt.data$v[xt.data$v==Inf] = 999
    head(xt.data)
    xt.data <- xt.data[order(xt.data[, 2]), ]
    N=dim(xt.data)[1]
    tt=xt.data[,2]
    delta=xt.data[,3]
    x.left=xt.data[,4]
    x.right=xt.data[,5]
    z=as.matrix(xt.data[,6:(5+p)])      
    colnames(xt.data)[6:(5+p)] = paste0("z", 1:p)
    
    # R: the observation indicator over time (N by N matrix)
    # Y: at risk indicator
    # D: counting process for the secondary event
    R=Y=D=matrix(0,N,N) 
    for(i in (1:N)){ 
      Y[i,]=as.numeric(tt<=tt[i])
      D[i,]=as.numeric(tt>=x.right[i])
      R[i,]=as.numeric(tt<=x.left[i])+as.numeric(tt>=x.right[i]) 
    }
    
    beta=beta.new=matrix(c(0.45,0.5,-0.25,1),p+1,1)
    beta.new.EM=rep(0,p+1)
    maxite = 1000
    Dout.EM=rbind(D,exp(D*beta.new.EM[p+1]),D*exp(D*beta.new.EM[p+1]))  
    
    fit=survreg(as.formula(paste0("Surv(u,v,type='interval2')~",paste0("z",1:p,collapse="+"))),data=xt.data,dist="weibull")
    alpha=1/fit$scale
    gamma=-fit$coef[2:(p+1)]*alpha
    h0=-fit$coef[1]*alpha
    for(mite in 1:maxite){
      beta.EM=beta.new.EM 
      cumlambda0.hat.EM=lambda0.EM(beta.EM,Dout.EM)  
      Dout.EM=expD.EM(beta.EM,cumlambda0.hat.EM)
      Sout.EM=Sobs.EM(beta.EM,xt.data,Dout.EM)
      score.hat.EM=score.EM(beta.EM,Sout.EM,Dout.EM,cumlambda0.hat.EM)
      update.EM=solve(hessian.EM(beta.EM,Sout.EM,Dout.EM,cumlambda0.hat.EM))%*%score.hat.EM
      beta.new.EM=beta.EM+update.EM
      if(max(abs(update.EM))<0.00005) break;
    }
    cat("EM", beta.new.EM, "\n")
    pointestimate[(iter-1)*INSTANCES_PER_JOB+iter_in,] = beta.new.EM
    cov.prob.out = tryCatch(covProbf(),
                            error = function(e) {
                              message("covProbf failed: ", conditionMessage(e))
                              list(p = rep(NA, p+1), var = rep(NA, p+1))
                            })
    htsis[(iter-1)*INSTANCES_PER_JOB+iter_in,] = c(cov.prob.out$p,cov.prob.out$var)
  }
}

true_ab     <- c(true_alpha, true_beta)
param_names <- c(paste0("alpha[", 1:p, "]"), "beta")

result <- rbind(
  Bias = colMeans(pointestimate) - true_ab,
  SE   = apply(pointestimate, 2, sd),
  SEE  = sqrt(colMeans(htsis[, (p+2):(2*(p+1))])),
  CP   = colMeans(htsis[, 1:(p+1)])
)
colnames(result) <- param_names
print(round(result, 4))
