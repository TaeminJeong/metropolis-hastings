#### GP MH ####
library(MASS)
y <- seq(-4.9,5,0.1) ### generate the data ###
N = length(y)
sigma2_true <- 2
rho_true <- 5
dists <- as.matrix(dist(y))
Sigma_true <- sigma2_true*exp(-dists/rho_true)
x <- mvrnorm(1,rep(0,N),Sigma_true)
# plot(y,x)

## posterior function
log_g_sigma2 <- function(sigma2,x,rho,theta,k,dists,N) {
  Sigma <- sigma2*exp(-dists/rho)
  return(-1/2* t(x) %*% solve(Sigma) %*% x - sigma2/theta + 
           (k-1-N/2)*log(sigma2))
}
log_g_rho <- function(sigma2,x,rho,theta,k,dists,N) {
  Sigma <- sigma2*exp(-dists/rho)
  return(-1/2* t(x) %*% solve(Sigma) %*% x - rho/theta + 
           (k-1)*log(rho)-1/2*log(det(Sigma)))
}


mh_sampl <- function(N, x, n_iter, sigma2_init, rho_init, sigma2_prop_sd,
                     rho_prop_sd, theta,k,dists) {
  sigma2_out <- numeric(n_iter)
  rho_out <- numeric(n_iter)
  accpt_cnt_sigma2 <- 0
  accpt_cnt_rho <- 0
  
  sigma2_now <- sigma2_init
  rho_now <- rho_init
  log_g_sigma2_now <- log_g_sigma2(sigma2=sigma2_now,x=x,rho=rho_now,
                                   theta=theta,k=k,dists=dists,N=N)
  log_g_rho_now <- log_g_rho(sigma2=sigma2_now,x=x,rho=rho_now,
                             theta=theta,k=k,dists=dists,N=N)
  
  for (i in 1:n_iter) {
    ## sigma2 sampling
    sigma2_cand <- rnorm(1, sigma2_now, sigma2_prop_sd)
    while(sigma2_cand <0){
      sigma2_cand <- rnorm(1, sigma2_now, sigma2_prop_sd)
    }
    log_g_sigma2_now <- log_g_sigma2(sigma2=sigma2_now,x=x,rho=rho_now,
                                      theta=theta,k=k,dists=dists,N=N)
    log_g_sigma2_cand <- log_g_sigma2(sigma2=sigma2_cand,x=x,rho=rho_now,
                                      theta=theta,k=k,dists=dists,N=N)
    log_alpha_sigma2 <- log_g_sigma2_cand - log_g_sigma2_now
    alpha <- exp(log_alpha_sigma2)
    
    if (runif(1) < alpha) {
      sigma2_now <- sigma2_cand
      accpt_cnt_sigma2 <- accpt_cnt_sigma2 + 1
    }
    
    sigma2_out[i] <- sigma2_now
    
    ## sigma2 sampling
    rho_cand <- rnorm(1, rho_now, rho_prop_sd)
    while(rho_cand <0){
      rho_cand <- rnorm(1, rho_now, rho_prop_sd)
    }
    log_g_rho_now <- log_g_rho(sigma2=sigma2_now,x=x,rho=rho_now,
                                theta=theta,k=k,dists=dists,N=N)
    log_g_rho_cand <- log_g_rho(sigma2=sigma2_now,x=x,rho=rho_cand,
                                      theta=theta,k=k,dists=dists,N=N)
    log_alpha_rho <- log_g_rho_cand - log_g_rho_now
    alpha <- exp(log_alpha_rho)
    
    if (runif(1) < alpha) {
      rho_now <- rho_cand
      accpt_cnt_rho <- accpt_cnt_rho + 1
    }
    
    rho_out[i] <- rho_now
  }
  
  return(list(sigma2 = sigma2_out,
              rho = rho_out,
              accpt_rate_sigma2 = accpt_cnt_sigma2/n_iter,
              accpt_rate_rho = accpt_cnt_rho/n_iter))
}

## run MCMC
sigma2_init=1
rho_init=1
n_iter=10000
sigma2_prop_sd=3 ## sd of proposal dist'n
rho_prop_sd=6    ## sd of proposal dist'n
theta=5  ## gamma dist'n(prior) scale parameter
k=1      ## gamma dist'n(prior) shape parameter
post_sampl=mh_sampl(N=N, x=x, n_iter=n_iter, sigma2_init=sigma2_init,
         rho_init=rho_init, sigma2_prop_sd=sigma2_prop_sd,
         rho_prop_sd=rho_prop_sd, theta=theta, k=k, dists=dists)

## result
# par(mfrow=c(2,2))
ts.plot(post_sampl$sigma2,main='traceplot',xlab='niter',ylab='sigma2')
ts.plot(post_sampl$rho,main='traceplot',xlab='niter',ylab='rho')
post_sampl$accpt_rate_sigma2
post_sampl$accpt_rate_rho

n_burnin=5000

plot(density(post_sampl$sigma2[(n_iter-n_burnin+1):n_iter]),main='densityplot',xlab='sigma2',ylab='Density') ## nburnin=5000

plot(density(post_sampl$rho[(n_iter-n_burnin+1):n_iter]),main='densityplot',xlab='rho',ylab='Density')

plot(seq(0,7,0.01),dgamma(seq(0,7,0.01),shape=1,scale=5),main='prior_sigma2',xlab='sigma2',ylab='Density')
plot(seq(0,7,0.01),dgamma(seq(0,7,0.01),shape=1,scale=5),main='prior_rho',xlab='rho',ylab='Density')
likelihood=c()
sigma2=seq(0.501,2,0.001)
for(i in 1:length(sigma2)){
  Sigma <- sigma2[i]*exp(-dists/rho_true)
  likelihood[i] <- 1/sqrt(det(Sigma)*(2*pi)^N)*exp(-1/2*t(x)%*%solve(Sigma)%*%x)
}
plot(sigma2,likelihood,main='likelihood_sigma2',xlab='sigma2',ylab='')
likelihood=c()
rho=seq(4,12,0.01)
for(i in 1:length(rho)){
  Sigma <- sigma2_true*exp(-dists/rho[i])
  likelihood[i] <- 1/sqrt(det(Sigma)*(2*pi)^N)*exp(-1/2*t(x)%*%solve(Sigma)%*%x)
}
plot(rho,likelihood,main='likelihood_rho',xlab='rho',ylab='')



### log nimble ####
library(nimble)
y <- c(10,40,30,20,40)
code <- nimbleCode({
  for(i in 1:N){
    y[i] ~ dpois(lam)
  }
  log(lam) ~ dnorm(0,var=1000)
})
N=5
data=list(y=y)
constants = list(N=N)
inits <- list(lam=-1)
examMCMC <- nimbleMCMC(code=code,constants=constants,
                       data=data,inits=inits,
                       niter=10000)
colnames(examMCMC)
plot(examMCMC[2000:10000])
exp(3.3)


#### GP nimble ####
library(nimble)
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    return(result)
  })
cExpcov <- compileNimble(expcov)
N=5
y <- c(1,2,3,6,10)
y <- y-mean(y)
code <- nimbleCode({
  rho ~ dgamma(shape=1,scale=5)
  sigma2 ~ dgamma(shape=1,scale=5)
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma2)
  y[1:N] ~ dmnorm(zeros[1:N], cov = cov[1:N, 1:N])
})
data=list(y=y)
constants = list(N=N,dists=as.matrix(dist(1:N)),
                 zeros=rep(0,N))
inits <- list(rho=1,sigma2=1)
examMCMC <- nimbleMCMC(code=code,constants=constants,
                       data=data,inits=inits,
                       niter=10000)
colnames(examMCMC)
ts.plot(examMCMC[,1])
ts.plot(examMCMC[,2])
rho=mean(examMCMC[,1])
sigma2=mean(examMCMC[,2])
k_tstar_tstar <- 
  cExpcov(dists = as.matrix(dist(1:N)),rho = rho,sigma2 = sigma2)
k_t_tstar <- sigma2*exp(-seq(5,1)/rho)
k_tstar_t <- sigma2*exp(-seq(1,5)/rho)
k_tt <- sigma2
mean <- t(k_t_tstar) %*% solve(k_tstar_tstar) %*% y
var <- k_tt - t(k_t_tstar) %*% solve(k_tstar_tstar) %*% k_tstar_t
sqrt(var)
plot(y)




