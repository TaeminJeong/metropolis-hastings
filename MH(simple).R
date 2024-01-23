log_g <- function(mu, n, x_bar) {
  mu2 <- mu^2
  
  return(n * (x_bar * mu - mu2 / 2.0) - log(1.0 + mu2))
}

mh_sampl <- function(n_data, x_bar, n_iter, mu_init, prop_sd) {
  mu_out <- numeric(n_iter)
  accpt_cnt <- 0
  
  mu_now <- mu_init
  log_g_now <- log_g(mu_now, n_data, x_bar)
  
  for (i in 1:n_iter) {
    mu_cand <- rnorm(1, mu_now, prop_sd)
    log_cand <- log_g(mu_cand, n_data, x_bar)
    log_alpha <- log_cand - log_g_now
    alpha <- exp(log_alpha)
    
    if (runif(1) < alpha) {
      mu_now <- mu_cand
      accpt_cnt <- accpt_cnt + 1
      log_g_now <- log_cand
    }
    
    mu_out[i] <- mu_now
  }
  
  return(list(mu = mu_out, accpt_rate = accpt_cnt/n_iter))
}

x <- c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)

x_bar <- mean(x)
n_data <- length(x)

library(coda)

post_sampl <- mh_sampl(n_data = n_data, x_bar = x_bar, n_iter = 1e4, mu_init = 0.0, prop_sd = 0.1)

traceplot(as.mcmc(post_sampl$mu))
post_sampl$accpt_rate

