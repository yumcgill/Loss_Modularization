library(rstan)
library(gtools)      
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

a_true   <- 1.0
b_true   <- 2.0
sigma    <- 1.0
sigma_x  <- 3
tau_z    <- 3.0
N        <- 100

n_sim    <- 1000
n_chains <- 2
n_iter   <- 3000
n_warmup <- 1000
N_SYN_SETS <- 1000     
alpha_dp   <- 2      
sig_dp     <- 0.5    

### Stan functions for Full Bayes and Cut Posteriors
stan_full <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] x;
  real<lower=0> sigma_x;
  real<lower=0> tau_z;
}
parameters {
  real a;
  real b;
  real<lower=0> sigma;
  vector[N] z;
}
model {
  a     ~ normal(0, 10);
  b     ~ normal(0, 10);
  sigma ~ normal(0, 5);
  z     ~ normal(0, tau_z);
  x     ~ normal(z, sigma_x);
  y     ~ normal(a + b * z, sigma);
}
"

stan_cut1 <- "
data {
  int<lower=1> N;
  vector[N] x;
  real<lower=0> sigma_x;
  real<lower=0> tau_z;
}
parameters {
  vector[N] z_star;
}
model {
  z_star ~ normal(0, tau_z);
  x      ~ normal(z_star, sigma_x);
}
"

stan_cut2 <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] z_fixed;
}
parameters {
  real a;
  real b;
  real<lower=0> sigma;
}
model {
  a     ~ normal(0, 10);
  b     ~ normal(0, 10);
  sigma ~ normal(0, 5);
  y     ~ normal(a + b * z_fixed, sigma);
}
"
extract_mu <- function(fit, N) {
  mu_draws <- rstan::extract(fit, pars = "mu")$mu
  list(
    mean = colMeans(mu_draws),
    lo   = apply(mu_draws, 2, quantile, 0.025),
    hi   = apply(mu_draws, 2, quantile, 0.975)
  )
}

extract_scalar <- function(fit, par) {
  mean(rstan::extract(fit, pars = par)[[par]])
}


#### DP-cut functions
#### Resample function using Polya urn
polya_resample <- function(n) {
  
  idx <- integer(n)
  counts <- rep(1, n)
  for (j in 1:n) {
    
    if (j == 1) {
      idx[j] <- sample.int(n, size = 1)
    } else {
      idx[j] <- sample.int(n,size = 1,prob = counts)
    }
    counts[idx[j]] <- counts[idx[j]] + 1
  }
  idx
}


dp_cut <- function(dat, K, alpha = 1, sig = 0.5) {
  
  x <- dat$x
  y <- dat$y
  n <- length(x)
  
  lm_fit  <- lm(y ~ x)
  a_hat   <- coef(lm_fit)[1]
  b_hat   <- coef(lm_fit)[2]
  
  a_draws <- NULL
  b_draws <- NULL
  
  for (k in seq_len(K)) {
    
    idx <- polya_resample(n)
    x_star <- x[idx]
    y_star  <- numeric(n)
    n_seen  <- n
    
    for (i in seq_len(n)) {
      prob_new <- alpha / (alpha + n_seen)   
      
      if (n_seen == 0 || runif(1) < prob_new) {
        y_star[i] <- rnorm(1, a_hat + b_hat * x_star[i], sig)
      } else {
        y_star[i] <- y[idx[i]]
      }
      n_seen <- n_seen + 1
    }
    ### Estimation for Module 1
    z_star_k <- (tau_z^2 / (tau_z^2 + sigma_x^2)) * x_star
    
    ### Estimation for Module 2
    ols_k   <- lm(y_star ~ z_star_k)
    a_draws <- c(a_draws, coef(ols_k)[1])
    b_draws <- c(b_draws, coef(ols_k)[2])
  }
  list(a = a_draws, b = b_draws)
}


mod_full <- stan_model(model_code = stan_full,  model_name = "full_bayes")
mod_cut1 <- stan_model(model_code = stan_cut1,  model_name = "cut_mod1")
mod_cut2 <- stan_model(model_code = stan_cut2,  model_name = "cut_mod2")


results <- vector("list", n_sim)
for (s in seq_len(n_sim)) {
  # Generate data
  z_true  <- rnorm(N, 0, tau_z)
  x       <- rnorm(N, z_true, sigma_x)
  mu_true <- a_true + b_true * z_true
  y       <- rnorm(N, mu_true, sigma)
  dat     <- list(y = y, x = x, z_true = z_true, mu_true = mu_true)
  
  dp <- dp_cut(dat, N_SYN_SETS, alpha = alpha_dp, sig = sig_dp)
  
  a_dp         <- mean(dp$a)
  b_dp         <- mean(dp$b)
  
  #Full Bayes 
  fit_full <- sampling(
    mod_full,
    data    = list(N = N, y = y, x = x,
                   sigma_x = sigma_x, tau_z = tau_z),
    chains  = n_chains, iter = n_iter, warmup = n_warmup,
    refresh = 0
  )
  a_full         <- extract_scalar(fit_full, "a")
  b_full         <- extract_scalar(fit_full, "b")
  
  #Cut posterior
  #Module 1: z* from x only
  fit_cut1 <- sampling(
    mod_cut1,
    data    = list(N = N, x = x, sigma_x = sigma_x, tau_z = tau_z),
    chains  = n_chains, iter = n_iter, warmup = n_warmup,
    refresh = 0
  )
  z_star_draws <- rstan::extract(fit_cut1, pars = "z_star")$z_star
  z_fixed      <- colMeans(z_star_draws)
  
  # Module 2: outcome with z fixed
  fit_cut2 <- sampling(
    mod_cut2,
    data    = list(N = N, y = y, z_fixed = z_fixed),
    chains  = n_chains, iter = n_iter, warmup = n_warmup,
    refresh = 0
  )
  a_cut         <- extract_scalar(fit_cut2, "a")
  b_cut         <- extract_scalar(fit_cut2, "b")
  
  results[[s]] <- data.frame(a_dp = a_dp, b_dp = b_dp, a_full = a_full, b_full = b_full, a_cut = a_cut,b_cut = b_cut)
}

df <- bind_rows(results)


