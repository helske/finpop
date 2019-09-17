chains <- 8
warmup <- 5000
iter <- 10000
k<-9
n<-203
library(rstan)
library(shinystan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
source("benchmark_data.R")
model <- stan_model("final_model.stan")

births <- readRDS("births_by_parish.rda")
deaths <- readRDS("deaths_by_parish.rda")

births[births==1] <- NA # better to model these, unlikely that there actually was just one death/birth
deaths[deaths==1] <- NA
have_obs_b <- colSums(!is.na(births))>k
have_obs_d <- colSums(!is.na(deaths))>k
births_trimmed <- t(births[, have_obs_b & have_obs_d])
deaths_trimmed <- t(deaths[, have_obs_b & have_obs_d])
(p <- nrow(births_trimmed))
births_obs <- births_trimmed
deaths_obs <- deaths_trimmed
births_obs[is.na(births_obs)] <- 0
deaths_obs[is.na(deaths_obs)] <- 0
tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800,1850, by = 5))
n_obs_pop <- length(tabellverket)
n_miss_pop <- n - n_obs_pop
i_obs_pop <- tabellverket - 1647
i_miss_pop <- (1:n)[!(time(pop_prior) %in% tabellverket)]
data_obs_pop <- pop_prior[time(pop_prior) %in% tabellverket]
data <- tibble::lst(n, p, births_obs, deaths_obs,
                    soldiers, n_miss_pop, n_obs_pop, data_obs_pop, i_miss_pop, i_obs_pop)

b <- births_trimmed
d <- deaths_trimmed
for (i in 1:n) {
  b[is.na(b[, i]), i] <- mean(b[, i], na.rm = TRUE)
  d[is.na(d[, i]), i] <- mean(d[, i], na.rm = TRUE)
}
b[is.nan(b)]<-d[is.nan(d)] <- 6

miss_b <- 0.03 * colMeans(is.na(births_trimmed)) * pop_prior
miss_d <- 0.03 * colMeans(is.na(deaths_trimmed)) * pop_prior
set.seed(123)
inits <- replicate(chains,list(
  ## sample around the posterior mean of pilot runs in order to speed the adaptation
  sigma_drift_b = runif(1, 0.001, 0.01), sigma_drift_d = runif(1, 0.001, 0.01),
  drift_b = runif(n, 0.001, 0.01), drift_d = runif(n, 0, 0.05),
  sigma_rw_b = runif(1, 0.01,0.1), sigma_rw_d = runif(1, 0.01,0.1),
  rate_bd = runif(1, 0.1, 0.4),
  rw_b = log(b)+rnorm(length(b), sd = 0.1), rw_d = log(d)+rnorm(length(d), sd = 0.1),
  miss_b = c(miss_b[apply(is.na(births_trimmed), 2, any)]),
  miss_d = c(miss_d[apply(is.na(deaths_trimmed), 2, any)]),
  mu = pop_prior + 1e4*rnorm(n),
  rate_mu = runif(1,0.1,0.5),
  mu_1647 = pop_prior[1] + 1e4*abs(rnorm(1)), 
  famine = runif(1, 1, 2),
  famine_ratio = runif(1, 0.7, 0.9),
  lc = c(MCMCpack::rdirichlet(1, c(0.5, 0.25, 0.25)*100)),
  lambda_n = runif(1,max(0.75,0.9*p/197),0.95),
  r_b = runif(1, 0.05, 0.15),
  r_d = runif(1, 0.05, 0.15), 
  mid_b_uc = runif(1,0.1,0.9), mid_d_uc = runif(1,0.1,0.9),
  sigma_pop = abs(rnorm(1, 1000, 500))
), simplify=FALSE) 


fit <- sampling(model, data = data, 
                control=list(adapt_delta=0.9, max_treedepth = 11),
                init = inits, 
                pars = c("sigma_drift_b", "sigma_drift_d", "drift_b", "drift_d", 
                         "sigma_rw_b", "sigma_rw_d", "rw_b", "rw_d", "rate_bd", "mu", 
                         "rate_mu", "mu_1647", "sigma_pop", "famine", "famine_ratio", 
                         "r_b", "r_d", "m_b", "m_d", "births", "deaths", "lambda_b", 
                         "lambda_d", "famine_correction", 
                         "mid_b_uc", "mid_d_uc", "lc", "lambda_n"), include = TRUE,
                chains = chains, cores = chains, warmup = warmup, 
                iter = iter, refresh = 10, save_warmup = FALSE)

saveRDS(fit, file="fit_final_16.rds")
check_hmc_diagnostics(fit)
print(fit, "rate_mu")

monitor(extract(fit, pars = c("sigma_drift_b", "sigma_drift_d",
                              "sigma_rw_b", "sigma_rw_d", "rate_bd", 
                              "rate_mu", "mu_1647", "sigma_pop", "famine_ratio", 
                              "r_b", "r_d", "m_b", "m_d", "famine_correction"), 
                permuted = FALSE ),warmup=0)

