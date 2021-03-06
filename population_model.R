chains <- 16
warmup <- 2500
iter <- 12500
k<-9
n<-203
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
source("benchmark_data.R")
model <- stan_model("population_model.stan")

births <- readRDS("births_by_parish.rda")
deaths <- readRDS("deaths_by_parish.rda")

births[births==1] <- NA # better to model these, unlikely that there actually was just one death/birth
deaths[deaths==1] <- NA
have_obs_b <- colSums(!is.na(births))>k
have_obs_d <- colSums(!is.na(deaths))>k
births_trimmed <- t(births[, have_obs_b & have_obs_d])
deaths_trimmed <- t(deaths[, have_obs_b & have_obs_d])
(p <- nrow(births_trimmed))

baptisms <- births_trimmed
burials <- deaths_trimmed
baptisms[is.na(baptisms)] <- 0
burials[is.na(burials)] <- 0
tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800,1850, by = 5))
n_obs_census <- length(tabellverket)
n_miss_census <- n - n_obs_census
i_obs_census <- tabellverket - 1647
i_miss_census <- (1:n)[!(time(pop_prior) %in% tabellverket)]
census <- pop_prior[time(pop_prior) %in% tabellverket]
data <- tibble::lst(n, p, baptisms, burials,
                    soldiers, n_miss_census, n_obs_census, 
                    census, i_miss_census, i_obs_census)

b <- births_trimmed
d <- deaths_trimmed
for (i in 1:n) {
  b[is.na(b[, i]), i] <- mean(b[, i], na.rm = TRUE)
  d[is.na(d[, i]), i] <- mean(d[, i], na.rm = TRUE)
}


miss_b <- 0.03 * colMeans(is.na(births_trimmed)) * pop_prior
miss_d <- 0.025 * colMeans(is.na(deaths_trimmed)) * pop_prior
set.seed(123)
inits <- replicate(chains, list(
  ## sample around the posterior means of pilot runs in order to speed the adaptation
  sigma_eta_b = runif(1, 0.1, 0.3), sigma_eta_d = runif(1, 0.1, 0.3),
  eta_b = runif(n, 0.001, 0.01), eta_d = runif(n, 0, 0.05),
  sigma_nu_b = runif(1, 0.1,0.4), sigma_nu_d = runif(1, 0.1,0.4),
  psi_bd = runif(1, 0.2, 0.4),
  nu_b = log(b)+rnorm(length(b), sd = 0.1), nu_d = log(d)+rnorm(length(d), sd = 0.1),
  miss_b = c(miss_b[apply(is.na(births_trimmed), 2, any)]),
  miss_d = c(miss_d[apply(is.na(deaths_trimmed), 2, any)]),
  mu = pop_prior + 1e4*rnorm(n),
  psi_mu = runif(1, 0.25, 1),
  mu_1647 = pop_prior[1] + 1e4*abs(rnorm(1)), 
  log_phi_d = rnorm(1, 2, 0.25),
  pi = rbeta(1, 0.775 * 200, (1 - 0.775) * 200),
  lc = c(MCMCpack::rdirichlet(1, c(10, 5, 5) * 10)),
  lambda_n = rbeta(1, p / 197.0 * 100, (1 - p / 197.0) * 100),
  r_b = runif(1, 0.05, 0.2),
  r_d = runif(1, 0.05, 0.2), 
  mid_b_uc = runif(1,0.2,0.7), mid_d_uc = runif(1,0.2,0.7),
  sigma_c = 1000 + abs(rnorm(1, 5000, 5000))
), simplify=FALSE) 

fit <- sampling(model, data = data, 
                control=list(adapt_delta=0.9, max_treedepth = 11),
                init = inits, 
                chains = chains, cores = chains, warmup = warmup, 
                iter = iter, refresh = 10, save_warmup = FALSE)

saveRDS(fit, file="fit_population_model.rds")
check_hmc_diagnostics(fit)
