library(rstan)
library(dplyr)

fit <- readRDS("fit_population_model.rds")

mu <- cbind(Year = 1647:1850,
            rbind(summary(fit,"mu_1647", use_cache = FALSE)[[1]][,c("2.5%","mean","97.5%", "se_mean", "n_eff", "Rhat")],
                  summary(fit,"mu", use_cache = FALSE)[[1]][,c("2.5%","mean","97.5%", "se_mean","n_eff", "Rhat")]))
write.csv(mu,file="mu.csv")

mu <- do.call("cbind", extract(fit, c("mu_1647", "mu")))
saveRDS(mu, file = "mu.rds")

lambda <- extract(fit, c("lambda_b", "lambda_d"))
saveRDS(lambda, file = "lambda.rds")

summary_mu <- apply(mu, 2, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))
summary_lambda_b <- apply(lambda$lambda_b, 2, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))
summary_lambda_d <- apply(lambda$lambda_d, 2, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))
saveRDS(summary_mu, file = "summary_mu.rds")
saveRDS(summary_lambda_d, file = "summary_lambda_d.rds")
saveRDS(summary_lambda_b, file = "summary_lambda_b.rds")

theta <- extract(fit, pars=c(
  "sigma_nu_b", "sigma_nu_d", "sigma_eta_b", "sigma_eta_d",
  "psi_bd", "psi_mu", "sigma_c",
  "mu_1647", "phi_d", "pi", 
  "phi_mu",
  "r_b", "r_d", "m_b", "m_d",
  "lc", "lambda_b[1]", "lambda_d[1]", "lambda_n"))
saveRDS(theta, file = "theta.rds")
  
  
sink("summary_stats.txt")
print(fit, pars = c(
  "sigma_nu_b", "sigma_nu_d", "sigma_eta_b", "sigma_eta_d",
  "psi_bd", "psi_mu", "sigma_c",
  "mu_1647", "phi_d", "pi", 
  "phi_mu",
  "r_b", "r_d", "m_b", "m_d",
  "lc", "lambda_b[1]", "lambda_d[1]", "lambda_n"), use_cache = FALSE)
sink()


# Data for posterior predictive check figure:

set.seed(1)
N <- 1000
n <- 203
source("benchmark_data.R")
enu_b <- exp(extract(fit, "nu_b")[[1]][1:N,,])
enu_d <- exp(extract(fit, "nu_d")[[1]][1:N,,])
psi_bd <- extract(fit, "psi_bd")[[1]][1:N]
l_b <- extract(fit, "lambda_b")[[1]][1:N,]
l_d <- extract(fit, "lambda_d")[[1]][1:N,]
fbd <- extract(fit, "phi_d")[[1]][1:N]
fr <- extract(fit, "pi")[[1]][1:N]
psi_mu <- extract(fit, "psi_mu")[[1]][1:N]
m0 <- extract(fit, "mu_1647")[[1]][1:N]
sigma_c <- extract(fit, "sigma_c")[[1]][1:N]
yrep <- matrix(NA, n+1, N)
brep <- matrix(NA, n, N)
drep <- matrix(NA, n, N)
mrep <- matrix(NA, n+1, N)
fcrep <- rep(NA, N)
for(i in 1:ncol(yrep)) {
  brep[, i] <- rgamma(n, psi_bd[i] * colSums(enu_b[i,,]), psi_bd[i])/l_b[i,]
  drep[-50, i] <- rgamma(n-1, psi_bd[i] * colSums(enu_d[i,,-50]), psi_bd[i])
  drep[50, i] <- rgamma(1, psi_bd[i] * sum(fbd[i]*enu_d[i,,50]), psi_bd[i])
  drep[, i] <- drep[, i] / l_d[i,]
  mrep[1, i] <- rgamma(1, psi_mu[i] * (m0[i] + brep[1,i] - drep[1,i] - soldiers[1]),
                       psi_mu[i])
  for(t in 2:50) { #1648-1696
    mrep[t,i] <- rgamma(1, psi_mu[i] * (mrep[t-1,i] + brep[t-1,i] -
                                          drep[t-1,i] - soldiers[t-1]), psi_mu[i]);
  }
  fcrep[i] <- (mrep[50,i] + brep[50,i] - soldiers[50] - fr[i] * mrep[49,i]) / drep[50,i]
  mrep[51,i] <- rgamma(1, psi_mu[i] * (mrep[50,i] + brep[50,i] - drep[50,i] *
                                         fcrep[i] - soldiers[50]), psi_mu[i]);
  for(t in 52:(n+1)) { #1648-1696
    mrep[t,i] <- rgamma(1, psi_mu[i] * (mrep[t-1,i] + brep[t-1,i] -
                                          drep[t-1,i] - soldiers[t-1]), psi_mu[i]);
  }
  yrep[, i] <- rnorm(n+1, mrep[, i], sigma_c[i])
}
saveRDS(yrep, file = "ppcheck_y.rds")

# Other results

idx <- sample(1:nrow(mu), size = 1000)

mu_sample <- ts(mu[idx,], start = 1647)
growth_rates <- colMeans(apply(window(mu_sample, end = 1690), 2, function(x) 100 * diff(x) / x[-length(x)]))

rates <- dplyr::case_when(
  growth_rates < 0 ~ "Negative growth",
  growth_rates > 0 & growth_rates < 0.5 ~ "Slow growth",
  growth_rates > 0.5 ~ "Fast growth")

df <- data.frame(population = c(mu_sample[1:44, ]), 
                 time = 1647:1690, replication = rep(1:ncol(mu_sample), each = 44),
                 growth = rep(factor(rates, levels = c("Fast growth", "Slow growth", "Negative growth"), 
                                     ordered = TRUE), each = 44))
growth_rates <- colMeans(apply(window(ts(t(mu), start = 1647), end = 1690), 2, 
                               function(x) 100 * diff(x) / x[-length(x)]))

sink("growth.txt")
prop.table(table(dplyr::case_when(
  growth_rates < 0 ~ "negative growth",
  growth_rates > 0 & growth_rates < 0.5 ~ "slow growth",
  growth_rates > 0.5 ~ "fast growth"))) * 100
sink()

###
sink("other_results.txt")
cbind(1647:1850,c(diff(apply(mu, 2, quantile, c(0.025, 0.975)))) / colMeans(mu) * 100)
difference <- apply(mu, 1, diff)
round(cbind(1648:1850, rowMeans(difference/t(mu[,-(n+1)]))*100, t(apply(difference/t(mu[,-(n+1)])*100,1, quantile, probs = c(0.025,0.975)))),1)

growth_rates <- colMeans(apply(mu[,1:44], 1, function(x) 100 * diff(x) / x[-length(x)]))

hist(growth_rates, breaks = 100, main = "Average annual growth rate 1647-1690")
round(c(mean = mean(growth_rates), quantile(growth_rates, c(0.025, 0.5, 0.975))), 3)

# when was the maximum in 1647-1699
x <- apply(mu[,1:52], 1, function(x) 1646 + c(which.min(x), which.max(x)))
rbind(mean = rowMeans(x), apply(x,1,quantile, prob= c(0.025,0.5, 0.975)))
# what was the value
(mean_mu <- mean(mu[,1674-1646]))

# find the peak population before famine, and when we got back there
x <- 1696 + apply(mu, 1, function(x) {
  which(x[51:204] > mean_mu)[1]
})
c(mean = mean(x), quantile(x, c(0.025,0.5, 0.975)))

x<-(mu[,51]/mu[,50]-1) * 100
c(mean = mean(x), quantile(x, c(0.025,0.5, 0.975)))

# find the proportion of highest and lowest point in 1690s
x <- apply(mu, 1, function(x) 1 - min(x[44:53])/max(x[44:53]))
c(mean = mean(x), quantile(x, c(0.025,0.5, 0.975)))


# find the proportion of highest and lowest point in 17th century
x <- apply(mu, 1, function(x) min(x[1:53])/max(x[1:53]))
c(mean = mean(x), quantile(x, c(0.025,0.5, 0.975)))


#1721
signif(c(mean = mean(mu[,75]), quantile(mu[,75], c(0.025,0.5, 0.975))),3)
#1700
signif(c(mean = mean(mu[,54]), quantile(mu[,54], c(0.025,0.5, 0.975))),3)
mean(mu[,54])-mean(mu[,75])

x <- apply(mu[,54:75], 1, function(x) (max(x) - min(x)))
round(c(mean = mean(x), quantile(x, c(0.025, 0.5, 0.975))))

x <- apply(mu[, 54:75], 1, function(x) 1699 + c(which.min(x), which.max(x)))
rbind(mean = rowMeans(x), apply(x,1,quantile, prob= c(0.025,0.5, 0.975)))


x <- apply(mu[, 54:75], 1, function(x) c(which.min(x)- which.max(x)))
c(mean = mean(x), quantile(x, prob= c(0.025,0.5, 0.975)))
mean(x < 0)*100
c(mean = mean(abs(x)), quantile(abs(x), prob= c(0.025,0.5, 0.975)))
#hist(x,breaks=100)
#ts.plot(ts(t(mu[which((x < 0))[1:10],54:75]),start=1700))
x <- apply(mu[,54:75], 1, function(x) (max(x) - min(x)))
round(c(mean = mean(x), quantile(x, c(0.025, 0.5, 0.975))))
#hist(mu[,75] - mu[,54], breaks = 100)
sink()
