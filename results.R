library(rstan)
library(ggplot2)
library(ggthemes)
library(patchwork) # needs to be installed from github: thomasp85/patchwork
library(extrafont)
library(scales)
library(readxl)
library(dplyr)

fit <- readRDS("fit_population_model.rds")
source("benchmark_data.R")
obs_births <- readRDS("births_by_parish.rda")
obs_deaths <- readRDS("deaths_by_parish.rda")


font_import(paths = NULL, recursive = TRUE, prompt = FALSE,pattern = NULL)
loadfonts()
bold_text <- element_text(face = "bold")
text <- element_text(size = 12)


mu <- cbind(Year = 1647:1850,
  rbind(summary(fit,"mu_1647", use_cache = FALSE)[[1]][,c("2.5%","mean","97.5%", "se_mean", "n_eff", "Rhat")],
    summary(fit,"mu", use_cache = FALSE)[[1]][,c("2.5%","mean","97.5%", "se_mean","n_eff", "Rhat")]))
write.csv(mu,file="mu.csv")


theta_monitor <- monitor(sims <- extract(fit, pars=c(
  "sigma_drift_b", "sigma_drift_d", "sigma_rw_b", "sigma_rw_d",
  "rate_bd", "rate_mu", "sigma_pop",
  "mu_1647", "famine", "famine_ratio", 
  "famine_correction",
  "r_b", "r_d", "m_b", "m_d",
  "lc", "lambda_b[1]", "lambda_d[1]", "lambda_n"), permuted = FALSE), warmup = 0)
saveRDS(theta_monitor, file = "theta_monitor.rds")
saveRDS(sims, file = "sims.rds")

mu <- extract(fit, c("mu_1647", "mu"), permuted = FALSE)
mu_monitor <- monitor(mu, warmup = 0)
saveRDS(mu_monitor, file = "mu_monitor.rds")

mu <- extract(fit, c("mu_1647", "mu"))
saveRDS(mu, file = "mu.rds")

n <- 203

parish_data <- data.frame(Year = 1648:1850, Count = 
    c(rowSums(obs_births, na.rm=TRUE), rowSums(obs_deaths, na.rm=TRUE)),
  Percentage = 100*c(rowSums(!is.na(obs_births)), rowSums(!is.na(obs_deaths)))/197,
  Series = rep(c("Baptisms", "Burials"), each = n))


p1 <- ggplot(parish_data, aes(Year, Percentage, group = Series, color = Series)) +
  geom_line() + scale_y_log10("% of parishes with records") +
  theme_bw() + scale_colour_few() + theme(text = text, axis.title=bold_text) +
  guides(color=FALSE)

p2 <- ggplot(parish_data, aes(Year, Count, group = Series, color = Series)) +
  geom_line() + scale_y_log10() + theme_bw() +
  scale_colour_few() + theme(text=text, axis.title=bold_text) +
  theme(legend.title=element_blank(),legend.position=c(0.85, 0.14))

svg("Fig1.svg", width = 16, height = 8)
p1 + p2
dev.off()

## Fig2
tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800, 1850, by = 5))
pop_prior[!(time(pop_prior) %in% tabellverket)] <- NA

N <- 1000
erw_b <- exp(extract(fit, "rw_b")[[1]][1:N,,])
erw_d <- exp(extract(fit, "rw_d")[[1]][1:N,,])
rate_bd <- extract(fit, "rate_bd")[[1]][1:N]
l_b <- extract(fit, "lambda_b")[[1]][1:N,]
l_d <- extract(fit, "lambda_d")[[1]][1:N,]
fbd <- exp(extract(fit, "famine")[[1]][1:N])
fr <- extract(fit, "famine_ratio")[[1]][1:N]
rate_mu <- extract(fit, "rate_mu")[[1]][1:N]
m0 <- extract(fit, "mu_1647")[[1]][1:N]
sigma_pop <- extract(fit, "sigma_pop")[[1]][1:N]
n <- 203
yrep <- matrix(NA, n+1, N)
brep <- matrix(NA, n, N)
drep <- matrix(NA, n, N)
mrep <- matrix(NA, n+1, N)
fcrep <- rep(NA, N)
set.seed(1)
for(i in 1:ncol(yrep)) {
  brep[, i] <- rgamma(n, rate_bd[i] * colSums(erw_b[i,,]), rate_bd[i])/l_b[i,]
  drep[-50, i] <- rgamma(n-1, rate_bd[i] * colSums(erw_d[i,,-50]), rate_bd[i])
  drep[50, i] <- rgamma(1, rate_bd[i] * sum(fbd[i]*erw_d[i,,50]), rate_bd[i])
  drep[, i] <- drep[, i] / l_d[i,]
  mrep[1, i] <- rgamma(1, rate_mu[i] * (m0[i] + brep[1,i] - drep[1,i] - soldiers[1]),
                       rate_mu[i])
  for(t in 2:50) { #1648-1696
    mrep[t,i] <- rgamma(1, rate_mu[i] * (mrep[t-1,i] + brep[t-1,i] -
                                           drep[t-1,i] - soldiers[t-1]), rate_mu[i]);
  }
  fcrep[i] <- (mrep[50,i] + brep[50,i] - soldiers[50] - fr[i] * mrep[49,i]) / drep[50,i]
  mrep[51,i] <- rgamma(1, rate_mu[i] * (mrep[50,i] + brep[50,i] - drep[50,i] *
                                          fcrep[i] - soldiers[50]), rate_mu[i]);
  for(t in 52:(n+1)) { #1648-1696
    mrep[t,i] <- rgamma(1, rate_mu[i] * (mrep[t-1,i] + brep[t-1,i] -
                                           drep[t-1,i] - soldiers[t-1]), rate_mu[i]);
  }
  yrep[, i] <- rnorm(n+1, mrep[, i], sigma_pop[i])
}

save(erw_b,
     erw_d,
     rate_bd,
     l_b,
     l_d,
     fbd,
     fr,
     rate_mu,
     m0,
     sigma_pop,
     yrep,
     brep,
     drep,
     mrep,
     fcrep, file = "ppcheck.rds") # Note this file is too big for Github

df0 <- data.frame(census = c(NA,pop_prior), time = 1647:1850)
df <- data.frame(population = c(yrep), time = 1647:1850, replication = rep(1:N, each = n+1))
p <- ggplot(df, aes(x = time, y = population)) +
  geom_line(aes(colour = "Posterior predictive sample",group = replication), alpha=0.03) +
  geom_point(data=df0, aes(x = time, y = census, colour = "Census")) +
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) +
  theme(legend.position=c(0.9, 0.1))

svg("Fig2.svg", width = 16, height = 8)
p
dev.off()

lambda_data <- data.frame(
  Year = 1648:1850, 
  Lambda = rep(c("Births", "Deaths"), each = n),
  Value = summary(fit,pars=c("lambda_b", "lambda_d"), use_cache = FALSE)[[1]][,"mean"], 
  lwr = summary(fit,pars=c("lambda_b", "lambda_d"), use_cache = FALSE)[[1]][,"2.5%"],
  upr = summary(fit,pars=c("lambda_b", "lambda_d"), use_cache = FALSE)[[1]][,"97.5%"])

p <- ggplot(lambda_data, aes(x = Year, y = Value, colour = Lambda)) + 
  geom_line() +
  geom_line(aes(y=lwr, colour = Lambda), linetype = "dashed") + 
  geom_line(aes(y=upr, colour = Lambda), linetype = "dashed") +
  scale_y_continuous("Value") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 

svg("Fig3.svg",width = 16, height = 8)
p
dev.off()

lambda_monitor <- monitor(extract(fit, c("lambda_b", "lambda_d"), permuted = FALSE), warmup = 0)
saveRDS(lambda_monitor, file = "lambda_monitor.rds")

musum <- rbind(
  summary(fit,"mu_1647", use_cache = FALSE)[[1]][,c("2.5%","mean","97.5%")],
  summary(fit,"mu", use_cache = FALSE)[[1]][,c("2.5%","mean","97.5%")])

pop_data <- data.frame(Year = 1647:1850, 
  Census = c(NA, pop_prior), 
  Mean = musum[,"mean"],
  lwr = musum[,"2.5%"],
  upr = musum[,"97.5%"])
p1 <- ggplot(pop_data[pop_data$Year <= 1750,], aes(x = Year)) + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.85, 0.1)) 

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data, aes(x = Year)) + 
  geom_rect(xmin=xlim[1], xmax=xlim[2], ymin=ylim[1], ymax=ylim[2], 
    color = "grey70", alpha=0, linetype="dashed") + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous(NULL,labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.85, 0.1)) 


p <- p1 + p2 + plot_layout(ncol=2)
svg("Fig4.svg",width = 16, height = 8)
p
dev.off()

###

mu <- cbind(extract(fit, "mu_1647")$mu, extract(fit, "mu")$mu)

difference <- apply(mu, 1, diff)

diff_data <- data.frame(Year = 1648:1850, 
  Prior = rowSums(obs_births,na.rm=TRUE) - rowSums(obs_deaths,na.rm=TRUE) - soldiers, 
  Mean = apply(difference, 1, mean),
  lwr = apply(difference, 1, quantile, 0.025),
  upr = apply(difference, 1, quantile, 0.975))

p <- ggplot(diff_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Parish records")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Yearly difference in population", label=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 


svg("Fig5.svg",width = 16, height = 8)
p
dev.off()


###

prior <- readxl::read_xlsx("Aiemmat estimaatit.xlsx")[[2]]
pop_data <- data.frame(Year = 1647:1850, 
                       Census = c(NA, pop_prior), 
                       Literature = c(rep(NA, 43), prior),
                       Mean = musum[,"mean"],
                       lwr = musum[,"2.5%"],
                       upr = musum[,"97.5%"])
cols <-few_pal("Medium")(3)
names(cols) <- c("Census", "This study", "Earlier literature")
p1 <- ggplot(pop_data, aes(x = Year)) + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  geom_line(aes(y = Literature, colour = "Earlier literature")) +
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_manual(values = cols) + 
  theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.85, 0.1)) 

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range


pop_data <- data.frame(Year = 1648:1850, 
  Prior = c(rep(NA, 43), diff(prior)/prior[-length(prior)]*100), 
  Mean = rowMeans(difference/t(mu[,-(n+1)]))*100,
  lwr = apply(difference/t(mu[,-(n+1)])*100,1,quantile, 0.025),
  upr = apply(difference/t(mu[,-(n+1)])*100,1,quantile, 0.975))

p2 <- ggplot(pop_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Earlier literature")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Growth rate") + scale_x_continuous("Year") + theme_bw() +
  coord_cartesian(ylim = c(-22, 4)) +
  scale_colour_manual(values = cols) + 
  theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 

svg("Fig6.svg", width = 16, height = 8)
p1 + p2
dev.off()


mu <- readRDS("mu.rds")
idx <- sample(1:nrow(mu[[1]]), size = 1000)

mu_sample <- ts(t(cbind(mu$mu_1647[idx], mu$mu[idx,])), start = 1647)
growth_rates <- colMeans(apply(window(mu_sample, end = 1690), 2, function(x) 100 * diff(x) / x[-length(x)]))

rates <- dplyr::case_when(
  growth_rates < 0 ~ "Negative growth",
  growth_rates > 0 & growth_rates < 0.5 ~ "Slow growth",
  growth_rates > 0.5 ~ "Fast growth")

df <- data.frame(population = c(mu_sample[1:44, ]), 
                 time = 1647:1690, replication = rep(1:ncol(mu_sample), each = 44),
                 growth = 
                   rep(factor(rates, levels = c("Fast growth", "Slow growth", "Negative growth"), 
                              ordered = TRUE), each = 44))

mu <- ts(t(cbind(mu$mu_1647, mu$mu)), start = 1647)
growth_rates <- colMeans(apply(window(mu, end = 1690), 2, function(x) 100 * diff(x) / x[-length(x)]))
negative_growth <- which(growth_rates < 0)
slow_growth <- which(growth_rates > 0 & growth_rates < 0.5)
fast_growth <- which(growth_rates > 0.5)
sink("growth.txt")
prop.table(table(dplyr::case_when(
  growth_rates < 0 ~ "negative growth",
  growth_rates > 0 & growth_rates < 0.5 ~ "slow growth",
  growth_rates > 0.5 ~ "fast growth"))) * 100
sink()
means <- data.frame(population = c(rowMeans(mu[1:44,negative_growth]), 
                                   rowMeans(mu[1:44,slow_growth]),
                                   rowMeans(mu[1:44,fast_growth])),
                    time = 1647:1690, replication = rep(1:3, each = 44),
                    growth = rep(factor(c("Negative growth", "Slow growth", "Fast growth")), each = 44))

cols <-few_pal("Medium")(3)
names(cols) <- c("Fast growth", "Slow growth", "Negative growth")

p <-ggplot(df, aes(x = time, y = population, colour = growth, group = replication)) + 
  geom_line(alpha=0.25) + geom_line(data = means, size = 2) + 
  scale_y_continuous("Population",labels = comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_manual(values = cols, name = "Growth rate") + theme(text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.1, 0.85)) 

svg("Fig7.svg", width = 16, height = 8)
p
dev.off()


