setwd("C:/repos/finpop/new_division/")
source("../benchmark_data.R")
births <- readRDS("births_by_parish.rda")
deaths <- readRDS("deaths_by_parish.rda")

library(rstan)
fit <- readRDS("../revision2019/fit2903.rds")
N <- 1000
erw_b <- exp(extract(fit, "rw_b")[[1]][1:N,,])
erw_d <- exp(extract(fit, "rw_d")[[1]][1:N,,])
rate_b <- extract(fit, "rate_b")[[1]][1:N]
rate_d <- extract(fit, "rate_d")[[1]][1:N]
l_b <- extract(fit, "gamma_b")[[1]][1:N,]
l_d <- extract(fit, "gamma_d")[[1]][1:N,]
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
for(i in 1:ncol(yrep)) {
  brep[, i] <- rgamma(n, rate_b[i] * colSums(erw_b[i,,]), rate_b[i])/l_b[i,]
  drep[-50, i] <- rgamma(n-1, rate_d[i] * colSums(erw_d[i,,-50]), rate_d[i])
  drep[50, i] <- rgamma(1, rate_d[i] * sum(fbd[i]*erw_d[i,,50]), rate_d[i])
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
  yrep[, i] <- rnorm(n+1, mrep[, i], 1e3*sigma_pop[i])
}

summary(fcrep)
mean(colSums(is.na(yrep))>0)*100
yrep[is.na(yrep)] <- 0

tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800,1830, by = 5), 1831:1850)
pop_prior[!(time(pop_prior) %in% tabellverket)] <- NA

library(ggplot2)
library(ggthemes)
library(scales)
library(patchwork)
library(extrafont)
#font_import(paths = NULL, recursive = TRUE, prompt = TRUE,pattern = NULL)
loadfonts()
bold_text <- element_text(face = "bold")
text <- element_text(size = 12)
df0 <- data.frame(census = c(NA,pop_prior), time = 1647:1850)
df <- data.frame(population = c(yrep), time = 1647:1850, replication = rep(1:N, each = n+1))
ggplot(df, aes(x = time, y = population)) + 
  geom_line(aes(colour = "Posterior predictive sample",group = replication), alpha=0.03) +
  geom_point(data=df0, aes(x = time, y = census, colour = "Census")) +
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.9, 0.1)) 

####
  