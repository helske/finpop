
source("benchmark_data.R")
obs_births <- readRDS("births_by_parish.rda")
obs_deaths <- readRDS("deaths_by_parish.rda")

library(rstan)
fit <- readRDS("fit_final.rds")
get_stanmodel(fit)
check_hmc_diagnostics(fit)

mu <- cbind(Year = 1647:1850,
  rbind(summary(fit,"mu_1647")[[1]][,c("2.5%","mean","97.5%", "se_mean", "n_eff", "Rhat")],
    summary(fit,"mu")[[1]][,c("2.5%","mean","97.5%", "se_mean","n_eff", "Rhat")]))
write.csv(mu,file="mu.csv")

#sigma_rw, rate_bd, rate_mu. _sigma_pop_!!
print(fit, pars=c(
  "sigma_drift_b", "sigma_drift_d", "sigma_rw_b", "sigma_rw_d",
  "rate_bd", "rate_mu", "sigma_pop",
  "mu_1647", "famine", "famine_ratio", 
  "famine_correction",
  "r_b", "r_d", "m_b", "m_d",
  "lc", "lambda_b[1]", "lambda_d[1]", "lambda_n"))

monitor(extract(fit, pars=c(
  "sigma_drift_b", "sigma_drift_d", "sigma_rw_b", "sigma_rw_d",
  "rate_bd", "rate_mu", "sigma_pop",
  "mu_1647", "famine", "famine_ratio", 
  "famine_correction",
  "r_b", "r_d", "m_b", "m_d",
  "lc", "lambda_b[1]", "lambda_d[1]", "lambda_n"), permuted = FALSE), warmup = 0)


tail(sort(summary(fit)[[1]][, c("Rhat")]),15)
head(sort(summary(fit)[[1]][, c("n_eff")]),15)

n <- 203

library(ggplot2)
library(ggthemes)
library(patchwork)
library(extrafont)
library(scales)

#font_import(paths = NULL, recursive = TRUE, prompt = TRUE,pattern = NULL)
loadfonts()
bold_text <- element_text(face = "bold")
text <- element_text(size = 12)

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

svg("Fig1.svg", width = 8.4,height=4)
p1 + p2
dev.off()


## Fig2


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

summary(fcrep)
mean(colSums(is.na(yrep))>0)*100
yrep[is.na(yrep)] <- 0

tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800, 1850, by = 5))
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



db <- exp(extract(fit,pars=c("drift_b"))[[1]])
dd <- exp(extract(fit,pars=c("drift_d"))[[1]])

drift_data <- data.frame(
  Year = 1648:1850, 
  Series = rep(c("Births", "Deaths"), each = n),
  Value = c(colMeans(db), colMeans(dd)), 
  lwr = c(apply(db,2,quantile,0.025), apply(dd,2,quantile,0.025)),
  upr = c(apply(db,2,quantile,0.975), apply(dd,2,quantile,0.975)))

p <- ggplot(drift_data, aes(x = Year, y = Value, colour = Series)) + 
  geom_line() +
  geom_line(aes(y=lwr, colour = Series), linetype = "dashed") + 
  geom_line(aes(y=upr, colour = Series), linetype = "dashed") +
  scale_y_continuous("Growth rate") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) + facet_wrap(~Series)

library(rstan)

lambda_data <- data.frame(
  Year = 1648:1850, 
  Lambda = rep(c("Births", "Deaths"), each = n),
  Value = summary(fit,pars=c("lambda_b", "lambda_d"))[[1]][,"mean"], 
  lwr = summary(fit,pars=c("lambda_b", "lambda_d"))[[1]][,"2.5%"],
  upr = summary(fit,pars=c("lambda_b", "lambda_d"))[[1]][,"97.5%"])

p <- ggplot(lambda_data, aes(x = Year, y = Value, colour = Lambda)) + 
  geom_line() +
  geom_line(aes(y=lwr, colour = Lambda), linetype = "dashed") + 
  geom_line(aes(y=upr, colour = Lambda), linetype = "dashed") +
  scale_y_continuous("Value") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 

svg("../final/Fig3.svg", width = 6,height=4)
p
dev.off()

pop_prior[!(time(pop_prior) %in% tabellverket)] <- NA
mu <- rbind(
  summary(fit,"mu_1647")[[1]][,c("2.5%","mean","97.5%")],
  summary(fit,"mu")[[1]][,c("2.5%","mean","97.5%")])

prior <- readxl::read_xlsx("Aiemmat estimaatit.xlsx")[[2]]
pop_data <- data.frame(Year = 1647:1850, 
  Census = c(NA, pop_prior), 
  #Literature = c(rep(NA, 43), prior),
  Mean = mu[,"mean"],
  lwr = mu[,"2.5%"],
  upr = mu[,"97.5%"])
library(scales)
p1 <- ggplot(pop_data[pop_data$Year <= 1750,], aes(x = Year)) + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
 # geom_line(aes(y = Literature, colour = "Earlier literature")) +
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
 # geom_line(aes(y = Literature, colour = "Earlier literature")) +
  scale_y_continuous(NULL,labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.85, 0.1)) 


p <- p1 + p2 + plot_layout(ncol=2)
ggsave(p, file = "Fig4.svg", width = 10,height=5)

###

b <- extract(fit, "births", permuted = FALSE)
d <- extract(fit, "deaths", permuted = FALSE)

difference <- b - d

diff_data <- data.frame(Year = 1648:1850, 
  Prior = rowSums(obs_births,na.rm=TRUE) - rowSums(obs_deaths,na.rm=TRUE), 
  Mean = apply(difference, 3, mean),
  lwr = apply(difference, 3, quantile, 0.025),
  upr = apply(difference, 3, quantile, 0.975))

p <- ggplot(diff_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Parish records")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Yearly difference in births and deaths", label=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 


svg("Fig5.svg", width = 8.4,height=4)
p
dev.off()


###
m <- summary(fit,pars="mu")[[1]][,"mean"]
pop_data <- data.frame(Year = 1649:1850, 
  Prior = diff(pop_prior)/pop_prior[-n]*100, 
  Mean = diff(m)/m[-n]*100)

p <- ggplot(pop_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Reference")) +
  geom_line(aes(y = Mean, colour = "This study"))+
  scale_y_continuous("Growth rate") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 
svg("paper/Fig6.svg", width = 6,height=4)
p
dev.off()

