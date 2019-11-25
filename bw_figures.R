library(rstan)
library(ggplot2)
library(ggthemes)
library(patchwork) # needs to be installed from github: thomasp85/patchwork
library(extrafont)
library(scales)
library(readxl)
library(dplyr)
library(svglite)

source("benchmark_data.R")
obs_births <- readRDS("births_by_parish.rda")
obs_deaths <- readRDS("deaths_by_parish.rda")


#font_import(paths = NULL, recursive = TRUE, prompt = FALSE,pattern = NULL)
loadfonts()
bold_text <- element_text(face = "bold", size = 12, family= "Palatino Linotype")
text <- element_text(size = 12, family= "Palatino Linotype")


n <- 203
parish_data <- data.frame(Year = 1648:1850, Count = 
                            c(rowSums(obs_births, na.rm=TRUE), rowSums(obs_deaths, na.rm=TRUE)),
                          Percentage = 100*c(rowSums(!is.na(obs_births)), rowSums(!is.na(obs_deaths)))/197,
                          Series = rep(c("Baptisms", "Burials"), each = n))


p1 <- ggplot(parish_data, aes(Year, Percentage, group = Series, color = Series)) +
  geom_line() + scale_y_log10("% of Parishes with Records") +
  theme_bw() +  scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.805, 0.13))

p2 <- ggplot(parish_data, aes(Year, Count, group = Series, color = Series)) +
  geom_line() + scale_y_log10() + 
  theme_bw() + scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.805, 0.13))

p <- p1 + p2
ggsave("Fig1.svg", p, width = 174, height = 174/2, unit = "mm")


## Fig2
tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800, 1850, by = 5))
pop_prior[!(time(pop_prior) %in% tabellverket)] <- NA

yrep <- readRDS("ppcheck_y.rds")

df0 <- data.frame(census = c(NA, pop_prior), time = 1647:1850)
df <- data.frame(population = c(yrep), time = 1647:1850, replication = rep(1:N, each = n + 1))

p <- ggplot(df, aes(x = time, y = population)) +
  geom_line(aes(colour = "Posterior Predictive Sample", group = replication), alpha = 0.03) +
  geom_point(data = df0, aes(x = time, y = census, colour = "Census")) +
  scale_y_continuous("Population",labels = comma) + scale_x_continuous("Year") + theme_bw() +
  theme_bw() + scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.position=c(0.9, 0.1))

svg("Fig2.svg", p,  width = 174, height = 174/2, unit = "mm")

lambda <- readRDS("lambda.rds")

lambda_data <- data.frame(
  Year = 1648:1850, 
  Lambda = rep(c("Births", "Deaths"), each = n),
  Value = colMeans(lambda), 
  lwr = apply(lambda, 2, quantile, 0.025),
  upr = apply(lambda, 2, quantile, 0.975))

p <- ggplot(lambda_data, aes(x = Year, y = Value, colour = Lambda)) + 
  geom_line() +
  geom_line(aes(y=lwr, colour = Lambda), linetype = "dashed") + 
  geom_line(aes(y=upr, colour = Lambda), linetype = "dashed") +
  scale_y_continuous("Value") + scale_x_continuous("Year") + 
  theme_bw() + scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.position=c(0.88, 0.16)) 

svg("Fig3.svg", p, width = 174, height = 174/2, unit = "mm")

mu <- readRDS("mu.rds")
musum <- t(apply(mu, 2, function(x) c(
  quantile(x, 0.025), 
  mean = mean(x), 
  quantile(x, 0.975))))


pop_data <- data.frame(Year = 1647:1850, 
                       Census = c(NA, pop_prior), 
                       Mean = colMeans(mu),
                       lwr = apply(mu, 2, quantile, 0.025),
                       upr = apply(mu, 2, quantile, 0.975))

p1 <- ggplot(pop_data[pop_data$Year <= 1750,], aes(x = Year)) + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + 
  theme_bw() + scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.position=c(0.85, 0.1)) 

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data, aes(x = Year)) + 
  geom_rect(xmin = xlim[1], xmax = xlim[2], 
            ymin = ylim[1], ymax = ylim[2], 
            color = "grey70", alpha = 0, linetype="dashed") + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous(NULL,labels=comma) + scale_x_continuous("Year") + 
  theme_bw() + scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.position=c(0.85, 0.1)) 


p <- p1 + p2 + plot_layout(ncol=2)
svg("Fig4.svg",width = 16, height = 8)
p
dev.off()

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
  scale_y_continuous("Yearly difference in population", label=comma) + scale_x_continuous("Year") + 
  theme_bw() + scale_colour_grey() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.position=c(0.88, 0.16)) 


svg("Fig5.svg",width = 16, height = 8)
p
dev.off()

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
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + 
  theme_bw() + scale_colour_grey() + 
  theme(legend.title=element_blank(), text = text, axis.title = bold_text, axis.text = text) +
  #scale_colour_manual(values = cols) + 
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
