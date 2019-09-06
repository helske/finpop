source("../benchmark_data.R")
obs_births <- readRDS("births_by_parish.rda")
obs_deaths <- readRDS("deaths_by_parish.rda")
n <- 203
install.packages("devtools")
devtools::install_github("thomasp85/patchwork")
install.packages(c("ggplot2", "ggthemes", "extrafont"))
library(ggplot2)
library(ggthemes)
library(patchwork)
library(extrafont)
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

svg("../revision2019/Fig1.svg", width = 8.4,height=4)
p1 + p2
dev.off()

#fit <- readRDS("paper/results.rda")

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

svg("../revision2019/figures2903/Fig2.svg", width = 6,height=4)
p
dev.off()

tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800,1830, by = 5), 1831:1850)
pop_prior[!(time(pop_prior) %in% tabellverket)] <- NA
mu <- rbind(
  summary(fit,"mu_1647")[[1]][,c("2.5%","mean","97.5%")],
  summary(fit,"mu")[[1]][,c("2.5%","mean","97.5%")])

prior <- readxl::read_xlsx("../revision2019/Aiemmat estimaatit.xlsx")[[2]]
pop_data <- data.frame(Year = 1647:1850, 
  Census = c(NA, pop_prior), 
  Literature = c(rep(NA, 43), prior),
  Mean = mu[,"mean"],
  lwr = mu[,"2.5%"],
  upr = mu[,"97.5%"])
library(scales)
p1 <- ggplot(pop_data[pop_data$Year <= 1750,], aes(x = Year)) + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  geom_line(aes(y = Literature, colour = "Earlier literature")) +
  scale_y_continuous("Population",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.9, 0.1)) 

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data, aes(x = Year)) + 
  geom_rect(xmin=xlim[1], xmax=xlim[2], ymin=ylim[1], ymax=ylim[2], 
    color = "grey70", alpha=0, linetype="dashed") + 
  geom_point(aes(y = Census, colour = "Census")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  geom_line(aes(y = Literature, colour = "Earlier literature")) +
  scale_y_continuous(NULL,labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.9, 0.1)) 


p <- p1 + p2 + plot_layout(ncol=2)
ggsave(p, file = "../revision2019/figures2903/Fig3.pdf", width = 10,height=5)

###

births <- summary(fit,"births")[[1]][,c("2.5%","mean","97.5%")]

birth_data <- data.frame(Year = 1648:1850, 
  Prior = rowSums(obs_births,na.rm=TRUE), 
  Mean = births[,"mean"],
  lwr = births[,"2.5%"],
  upr = births[,"97.5%"])


p1 <- ggplot(birth_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Parish records")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Births",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  guides(color=FALSE)

deaths <- summary(fit,"deaths")[[1]][,c("2.5%","mean","97.5%")]

death_data <- data.frame(Year = 1648:1850, 
  Prior = rowSums(obs_deaths,na.rm=TRUE), 
  Mean = deaths[,"mean"],
  lwr = deaths[,"2.5%"],
  upr = deaths[,"97.5%"])


p2 <- ggplot(death_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Parish records")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Deaths",labels=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  guides(color=FALSE)


b <- extract(fit, "births", permuted = FALSE)
d <- extract(fit, "deaths", permuted = FALSE)

difference <- b - d

diff_data <- data.frame(Year = 1648:1850, 
  Prior = rowSums(obs_births,na.rm=TRUE) - rowSums(obs_deaths,na.rm=TRUE), 
  Mean = apply(difference, 3, mean),
  lwr = apply(difference, 3, quantile, 0.025),
  upr = apply(difference, 3, quantile, 0.975))

p3 <- ggplot(diff_data, aes(x = Year)) + 
  geom_line(aes(y = Prior, colour = "Parish records")) +
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Yearly difference in births and deaths", label=comma) + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.16)) 


svg("../revision2019/figures2903/Fig3parish.svg", width = 8.4,height=4)
p <- (p1 | p2 ) / p3
dev.off()




miss_b <- summary(fit,"miss_b")[[1]][,c("2.5%","mean","97.5%")]
miss_d <- summary(fit,"miss_d")[[1]][,c("2.5%","mean","97.5%")]
ts.plot(miss_b)
missb_data <- data.frame(Year = 1648:1850, 
  Mean = miss_b[,"mean"],
  lwr = miss_b[,"2.5%"],
  upr = miss_b[,"97.5%"])


p1 <- ggplot(missb_data, aes(x = Year)) + 
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Births") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  guides(color=FALSE)

missd_data <- data.frame(Year = 1648:1850, 
  Mean = miss_d[,"mean"],
  lwr = miss_d[,"2.5%"],
  upr = miss_d[,"97.5%"])


p2 <- ggplot(missd_data, aes(x = Year)) + 
  geom_line(aes(y = Mean, colour = "This study")) +
  geom_line(aes(y = lwr, colour = "This study"), linetype = "dashed") + 
  geom_line(aes(y = upr, colour = "This study"), linetype = "dashed") +
  scale_y_continuous("Deaths") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  guides(color=FALSE)
p1 + p2

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
svg("paper/Fig4.svg", width = 6,height=4)
p
dev.off()

## Fig5
mu_rep <- extract(fit, pars="mu_rep")[[1]]
pp_data <- data.frame(Year = 1648:1850,
  Population = c(t(mu_rep)), 
  replication = rep(1:(length(mu_rep)/203), each = 203))
library(dplyr)
pp_data %>% group_by(Year) %>% 
  summarise(mean=mean(Population), 
    q1l=quantile(Population, probs=0.025),
    q2l=quantile(Population, probs=0.125),
    q3l=quantile(Population, probs=0.25),
    q3u=quantile(Population, probs=0.75),
    q2u=quantile(Population, probs=0.875),
    q1u=quantile(Population, probs=0.975)) -> new_data

p <- ggplot(new_data, aes(x = Year, y = mean)) + 
  geom_ribbon(aes(ymin=q1l, ymax=q1u),
    fill = few_pal(palette = "medium")(2)[2],
    alpha=0.2,       #transparency
    linetype=0,      #solid, dashed or other line types
    colour="grey", #border line color
    size=0.1) +    #fill color
  geom_ribbon(aes(ymin=q2l, ymax=q2u),
    fill = few_pal(palette = "medium")(2)[2],
    alpha=0.2,       #transparency
    linetype=0,      #solid, dashed or other line types
    colour="grey", #border line color
    size=0.1) +    #fill color
  geom_ribbon(aes(ymin=q3l, ymax=q3u),
    fill = few_pal(palette = "medium")(2)[2],
    alpha=0.2,       #transparency
    linetype=0,      #solid, dashed or other line types
    colour="grey", #border line color
    size=0.1) +    #fill color
  geom_line(aes(colour = "Replication")) +
  geom_line(aes(y = pop_prior, colour = "Reference")) +
  scale_y_continuous("Population") + scale_x_continuous("Year")  + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.88, 0.1)) + guides(fill=FALSE)
p

svg("paper/Fig5.svg", width = 6,height=4)
p
dev.off()


m1700 <- readRDS("m1700.rda")
m1750 <- readRDS("m1750.rda")
f <- factor(c("midpoint 1700", "midpoint 1750", "midpoint 1800", "Reference"))
levels(f)<-rev(levels(f))
pop_data <- data.frame(Year = 1648:1850, 
  series = rep(rev(f),each=203),
  Mean = c(m1700[,1], m1750[,1], summary(fit,pars="mu")[[1]][,"mean"], pop_prior),
  lwr = c(m1700[,2], m1750[,2], summary(fit,pars="mu")[[1]][,"2.5%"], rep(NA,203)),
  upr = c(m1700[,3], m1750[,3], summary(fit,pars="mu")[[1]][,"97.5%"], rep(NA,203)))


p1 <- ggplot(pop_data[pop_data$Year <= 1750,], aes(x = Year, y = Mean, colour = series)) + 
  geom_line() +
  geom_line(aes(y = lwr), linetype = "dashed") + 
  geom_line(aes(y = upr), linetype = "dashed") +
  scale_y_continuous("Population") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  guides(color=FALSE)

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data, aes(x = Year, y = Mean, colour = series)) + 
  geom_rect(xmin=xlim[1], xmax=xlim[2], ymin=ylim[1], ymax=ylim[2], 
    color = "grey", alpha=0, linetype="dotted") + 
  geom_line() +
  geom_line(aes(y = lwr), linetype = "dashed") + 
  geom_line(aes(y = upr), linetype = "dashed") +
  scale_y_continuous("Population") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.2, 0.8)) 

svg("paper/Fig6.svg", width = 8.4,height=4)
p1 + p2
dev.off()


bold_text <- element_text(face = "bold")
text <- element_text(size = 12)

m <- summary(readRDS("paper/results.rda"),pars="mu")[[1]]
m10000 <- readRDS("m10000.rda")
m50000 <- readRDS("m50000.rda")
f <- factor(c("Reference", "sigma 30000", "sigma 10000", "sigma 50000"), levels = c("Reference", "sigma 30000", "sigma 10000", "sigma 50000"))

pop_data <- data.frame(Year = 1648:1850, 
  series = rep(f,each=203),
  Mean = c(pop_prior, m[,"mean"],  m10000[,1], m50000[,1]),
  lwr = c(rep(NA,203), m[,"2.5%"], m10000[,2], m50000[,2]),
  upr = c(rep(NA,203), m[,"97.5%"], m10000[,3], m50000[,3]))


p1 <- ggplot(pop_data[pop_data$Year <= 1750,], aes(x = Year, y = Mean, colour = series)) + 
  geom_line() +
  geom_line(aes(y = lwr), linetype = "dashed") + 
  geom_line(aes(y = upr), linetype = "dashed") +
  scale_y_continuous("Population") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  guides(color=FALSE)

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data, aes(x = Year, y = Mean, colour = series)) + 
  geom_rect(xmin=xlim[1], xmax=xlim[2], ymin=ylim[1], ymax=ylim[2], 
    color = "grey", alpha=0, linetype="dotted") + 
  geom_line() +
  geom_line(aes(y = lwr), linetype = "dashed") + 
  geom_line(aes(y = upr), linetype = "dashed") +
  scale_y_continuous("Population") + scale_x_continuous("Year") + theme_bw() +
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text) + 
  theme(legend.position=c(0.2, 0.8)) 



svg("paper/FigA2.1.svg", width = 8.4,height=4)
p1 + p2
dev.off()

bold_text <- element_text(face = "bold")
text <- element_text(size = 12)
dat <- data.frame(Year = 1648:1850, 
  SD = summary(fit,pars="mu")[[1]][,"sd"])
p <- ggplot(dat, aes(Year, SD)) + geom_line() + theme_bw() +  scale_y_continuous(expression(bold("Posterior standard deviation of "~mu[t]))) + 
  scale_colour_few() + theme(legend.title=element_blank(), text=text, axis.title=bold_text)
p

svg("paper/FigA2.2.svg", width = 8.4,height=4)
p
dev.off()
