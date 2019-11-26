library(rstan)
library(ggplot2)
library(ggthemes)
library(patchwork) # needs to be installed from github: thomasp85/patchwork
library(extrafont)
library(scales)
library(readxl)
library(dplyr)

width <- 174

source("benchmark_data.R")
obs_births <- readRDS("births_by_parish.rda")
obs_deaths <- readRDS("deaths_by_parish.rda")


#font_import(paths = NULL, recursive = TRUE, prompt = FALSE,pattern = NULL)
loadfonts(device = "pdf")
bold_text <- element_text(face = "bold", size = 12, family= "Arial")
text <- element_text(size = 12, family= "Arial")


n <- 203
parish_data <- data.frame(Year = 1648:1850, Count = 
                            c(rowSums(obs_births, na.rm=TRUE), rowSums(obs_deaths, na.rm=TRUE)),
                          Percentage = 100*c(rowSums(!is.na(obs_births)), rowSums(!is.na(obs_deaths)))/197,
                          Series = rep(c("Baptisms", "Burials"), each = n))


p1 <- ggplot(parish_data, aes(Year, Percentage, group = Series, color = Series)) +
  geom_line() + scale_y_log10("% of Parishes with Records") +
  theme_bw() + scale_colour_grey(start = 0, end = 0.7) + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.795, 0.13))

p2 <- ggplot(parish_data, aes(Year, Count, group = Series, color = Series)) +
  geom_line() + scale_y_log10() + 
  theme_bw() + scale_colour_grey(start = 0, end = 0.7) + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.795, 0.13))

p <- p1 + p2
ggsave("Fig1.pdf", p, width = width, height = width / 2, unit = "mm", device = cairo_pdf)


## Fig2
tabellverket <- c(1749, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1780, seq(1800, 1850, by = 5))
pop_prior[!(time(pop_prior) %in% tabellverket)] <- NA

yrep <- readRDS("ppcheck_y.rds")


df <- data.frame(population = c(c(NA, pop_prior), c(yrep)) / 1000, 
                 time = 1647:1850,
                 type = rep(c("Census", "Posterior Predictive Sample"), times = c(n+1, (n+1) * ncol(yrep))),
                 replication = rep(factor(0:ncol(yrep)), each = n + 1))


p <- ggplot(df, aes(x = time, y = population, group = replication)) + 
  geom_line(aes(linetype = type, colour = type), alpha = 0.01) +
  geom_point(aes(shape = type, colour = type, fill = type)) +
  scale_linetype_manual(breaks = c("Census", "Posterior Predictive Sample"), 
                        values = c("blank", "solid")) + 
  scale_shape_manual(breaks = c("Census", "Posterior Predictive Sample"), 
                     values = c(21, NA)) + 
  scale_fill_manual(breaks = c("Census", "Posterior Predictive Sample"), 
                    values = c("white", "black"), aesthetics = "fill") + 
  scale_colour_manual(breaks = c("Census", "Posterior Predictive Sample"), 
                      values = c("black", "black"), aesthetics = "colour") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_y_continuous("Population (in thousands)", labels = comma) +
  scale_x_continuous("Year", limits = c(1647, 1850), expand = c(0.01, 0.01)) +
  theme_bw() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.817, 0.13)) +
  theme(plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt"))


ggsave("Fig2.pdf", p,  width = width, height = width / 2, unit = "mm", device = cairo_pdf)

lambda_b <- readRDS("summary_lambda_b.rds")
lambda_d <- readRDS("summary_lambda_d.rds")

lambda_data <- data.frame(
  Year = 1648:1850, 
  Lambda = rep(c("Births", "Deaths"), each = n),
  Value = c(lambda_b["mean",], lambda_d["mean",]), 
  lwr = c(lambda_b["2.5%",], lambda_d["2.5%",]),
  upr = c(lambda_b["97.5%",], lambda_d["97.5%",]))

p <- ggplot(lambda_data, aes(x = Year, y = Value, colour = Lambda)) + 
  geom_line() +
  geom_line(aes(y=lwr, colour = Lambda), linetype = "dashed") + 
  geom_line(aes(y=upr, colour = Lambda), linetype = "dashed") +
  scale_y_continuous(expression(bold("Value of " ~ lambda))) + 
  scale_x_continuous("Year", expand = c(0.01, 0.01)) + 
  theme_bw() + scale_colour_grey(start = 0, end = 0.7) + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.92, 0.13)) +
  theme(plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt"))

ggsave("Fig3.pdf", p, width = width, height = width / 2, unit = "mm", device = cairo_pdf)

mu <- readRDS("summary_mu.rds")

pop_data <- data.frame(year = 1647:1850, 
                       value = c(NA, pop_prior,  mu["mean",], mu["2.5%",], mu["97.5%",]) / 1000,
                       series = rep(factor(c("Census", "Posterior Mean", 
                                             "95% Posterior Interval", "95% Posterior Interval"),
                                           levels = c("Census", "Posterior Mean", 
                                                      "95% Posterior Interval"), 
                                           ordered = TRUE), each = n + 1),
                       group = rep(factor(1:4), each = n + 1))

p1 <- ggplot(pop_data[pop_data$year <= 1750,], aes(x = year, y = value, colour = series)) + 
  geom_line(aes(linetype = series, group = group)) + 
  geom_point(aes(shape = series, fill = series)) +
  scale_y_continuous("Population (in thousands)", labels = comma) + 
  scale_x_continuous("Year") + 
  theme_bw() +  
  scale_linetype_manual(values = c("Census" = "blank", "Posterior Mean" = "solid", 
                                   "95% Posterior Interval" = "dashed")) +
  scale_shape_manual(values = c("Census" = 21, "Posterior Mean" = NA, 
                                "95% Posterior Interval" = NA)) +
  scale_colour_manual(values = c("Census" = "white", "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "fill") + 
  scale_colour_manual(values = c("Census" = "black", "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "colour") + 
  # scale_alpha_manual(values = c("Census" = 1, "Posterior Mean" = 1, 
  #                               "95% Posterior Interval" = 0.5)) + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.position = "none")

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data, aes(x = year, y = value, colour = series)) + 
  geom_rect(xmin = xlim[1], xmax = xlim[2], 
            ymin = ylim[1], ymax = ylim[2], 
            color = "grey90", alpha = 0, linetype="dotted") + 
  geom_line(aes(linetype = series, group = group)) + 
  geom_point(aes(shape = series, fill = series)) +
  scale_linetype_manual(values = c("Census" = "blank", "Posterior Mean" = "solid", 
                                   "95% Posterior Interval" = "dashed")) +
  scale_shape_manual(values = c("Census" = 21, "Posterior Mean" = NA, 
                                "95% Posterior Interval" = NA)) +
  scale_colour_manual(values = c("Census" = "white", "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "fill") + 
  scale_colour_manual(values = c("Census" = "black", "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "colour") + 
  # scale_alpha_manual(values = c("Census" = 1, "Posterior Mean" = 1, 
  #                               "95% Posterior Interval" = 0.5)) + 
  theme_bw() + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.35, 0.825)) +
  scale_y_continuous(name = NULL, labels = comma) + 
  scale_x_continuous("Year", expand = c(0.01, 0.01)) 



p <- p1 + p2 + plot_layout(ncol=2) + 
  plot_annotation(theme = theme(plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt")))

ggsave("Fig4.pdf",  p, width = width, height = width / 2, unit = "mm", device = cairo_pdf)

difference <- apply(mu, 1, diff)

diff_data <- data.frame(year = 1648:1850, 
                        value = c(rowSums(obs_births,na.rm=TRUE) - rowSums(obs_deaths,na.rm=TRUE) - soldiers,  
                                  apply(difference, 1, mean), apply(difference, 1, quantile, 0.025),
                                  apply(difference, 1, quantile, 0.975)) / 1000,
                        series = rep(factor(
                          c("Parish Records", "Posterior Mean", 
                            "95% Posterior Interval", "95% Posterior Interval"),
                          levels = c("Parish Records", "Posterior Mean", 
                                     "95% Posterior Interval"),
                          ordered = TRUE), each = n),
                        group = rep(factor(1:4), each = n))

p <- ggplot(diff_data, aes(x = year, y = value)) + 
  geom_line(aes(group = group, colour = series, linetype = series)) +
  scale_y_continuous("Population Difference (in thousands)", minor_breaks = NULL,
                     breaks = seq(-100, 20, by = 20), labels = seq(-100, 20, by = 20)) + 
  scale_x_continuous("Year", expand = c(0.01, 0.01)) + 
  theme_bw() + 
  scale_linetype_manual(values = c("Parish Records" = "solid", 
                                   "Posterior Mean" = "solid", 
                                   "95% Posterior Interval" = "dashed")) +
  scale_colour_manual(values = c("Parish Records" = "grey70", "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "colour") + 
  theme(text = text, axis.title = bold_text, axis.text = text) + 
  theme(plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt")) +
  theme(legend.title = element_blank(), legend.position=c(0.848, 0.171))


ggsave("Fig5.pdf", p, width = width, height = width / 2, unit = "mm", device = cairo_pdf)


prior <- readxl::read_xlsx("Aiemmat estimaatit.xlsx")[[2]]
pop_data1 <- data.frame(year = 1647:1850, 
                        value = c(NA, pop_prior, 
                                  rep(NA, 43), prior, 
                                  mu["mean",], mu["2.5%",], mu["97.5%",]) / 1000,
                        series = rep(factor(
                          c("Census", 
                            "Earlier Literature",
                            "Posterior Mean", 
                            "95% Posterior Interval", 
                            "95% Posterior Interval"),
                          levels = c("Census", 
                                     "Earlier Literature",
                                     "Posterior Mean", 
                                     "95% Posterior Interval"), ordered = TRUE), each = n + 1),
                        group = rep(factor(1:5), each = n + 1))


pop_data2 <-  data.frame(year = 1648:1850, 
                         value = 100 * c(rep(NA, 43), diff(prior)/prior[-length(prior)],
                                         rowMeans(difference/t(mu[,-(n+1)])),
                                         apply(difference/t(mu[,-(n+1)]),1,quantile, 0.025),
                                         apply(difference/t(mu[,-(n+1)]),1,quantile, 0.975)),
                         series = rep(factor(
                           c("Earlier Literature",
                             "Posterior Mean", 
                             "95% Posterior Interval", 
                             "95% Posterior Interval"),
                           levels = c("Earlier Literature",
                                      "Posterior Mean", 
                                      "95% Posterior Interval"), ordered = TRUE), each = n),
                         group = rep(factor(1:4), each = n))


p1 <- ggplot(pop_data1, aes(x = year, y = value, colour = series)) + 
  geom_line(aes(linetype = series, group = group)) +
  geom_point(aes(fill = series, shape = series)) +
  scale_y_continuous("Population (in thousands)", labels = comma) +
  scale_x_continuous("Year", expand = c(0.01, 0.01)) + 
  theme_bw() + 
  scale_linetype_manual(values = c("Census" = "blank", 
                                   "Earlier Literature" = "solid",
                                   "Posterior Mean" = "solid", 
                                   "95% Posterior Interval" = "dashed")) +
  scale_shape_manual(values = c("Census" = 21, 
                                "Earlier Literature" = NA,
                                "Posterior Mean" = NA, 
                                "95% Posterior Interval" = NA)) +
  scale_colour_manual(values = c("Census" = "white", 
                                 "Earlier Literature" = "black",
                                 "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "fill") + 
  scale_colour_manual(values = c("Census" = "black", 
                                 "Earlier Literature" = "grey70",
                                 "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "colour") + 
  theme(text = text, axis.title = bold_text, axis.text = text) +
  theme(legend.title = element_blank(), legend.position=c(0.35, 0.785)) 

xlim <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
ylim <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

p2 <- ggplot(pop_data2, aes(x = year, y = value, colour = series)) + 
  geom_line(aes(linetype = series, group = group)) +
  scale_y_continuous("Growth rate (%)") + 
  scale_x_continuous("Year", expand = c(0.01, 0.01)) + 
  scale_linetype_manual(values = c("Earlier Literature" = "solid",
                                   "Posterior Mean" = "solid", 
                                   "95% Posterior Interval" = "dashed")) +
  scale_colour_manual(values = c("Earlier Literature" = "grey70",
                                 "Posterior Mean" = "black", 
                                 "95% Posterior Interval" = "black"), aesthetics = "colour") + 
  theme_bw() +
  theme(text = text, axis.title = bold_text, axis.text = text) +
  coord_cartesian(ylim = c(-22, 4)) +
  theme(legend.title = element_blank(), legend.position=c(0.65, 0.175))


p <- p1 + p2 + plot_annotation(theme = theme(plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt")))

ggsave("Fig6.pdf", p, width = width, height = width / 2, unit = "mm", device = cairo_pdf)

