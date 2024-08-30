#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');
if (!require('sjPlot', quietly = T)) install.packages('sjPlot');
if (!require('RColorBrewer', quietly = T)) install.packages('RColorBrewer');
if (!require('dplyr', quietly = T)) install.packages('dyplr');
if (!require('readr', quietly = T)) install.packages('readr');
if (!require('dunn.test', quietly = T)) install.packages('cowplot');

library(ggplot2)
library(sjPlot)
library(RColorBrewer)
library(dplyr)
library(readr)
library(cowplot)
library(dunn.test)

########################################################################
rm(list=ls())

t_100GR_low <- read.table("fix_100GR_1e-8.txt", header = FALSE, col.names = "time")
t_100GR_low$m = "0.01"
t_10GR_low <- read.table("fix_10GR_1e-8.txt", header = FALSE, col.names = "time")
t_10GR_low$m = "0.1"
t_50GR_low <- read.table("fix_50GR_1e-8.txt", header = FALSE, col.names = "time")
t_50GR_low$m = "0.05"
t_EGR_low <- read.table("fix_EGR_1e-8.txt", header = FALSE, col.names = "time")
t_EGR_low$m = "1"
t_low <- rbind(t_100GR_low, t_50GR_low, t_10GR_low, t_EGR_low)

t_100GR_mid <- read.table("fix_100GR_5e-8.txt", header = FALSE, col.names = "time")
t_100GR_mid$m = "0.01"
t_10GR_mid <- read.table("fix_10GR_5e-8.txt", header = FALSE, col.names = "time")
t_10GR_mid$m = "0.1"
t_50GR_mid <- read.table("fix_50GR_5e-8.txt", header = FALSE, col.names = "time")
t_50GR_mid$m = "0.05"
t_EGR_mid <- read.table("fix_EGR_5e-8.txt", header = FALSE, col.names = "time")
t_EGR_mid$m = "1"
t_mid <- rbind(t_100GR_mid, t_50GR_mid, t_10GR_mid, t_EGR_mid)

t_100GR_high <- read.table("fix_100GR_1e-7.txt", header = FALSE, col.names = "time")
t_100GR_high$m = "0.01"
t_10GR_high <- read.table("fix_10GR_1e-7.txt", header = FALSE, col.names = "time")
t_10GR_high$m = "0.1"
t_50GR_high <- read.table("fix_50GR_1e-7.txt", header = FALSE, col.names = "time")
t_50GR_high$m = "0.05"
t_EGR_high <- read.table("fix_EGR_1e-7.txt", header = FALSE, col.names = "time")
t_EGR_high$m = "1"
t_high <- rbind(t_100GR_high, t_50GR_high, t_10GR_high, t_EGR_high)

########################################################################

t_low_plot = ggplot(t_low, aes(x = m, y = time - 2000, fill = m)) +
  geom_boxplot() +
  labs(x = expression(m),
       y = "Fixation time (generations)",
       color=expression(m),
       fill=expression(m),
       linetype=expression(m)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  ggtitle(expression(paste(rho ~ " = 1e-8"))) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size=16), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  guides(color = "none", fill = "none") 
t_low_plot

t_mid_plot = ggplot(t_mid, aes(x = m, y = time - 2000, fill = m)) +
  geom_boxplot() +
  labs(x = expression(m),
       y = "Fixation time (generations)",
       color=expression(m),
       fill=expression(m),
       linetype=expression(m)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  ggtitle(expression(paste(rho ~ " = 5e-8 (control)"))) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5), 
    legend.position = "top"
  ) 
t_mid_plot

t_high_plot = ggplot(t_high, aes(x = m, y = time - 2000, fill = m)) +
  geom_boxplot() +
  labs(x = expression(m),
       y = "Fixation time (generations)",
       color=expression(m),
       fill=expression(m),
       linetype=expression(m)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  ggtitle(expression(paste(rho ~ " = 1e-7"))) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5), 
  ) +
  guides(color = "none", fill = "none") 
t_high_plot

# Extract the legend from the main plot
legend_only <- get_legend(t_mid_plot)

# Create the grid of plots using plot_grid
combined_plot <- plot_grid(
  t_low_plot + theme(plot.margin = margin(l = 23)),
  t_mid_plot + theme(legend.position = "none", plot.margin = margin(l = 23)),
  t_high_plot + theme(plot.margin = margin(l = 23)), 
  NULL, legend_only, NULL,
  rel_heights = c(3,0.5),
  rel_widths = c(1.2,1,1),
  ncol = 3, nrow = 2, labels = c("A","B","C"), label_size = 20
)
combined_plot
png('rhovar_timefix.png',width=16,height=8,units="in",bg = "white", res=300)
combined_plot
dev.off()


# Statistical tests 

kruskal.test(time -2000 ~ factor(m), data = t_low)
# p-value = 0.4169 ==> No significiant != between some groups for low rho

kruskal.test(time -2000 ~ factor(m), data = t_mid)
# p-value = 0.01356 ==> significiant != between some groups for ctrl rho : 0.01 and 1 
dunn.test(t_mid$time - 2000, factor(t_mid$m), method = "bonferroni")

kruskal.test(time -2000 ~ factor(m), data = t_high)
# p-value = 0.1473 ==> No significiant != between some groups for high rho