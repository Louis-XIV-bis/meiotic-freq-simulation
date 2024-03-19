#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@universite-paris-saclay.fr
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');
if (!require('readr', quietly = T)) install.packages('readr');
if (!require('tibble', quietly = T)) install.packages('tibble');
if (!require('dplyr', quietly = T)) install.packages('dplyr');
if (!require('sjPlot', quietly = T)) install.packages('sjPlot');
if (!require('RColorBrewer', quietly = T)) install.packages('RColorBrewer');
if (!require('cowplot', quietly = T)) install.packages('cowplot');
if (!require('rstatix', quietly = T)) install.packages('rstatix');

library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(sjPlot)
library(cowplot)
library(rstatix)

rm(list=ls())

# Read the CSV files into data frames
col_names <- c("submodel", "pi")  # Specify your column names here

nosln_100GR <- read.csv("nosln/pi_gw_100GR.csv", header = TRUE)
nosln_EGR <- read.csv("nosln/pi_gw_EGR.csv", header = TRUE)
colnames(nosln_100GR) <- col_names
colnames(nosln_EGR) <- col_names

merged_nosln <- rbind(nosln_100GR, nosln_EGR)

sln_100GR <- read.csv("sln/pi_gw_100GR.csv", header = TRUE)
sln_EGR <- read.csv("sln/pi_gw_EGR.csv", header = TRUE)
colnames(sln_100GR) <- col_names
colnames(sln_EGR) <- col_names

merged_sln <- rbind(sln_100GR, nosln_EGR)

########################################################################
# Statistical tests 
shapiro.test(nosln_100GR$pi) 
# p-value = 0.695 => gaussian distribution
hist(as.vector(nosln_100GR$pi), probability = TRUE)

shapiro.test(nosln_EGR$pi) 
# p-value = 0.2165 => gaussian distribution
hist(as.vector(nosln_EGR$pi), probability = TRUE)

t.test(nosln_100GR$pi, nosln_EGR$pi)
# p-value = 0.0008916 ==> difference
summary(nosln_100GR$pi)
summary(nosln_EGR$pi)

shapiro.test(sln_100GR$pi) 
# p-value = 0.3689 => gaussian distribution
hist(as.vector(sln_100GR$pi), probability = TRUE)

shapiro.test(sln_EGR$pi) 
# p-value = 0.7173 => gaussian distribution
hist(as.vector(sln_EGR$pi), probability = TRUE)

t.test(sln_100GR$pi, sln_EGR$pi, alternative="less")
# p-value = < 2.2e-16 ==> 100GR < EGR

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(merged_sln$pi, merged_nosln$pi)),
                   max(c(merged_sln$pi, merged_nosln$pi)))

# Create color palette 
colors <- c("#1B9E77","#E7298A")

pi_nosln = ggplot(merged_nosln, aes(x = factor(submodel), y = pi, group=submodel, fill=factor(submodel))) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  labs(x = "Meiotic frequency",
       y = expression(pi ~ "(branch length)"))+ 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  ggtitle("Without selection") +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
pi_nosln

pi_sln = ggplot(merged_sln, aes(x = factor(submodel), y = pi, group=submodel, fill=factor(submodel))) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  labs(x = "Meiotic frequency",
       y = expression(pi ~ "(branch length)"))+
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  ggtitle("With selection") +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
pi_sln

# Arrange the plots side by side using plot_grid
combined_plot <- plot_grid(pi_nosln, pi_sln + theme(plot.margin = margin(l = 23)),
                           labels = "AUTO",
                           label_size = 20,
                           ncol = 2,
                           align = "h"  # Horizontal alignment
)
combined_plot

ggsave("ctrl_wgpi.png", plot = combined_plot, width = 16, height = 6, units = "in")

rm(list=ls())

