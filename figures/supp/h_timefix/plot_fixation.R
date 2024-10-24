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

library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(sjPlot)
library(RColorBrewer)
library(cowplot)

########################################################################
# Function that upload and format the data

rm(list=ls())

# Subset data for chr1
data_t = read_csv('../../../data/t_merged.csv') %>%
  filter(h == 0.5 & s == 0.05 & rho != '0,00000005') %>% 
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)
data_t

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(0, max(data_t$ymax))

# Create color palette 
num_conditions <- length(unique(data_t$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Assuming `colors` is a predefined vector of colors
plot_t = data_t %>% 
  ggplot(aes(x = alpha, y = mean, group = interaction(alpha, rho))) +
  geom_point(aes(shape = factor(rho), color = factor(alpha)), 
             size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = factor(alpha)), 
                width = 0.1, linewidth = 1.3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(15, 16, 17)) +  # Shapes: dot, square, triangle
  scale_color_manual(values = colors) +
  labs(x = expression(alpha), 
       y = expression("Time to fixation (generations)"), 
       shape = expression(rho), 
       color = expression(alpha)) + 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  scale_y_continuous(limits = y_axis_limits)
plot_t

ggsave("rho_timefix.png", plot = plot_t, width = 12, height = 10, units = "in")
