#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');
if (!require('readr', quietly = T)) install.packages('readr');
if (!require('tibble', quietly = T)) install.packages('tibble');
if (!require('dplyr', quietly = T)) install.packages('dplyr');
if (!require('sjPlot', quietly = T)) install.packages('sjPlot');
if (!require('RColorBrewer', quietly = T)) install.packages('RColorBrewer');

library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(sjPlot)
library(RColorBrewer)

########################################################################

rm(list=ls())
data = read.csv('../results/pi_neutral.csv')

# Add the alpha column
data = data %>% tibble() %>%
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)
data  
  
data_chr1 = data %>% 
  filter(alpha %in% c(0.01, 1)) %>% 
  filter(window != 'chr2') %>% 
  mutate(window = as.numeric(window) / 1000000)

data_chr2 = data %>% 
  filter(alpha %in% c(0.01, 1)) %>% 
  filter(window == 'chr2') 

########################################################################

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(data$ymin, data$ymin)),
                   max(c(data$ymax, data$ymax)))

# Create color palette 
num_conditions <- length(unique(data$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
pi_chr1 = ggplot(data_chr1, aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha=0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi),
       color = expression(alpha),
       fill = expression(alpha),
       linetype = expression(alpha)) +
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22)
  ) +
  scale_y_continuous(limits = y_axis_limits)
pi_chr1

ggsave("pi_chr1.png") 

pi_chr2 = ggplot(data_chr2, aes(x = alpha, y = mean, color = factor(alpha))) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1, size = 1.3) +
  labs(x = expression(alpha), 
       y = expression("average " ~ pi ~ "(chromosome 2)"))+ 
  theme_light() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)) + 
  scale_y_continuous(limits = y_axis_limits) + 
  guides(color = "none", fill = "none") + 
  scale_color_manual(values = colors)  # Applying the predefined colors
pi_chr2

ggsave("pi_chr2.png") 

rm(list=ls())

