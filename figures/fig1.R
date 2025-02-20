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

# Neutral dataset 
data_neutral = read_csv('../data/pi_merged.csv') %>%
  filter(s == 0) %>% 
  filter(window != 'chr2') %>%   
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(alpha %in% c(0.01, 1)) %>% 
  mutate(ymin = mean - sd, ymax = mean + sd) %>% 
  mutate(window = as.numeric(window) / 1000000)
data_neutral

# Chromosome 1 w/ sweep
data_chr1 = read_csv('../data/pi_merged.csv') %>%
  filter((s == 0.05 & h == 0.5 & rho == '5e-08')) %>% 
  filter(window != 'chr2') %>%   
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(alpha %in% c(0.01, 1)) %>% 
  mutate(ymin = mean - sd, ymax = mean + sd) %>% 
  mutate(window = as.numeric(window) / 1000000) 
data_chr1

# Full chromosome 2 positions
data_chr2 = read_csv('../data/pi_full_chr2_ctrl.csv') %>%
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(alpha %in% c(0.01, 1)) %>% 
  mutate(ymin = mean - sd, ymax = mean + sd) %>% 
  mutate(window = as.numeric(window) / 1000000)
data_chr2

data_sln = bind_rows(data_chr1, data_chr2)

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(-1e-5, 8e-05)

# Create color palette 
colors <- c("#1B9E77","#E7298A")

# Using the extracted color palette for both geom_line and geom_ribbon
pi_neutral = data_neutral %>%
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha=0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "along chromosome 1 (x 10"^-5*")")) +
theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 15)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(
    limits = y_axis_limits,
    labels = scientific_formatter
  )
pi_neutral

# Using the extracted color palette for both geom_line and geom_ribbon
pi_sln = data_sln %>%
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha=0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "along the genome"),
       color = expression(alpha),
       fill = expression(alpha),
       linetype = expression(alpha)) + 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(limits = y_axis_limits) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "black")
pi_sln

# Adding labels 
plot_sln = ggdraw() +
  draw_plot(pi_sln, width = 1, height = 1) +
  draw_text("Chromosome 1", x = 0.24, y = 0.92, size = 18, color = "black") +
  draw_text("Chromosome 2", x = 0.61, y = 0.92, size = 18, color = "black")
plot_sln

combined_plot <- plot_grid(
  pi_neutral,          # First plot
  plot_sln,            # Second plot
  ncol = 2,           
  rel_widths = c(1.3, 2.6), 
  labels = c("A", "B"), label_size = 20
  
)
combined_plot


png('fig1.png', width=14, height=6, units="in",bg = "white", res=600)
combined_plot
dev.off()
