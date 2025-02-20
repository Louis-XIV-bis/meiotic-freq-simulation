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

# Subset data for ctrl and neutral condition and for alpha = 0.01 and 1
data_var = read_csv('../../../data/pi_merged.csv') %>%
  filter((s == 0.05 & h == 0.5 & rho == '5e-08')) %>% 
  filter(window != 'chr2') %>%   
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(alpha %in% c(0.01, 1)) %>% 
  mutate(ymin = mean - sd, ymax = mean + sd) %>% 
  mutate(window = as.numeric(window) / 1000000)
data_var

data_fix = read_csv('pi_nogaussian_alpha.csv') %>%
  filter((s == 0.05 & h == 0.5 & rho == '5e-08') | (s == 0)) %>% 
  filter(window != 'chr2') %>%   
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(alpha %in% c(0.01, 1)) %>% 
  mutate(ymin = mean - sd, ymax = mean + sd) %>% 
  mutate(window = as.numeric(window) / 1000000)
data_fix

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(-1e-5, 8e-05)

# Create color palette 
colors <- c("#1B9E77","#E7298A")
pi_fix = data_fix %>%
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
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
    axis.text.y = element_text(size = 18)
    ) +
  guides(color = "none", fill = "none") +
  scale_y_continuous(limits = y_axis_limits)
pi_fix


pi_var = data_var %>%
   ggplot(aes(x = window, y = mean, group = alpha)) +
   geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
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
       axis.title.y = element_blank(),
      legend.title = element_text(size = 22),
       legend.text = element_text(size = 22),
       axis.text.x = element_text(size = 18),
       axis.text.y = element_blank()
     ) + 
   scale_y_continuous(limits = y_axis_limits)
pi_var

# Arrange the plots side by side using plot_grid
combined_plot <- plot_grid(pi_fix, pi_var + theme(plot.margin = margin(l = 23)),
                           labels = "AUTO",
                           label_size = 20,
                           ncol = 2,
                           align = "h",  # Horizontal alignment
                           rel_widths = c(1, 1)  # Adjust relative widths to account for the legend
)
combined_plot

# Adding labels 
plot = ggdraw() +
  draw_plot(combined_plot, width = 1, height = 1) +
  draw_text("Fixed intervals", x = 0.29, y = 0.92, size = 18, color = "black") +
  draw_text("Variable intervals", x = 0.71, y = 0.92, size = 18, color = "black")
plot

ggsave("figS2.png", plot = plot, width = 16, height = 6, units = "in")

rm(list=ls())
