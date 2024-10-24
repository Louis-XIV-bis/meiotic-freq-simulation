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
data_chr1 = read_csv('../data/pi_merged.csv') %>%
  filter(rho == '5e-08' & s == 0.05) %>% 
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd) %>% 
  filter(window != 'chr2') %>% 
  mutate(window = as.numeric(window) / 1000000)
data_chr1

# Subset data for chr2
data_chr2 = read_csv('../data/pi_merged.csv') %>%
  filter(window == 'chr2') %>% 
  filter(rho == '5e-08' & s == 0.05) %>% 
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)  
data_chr2

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(-1e-5, 8e-05)

# Create color palette 
num_conditions <- length(unique(data_chr1$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
plot_low_h_chr1 = data_chr1 %>%
  filter(h == 0.2) %>% 
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "along chromosome 1")) + 
  theme_light() + 
  ggtitle("h = 0.2") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_low_h_chr1

plot_ctrl_chr1 = data_chr1 %>% 
  filter(h == 0.5) %>% 
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       color = expression(alpha),
       fill = expression(alpha),
       linetype = expression(alpha)) + 
  theme_light() + 
  ggtitle("h = 0.5") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5)) + 
  scale_y_continuous(limits = y_axis_limits)
plot_ctrl_chr1

plot_high_h_chr1 = data_chr1 %>% 
  filter(h == 0.8) %>%
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)") +
  theme_light() + 
  ggtitle("h = 0.8") +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_high_h_chr1

# Assuming `colors` is a predefined vector of colors
plot_chr2 = data_chr2 %>% 
  ggplot(aes(x = alpha, y = mean, group = interaction(alpha, h))) +
  geom_point(aes(shape = factor(h), color = factor(alpha)), 
             size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = factor(alpha)), 
                width = 0.1, linewidth = 1.3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(15, 16, 17)) +  # Shapes: dot, square, triangle
  scale_color_manual(values = colors) +
  labs(x = expression(alpha), 
       y = expression("average " ~ pi ~ "(chromosome 2)"), 
       shape = "h", 
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
  scale_y_continuous(limits = y_axis_limits) +
  guides(color = "none")  # Remove legend for `alpha`

plot_chr2 

# Extract the legend from the main plot
all_components <- get_plot_component(plot_ctrl_chr1, "guide-box", return_all = TRUE)
legend_only <- all_components[[1]]  # Adjust the index if needed

# First row: 3 plots with the same relative widths
first_row <- plot_grid(
  plot_low_h_chr1 + theme(plot.margin = margin(l = 23)),
  plot_ctrl_chr1 + theme(legend.position = "none", plot.margin = margin(l = 23)),
  plot_high_h_chr1 + theme(plot.margin = margin(l = 23)),
  ncol = 3, rel_widths = c(1.2, 1, 1), labels = c("A", "B", "C"), label_size = 20
)

# Second row: plot_chr2 spanning two columns and legend_only in the third column
second_row <- plot_grid(
  plot_chr2 + theme(plot.margin = margin(l = 23, t = 23)),
  legend_only,
  ncol = 2, rel_widths = c(2.5, 1),
  labels = c("D"), label_size = 20
)

# Combine the two rows into one plot
combined_plot <- plot_grid(
  first_row,
  second_row,
  ncol = 1, rel_heights = c(3, 3, 0.5) 
)
combined_plot
ggsave("fig6.png", plot = combined_plot, width = 16, height = 12, units = "in",bg = "white")
rm(list=ls())
