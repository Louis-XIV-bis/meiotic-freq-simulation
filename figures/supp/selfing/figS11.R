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

rm(list=ls())

# Upload all dataset (different selfing rates)

s0p1 = read.csv(paste0("../../../script/selfing/0p1/results/pi_0p1.csv"))
s0p1$selfRate = 0.1

s0p5 = read.csv(paste0("../../../script/selfing/0p5/results/pi_0p5.csv"))
s0p5$selfRate = 0.5

s0p9 = read.csv(paste0("../../../script/selfing/0p9/results/pi_0p9.csv"))
s0p9$selfRate = 0.9

data = bind_rows(s0p1, s0p5, s0p9)
data

# Add the alpha column
data = data %>%
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd) %>%
  filter(window != 'chr2') %>% 
  mutate(window = as.numeric(window) / 1000000)
data  

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(-1e-5, 8e-05)

# Create color palette 
num_conditions <- length(unique(data$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
plot_low_self = data %>%
  filter(selfRate == 0.1) %>% 
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "along chromosome 1")) +
  theme_light() + 
  ggtitle("selfing rate = 0.1") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") +
  scale_y_continuous(limits = y_axis_limits)
plot_low_self

plot_mid_self = data %>% 
  filter(selfRate == 0.5) %>% 
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
  ggtitle("selfing rate = 0.5") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "bottom") + 
  scale_y_continuous(limits = y_axis_limits)
plot_mid_self

plot_high_self = data %>% 
  filter(selfRate == 0.9) %>%
  ggplot(aes(x = window, y = mean, group = alpha)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = alpha), alpha = 0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)") +
  theme_light() + 
  ggtitle("selfing rate = 0.9") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_high_self


# Extract the legend from the main plot
all_components <- get_plot_component(plot_mid_self, "guide-box", return_all = TRUE)
legend_only <- all_components[[3]]  # Adjust the index if needed

# First row: 3 plots with the same relative widths
first_row <- plot_grid(
  plot_low_self + theme(plot.margin = margin(l = 23)),
  plot_mid_self + theme(legend.position = "none", plot.margin = margin(l = 23)),
  plot_high_self + theme(plot.margin = margin(l = 23)),
  ncol = 3, rel_widths = c(1.2, 1, 1), labels = c("A", "B", "C"), label_size = 20
)

# Second row: plot_chr2 spanning two columns and legend_only in the third column
second_row <- plot_grid(
  legend_only,
  ncol = 1,
  label_size = 20
)
second_row

# Combine the two rows into one plot
combined_plot <- plot_grid(
  first_row,
  second_row,
  ncol = 1, rel_heights = c(3, 0.5) 
)
combined_plot

ggsave("FigS11.png", plot = combined_plot, width = 20, height = 12, units = "in",bg = "white")
rm(list=ls())