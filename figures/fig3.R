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

##### SIMU DATA ######
# Subset data for chr1
data_chr1 = read_csv('../data/pi_merged.csv') %>%
  filter(h == 0.5 & rho == '5e-08' & s > 0) %>% 
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(window != 'chr2')

# Subset data for chr2
data_chr2 = read_csv('../data/pi_merged.csv') %>%
  filter(window == 'chr2') %>% 
  filter(h == 0.5 & rho == '5e-08' & s > 0) %>% 
  mutate(alpha = as.factor(1 / GR))
data_chr2

# Compute pi/pi0(chr2) for each condition 
mean_s002 <- data_chr2 %>% filter(alpha == 1, s == 0.02) %>% pull(mean)
mean_s005 <- data_chr2 %>% filter(alpha == 1, s == 0.05) %>% pull(mean)
mean_s01 <- data_chr2 %>% filter(alpha == 1, s == 0.1) %>% pull(mean)

data_chr2 = data_chr2 %>% 
  mutate(mean_normalized = case_when(
    s == 0.02 ~ mean / mean_s002,  # For s = 0.02, divide by mean_s1
    s == 0.05 ~ mean / mean_s005,  # For s = 0.05, divide by mean_s2
    s == 0.1  ~ mean / mean_s01   # For s = 0.1, divide by mean_s3
  ),
  sd_normalized = case_when(
    s == 0.02 ~ sd / mean_s002,  # For s = 0.02, divide by mean_s1
    s == 0.05 ~ sd / mean_s005,  # For s = 0.05, divide by mean_s2
    s == 0.1  ~ sd / mean_s01   # For s = 0.1, divide by mean_s3
  ))
data_chr2

data_chr1 <- data_chr1 %>%
  mutate(
    mean_normalized = case_when(
      s == 0.02 ~ mean / mean_s002,  # For s = 0.02, divide by mean_s1
      s == 0.05 ~ mean / mean_s005,  # For s = 0.05, divide by mean_s2
      s == 0.1  ~ mean / mean_s01   # For s = 0.1, divide by mean_s3
    ),
    sd_normalized = case_when(
      s == 0.02 ~ sd / mean_s002,  # For s = 0.02, divide by mean_s1
      s == 0.05 ~ sd / mean_s005,  # For s = 0.05, divide by mean_s2
      s == 0.1  ~ sd / mean_s01   # For s = 0.1, divide by mean_s3
    ) 
  )
data_chr1

data_chr1$origin = "simulation"
data_chr1 = data_chr1 %>%
  mutate(window = as.numeric(window) / 1000000) %>%
  filter(window >= 0.5)

data_chr2$origin = "simulation"

###### THEORICAL DATA ####### 

# Subset data for chr1
data_chr1_eqn = read_csv('../data/theoretical.csv') %>%
  filter(h == 0.5, window != 'chr2') %>% 
  mutate(alpha = as.factor(alpha),
         mean_normalized = pi.pi0,
         window = as.numeric(window) + 0.5) %>%
  filter(alpha != 0.004)
data_chr1_eqn

# Subset data for chr2
data_chr2_eqn = read_csv('../data/theoretical.csv')  %>%
  filter(h == 0.5) %>% 
  mutate(alpha = as.factor(alpha),
         mean_normalized = pi.pi0 ) %>% 
  filter(window == 'chr2')
data_chr2_eqn

data_chr1_eqn$origin = "standard theory"
data_chr2_eqn$origin = "standard theory"

############# Merge dataset 
data_combined_chr1 <- bind_rows(
  data_chr1,
  data_chr1_eqn
)

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(-0.1, 1.30)

# Create color palette 
num_conditions <- length(unique(data_chr1_eqn$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Plot
plot_low_s_chr1 = data_combined_chr1 %>%
  filter(s == 0.02) %>%
  ggplot(aes(x = window, y = mean_normalized, group = interaction(alpha, origin), 
             color = alpha, linetype = origin)) +
  geom_line() +  
  scale_color_manual(values = colors, name = expression(alpha)) + 
  scale_linetype_manual(values = c("simulation" = "solid", "standard theory" = "dashed"), name = NULL, guide = NULL) +
  labs(x = "Sequence position (Mbp)", y = expression(pi ~ "/" ~ pi[0] ~ "(chromosome 1)")) + 
  theme_light() + 
  ggtitle("s = 0.02") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_low_s_chr1

plot_ctrl_chr1 = data_combined_chr1 %>%
  filter(s == 0.05) %>%
  ggplot(aes(x = window, y = mean_normalized, group = interaction(alpha, origin), 
             color = alpha, linetype = origin)) +
  geom_line() +  
  scale_color_manual(values = colors, name = expression(alpha)) + 
  scale_linetype_manual(values = c("simulation" = "solid", "standard theory" = "dashed"), name = NULL, guide = NULL) +
  labs(x = "Sequence position (Mbp)", y = expression(pi ~ "/" ~ pi[0] ~ "(chromosome 1)")) + 
  theme_light() + 
  ggtitle("s = 0.05") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y =element_blank(),,
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_ctrl_chr1

plot_high_s_chr1 = data_combined_chr1 %>%
  filter(s == 0.1) %>%
  ggplot(aes(x = window, y = mean_normalized, group = interaction(alpha, origin), 
             color = alpha, linetype = origin)) +
  geom_line() +  
  scale_color_manual(values = colors, name = expression(alpha)) + 
  scale_linetype_manual(
    values = c("simulation" = "solid", "standard theory" = "dashed"), 
    labels = c("simulation" = "simulation", "standard theory" = "standard\ntheory"), 
    name = NULL
  ) +
  labs(x = "Sequence position (Mbp)", y = expression(pi ~ "/" ~ pi[0] ~ "(chromosome 1)")) + 
  theme_light() + 
  ggtitle("s = 0.1") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.height = unit(30, "pt")  # Adjust spacing inside the linetype legend
  ) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) + 
  scale_y_continuous(limits = y_axis_limits)

plot_high_s_chr1


# First row: 3 plots with the same relative widths
plot <- plot_grid(
  plot_low_s_chr1 + theme(plot.margin = margin(l = 23)),
  plot_ctrl_chr1 + theme(plot.margin = margin(l = 23)),
  plot_high_s_chr1 + theme(plot.margin = margin(l = 23)),
  ncol = 3, rel_widths = c(1.17, 1, 1.4), labels = c("A", "B", "C"), label_size = 20
)
plot

ggsave("fig3.png", plot = plot, width = 16, height = 10, units = "in",bg = "white")
rm(list=ls())
