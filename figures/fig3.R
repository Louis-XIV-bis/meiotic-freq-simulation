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
if (!require('patchwork', quietly = T)) install.packages('patchwork');

library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(sjPlot)
library(RColorBrewer)
library(cowplot)
library(patchwork)

########################################################################
# Function that upload and format the data

rm(list=ls())

data_chr1 = read_csv('../data/pi_merged.csv') %>%
  filter(h == 0.5 & rho == '5e-08' & s > 0) %>% 
  mutate(alpha = as.factor(1 / GR)) %>%
  filter(window != 'chr2') %>%
  mutate(window = as.numeric(window) / 1000000)
data_chr1

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
data_chr2 = data_chr2 %>% select(-c(GR, rho, rho_scaled, mean, sd, sd_normalized))
data_chr2$origin <- "simulation"
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
data_chr1 = data_chr1 %>% select(-c(GR, rho, rho_scaled, mean, sd, sd_normalized))
data_chr1$origin <- "simulation"
data_chr1

# Read the data from file
data2 <- read_csv("../data/fig3_theoretical.csv") %>% as_tibble()
data2_chr1 = data2 %>%
  filter(window != 'chr2') %>% as_tibble() %>%
  mutate(h = 0.5, mean_normalized = pi.pi0) %>% 
  select(-pi.pi0)
data2_chr1$window = as.numeric(data2_chr1$window) 
data2_chr1$origin <- "theory"

data2_chr1 = data2_chr1 %>% 
  mutate(window = window + 0.5)
data2_chr1

data2_chr2 = data2 %>%
  filter(window == 'chr2') %>% as_tibble() %>%
  mutate(h = 0.5, mean_normalized = pi.pi0) %>% 
  select(-pi.pi0)
data2_chr2$origin <- "theory"
data2_chr2

############## 

chr1 = rbind(data_chr1, data2_chr1) 
chr1 = chr1 %>%
  filter(!is.na(mean_normalized), window > 0.5) 
chr1

data_chr2 = data_chr2 %>% mutate(window = 1.1)
data2_chr2 = data2_chr2 %>% mutate(window = 1.2)
chr2 = rbind(data_chr2, data2_chr2)
chr2 = chr2 %>%
  filter(!is.na(mean_normalized))
chr2

###########################  
# Combine datasets
combined_data <- bind_rows(chr1, chr2)
combined_data = combined_data %>%
  filter(alpha == 1)
combined_data

theo_vs_simu <- ggplot(combined_data, aes(x = window, y = mean_normalized, color = origin, shape = origin, group = origin)) +
  # Add grey band for chr2 (window 1.1 - 1.2)
  annotate("rect", xmin = 1.05, xmax = 1.25, ymin = -Inf, ymax = Inf, 
           fill = "grey80", alpha = 0.5) +  
  geom_line(data = combined_data %>% filter(!window %in% c(1.1, 1.2)), size = 1) +  # Line plot for most data
  geom_point(data = combined_data %>% filter(window %in% c(1.1, 1.2)), size = 4) +  # Points for special windows
  facet_grid(~s, scales = "free_x", labeller = labeller(s = function(x) paste0("s = ", x))) +  # Custom facet labels
  scale_x_continuous(
    breaks = c(seq(0, 1, by = 0.2), 1.15),  
    labels = c(as.character(seq(0, 1, by = 0.2)), "chr2")  
  ) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "/" ~ pi[0]),
       color = "Origin",
       shape = "Origin") +  # Make sure both aesthetics have the same title
  theme_minimal() +
  theme(
    legend.title = element_blank(),  
    text = element_text(size = 16),  
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 14),  
    strip.text = element_text(size = 16, face = "bold"),  
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = c(16, 17))),  # Ensure different shapes are in the legend
    shape = guide_legend(override.aes = list(color = c("#F8766D", "#00BFC4")))  # Ensure correct colors for shapes
  )

theo_vs_simu

#############
data_plot_alpha_1 = data_chr2 %>%
  filter(alpha != 1)

data_plot_alpha_2 = data2_chr2 %>%
  filter(alpha != 1)

data_plot_alpha = rbind(data_plot_alpha_1, data_plot_alpha_2)
data_plot_alpha

theo_simu_chr2 = ggplot(data_plot_alpha, aes(x = alpha, y = mean_normalized, color = origin)) +
  geom_point(aes(shape = origin), size = 5, alpha = 0.8) +  # Different shapes for origin
  facet_grid(~s, labeller = labeller(s = function(x) paste0("s = ", x))) +  # Facet per s
  labs(x = expression(alpha), y = expression(pi ~ "/" ~ pi[0] ~ "(chromosome 2)"), 
       color = "Origin", shape = "Origin") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    
    # Add separation lines between facets
    panel.spacing = unit(0.8, "lines"),  # Adjust spacing between facets
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)  # Thin black border around each facet
  )

theo_simu_chr2

# Label the plots as A and B
final_plot <- (theo_vs_simu / theo_simu_chr2) + 
  plot_annotation(tag_levels = "A")  # This automatically adds "A" and "B" labels

final_plot

ggsave('fig3.png', plot = final_plot, bg='white')
