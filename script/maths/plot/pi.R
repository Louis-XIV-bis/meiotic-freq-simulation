#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

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
input_string = "maths"

# Function that upload and format the data
data = read.csv(paste0("../", "results/pi_", input_string, ".csv"))

# Add the alpha column
data = data %>% tibble() %>%
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)
data  

data_chr2 = data %>% 
  filter(window == 'chr2')
data_chr2

# Step 1: Extract the mean value for alpha = 1 in data_chr2
mean_alpha_1 <- data_chr2 %>% 
  filter(alpha == 1) %>% 
  pull(mean)

data_chr2 = data_chr2 %>% 
  mutate(mean_normalized = mean / mean_alpha_1, 
         sd_normalized = sd / mean_alpha_1)

data_chr1 = data %>% 
  filter(window != 'chr2') %>% 
  mutate(
    window = as.numeric(window),
    mean_normalized = mean / mean_alpha_1,   # Normalized mean
    sd_normalized = sd / mean_alpha_1          # Normalized sd
  ) %>% 
  select(-ymin, -ymax, -sd) # remove intermediate columns
data_chr1

data_chr1 = data_chr1 %>% 
  filter(as.numeric(window) >= 500000) %>% 
  mutate(
    r = 0.5*(1 - exp(-2 * (as.numeric(rho_scaled) * abs((as.numeric(window) - 500000)))))  # r
  )  
data_chr1

########################################################################

# Create color palette 
num_conditions <- length(unique(data$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
pi_chr1 = ggplot(data_chr1, aes(x = r, y = mean_normalized, group = alpha)) +
  geom_ribbon(aes(ymin = mean_normalized - sd_normalized, ymax = mean_normalized + sd_normalized, fill = alpha), alpha=0.2) +
  geom_line(aes(color = alpha)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "recombination frequency",
       y = expression(pi ~ "/" ~ pi ~ "chr2 (" * alpha *  "=1)"),
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
  )
pi_chr1
ggsave(paste0("pi_chr1_", input_string, ".png"), plot = pi_chr1, width = 10, height = 6) 

pi_chr2 = ggplot(data_chr2, aes(x = alpha, y = mean_normalized, color = factor(alpha))) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = mean_normalized- sd_normalized, ymax = mean_normalized + sd_normalized), width = 0.1, linewidth = 1.3) +
  labs(x = expression(alpha), 
       y = expression(pi ~ "/" ~ pi ~ "chr2 (" * alpha *  "=1)")) + 
  theme_light() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)) + 
  guides(color = "none", fill = "none") + 
  scale_color_manual(values = colors)  # Applying the predefined colors
pi_chr2
ggsave(paste0("pi_chr2_", input_string, ".png"), plot = pi_chr2, width = 10, height = 6) 