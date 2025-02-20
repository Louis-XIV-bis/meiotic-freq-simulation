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
data = read.csv(paste0("../script/", input_string, "/results/pi_", input_string, ".csv"))

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

###### EQUATION 10 CHARLESWORTH ########

# Define constants
N <- 1000
s <- 0.1
h <- 0.5

# Define the function to calculate pi_pi0
calculate_pi_pi0 <- function(r, alpha) {
  gamma <- 2 * N * h * s
  Ts1 <- 2 * (1 / gamma) * (log(log(gamma) + 0.5772)) # Ts1 calculation
  ts1 <- Ts1 * 2 * N
  Ts2 <- (1 / (2 * N * (1-h) * s)) * log(4 * N * (1 -h) * s) # Ts2 calculation
  Tsm2 <- (1 / (2 * N * (1-h) * s)) * (log(log(2 * N * (1 - h) * s)) + log(log(2 * N * h * s)) - 0.5 * (1 / alpha) * h * s) # Tsm2 calculation
  
  # Equation for pi_pi0
  2 * r * (1 - r) + 1 / (2 * N * alpha) + alpha * ts1 * Tsm2 + (1 - alpha * ts1) * (Ts1 + Ts2)
}

# Generate data for plotting
r_values <- seq(0, 0.5, length.out = 50) # r values from 0 to 0.5
alpha_values <- c(0.004, 0.01) # Different alpha values

theory_data <- do.call(rbind, lapply(alpha_values, function(alpha) {
  data.frame(
    r = r_values,
    pi_pi0 = sapply(r_values, calculate_pi_pi0, alpha = alpha),
    alpha = as.factor(alpha),
    dataset = "Theoretical"
  ) %>% as_tibble()
}))

# Filter observed data to match alpha values of interest
observed_data <- data_chr1 %>% 
  filter(alpha %in% c(0.01, 0.004)) %>% 
  mutate(dataset = "Simulated", pi_pi0 = mean_normalized) # Rename mean_normalized to pi_pi0

# Combine theoretical and observed data
combined_data <- bind_rows(theory_data, observed_data)

# Plot the data
pi_pi0 = ggplot(combined_data, aes(x = r, y = pi_pi0, color = alpha, linetype = dataset, group = interaction(alpha, dataset))) +
  geom_line(linewidth = 1) +
  labs(
    x = "r",
    y = expression("average" ~ pi ~ "/" ~ pi[0] ~ "(chromosome 1)"),
    color = expression(alpha),
    linetype = "Dataset"
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 16),          # Increase overall text size
    axis.title = element_text(size = 18),    # Increase axis title size
    axis.text = element_text(size = 16),     # Increase axis tick label size
    legend.text = element_text(size = 16),   # Increase legend text size
    legend.title = element_text(size = 18)   # Increase legend title size
  )
pi_pi0

ggsave("fig4.png", plot = pi_pi0, width = 10, height = 6) 

