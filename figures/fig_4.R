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
  mutate(pi_pi0 = mean / mean_alpha_1, 
         sd_normalized = sd / mean_alpha_1)

data_chr1 = data %>% 
  filter(window != 'chr2') %>% 
  mutate(
    window = as.numeric(window),
    pi_pi0 = mean / mean_alpha_1,   # Normalized mean
    sd_normalized = sd / mean_alpha_1          # Normalized sd
  ) %>% 
  select(-ymin, -ymax, -sd) %>% # remove intermediate columns 
  filter(alpha == 0.004) %>%
  mutate(origin = "simulation")

data_chr1

data_chr1 = data_chr1 %>% 
  filter(as.numeric(window) >= 500000) %>% 
  mutate(
    r = 0.5*(1 - exp(-2 * (as.numeric(rho_scaled) * abs((as.numeric(window) - 500000)))))  # r
  ) 
data_chr1

###### DATA standard model ############

# Subset data for chr1
data_chr1_eqn = read_csv('../data/theoretical.csv') %>%
  filter(h == 0.5, window != 'chr2') %>% 
  mutate(alpha = as.factor(alpha),
         pi_pi0 = pi.pi0,
         rho_scaled = 1.25e-05,
         window = (as.numeric(window) + 0.5)*1000000) %>%
  filter(alpha == 0.004) %>% 
  mutate(
    r = 0.5*(1 - exp(-2 * (as.numeric(rho_scaled) * abs((as.numeric(window) - 500000)))))  # r
  ) 
data_chr1_eqn

data_chr1_eqn$origin = "standard theory"


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
alpha_values <- c(0.004) # Different alpha values

theory_data <- do.call(rbind, lapply(alpha_values, function(alpha) {
  data.frame(
    r = r_values,
    pi_pi0 = sapply(r_values, calculate_pi_pi0, alpha = alpha),
    alpha = as.factor(alpha),
    origin = "new_theory"
  ) %>% as_tibble() 
}))
theory_data

# Combine theoretical and observed data
combined_data <- bind_rows(theory_data, data_chr1,data_chr1_eqn)
combined_data

pi_pi0 = ggplot(combined_data, aes(x = r, y = pi_pi0, linetype = origin, group = origin)) +
  geom_line(linewidth = 1) +
  labs(
    x = "r",
    y = expression("average" ~ pi ~ "/" ~ pi[0] ~ "(chromosome 1)"),
  ) +
  scale_linetype_manual(
    values = c("simulation" = "solid",
               "standard theory" = "dashed",
               "new_theory" = "dotted"),
    labels = c("new_theory" = "new theory (single meiosis)",
               "standard theory" = "standard theory (rescaled r)",
               "simulation" = "simulation")
  ) +
  # Override legend aesthetics so the dashed line is more visible
  guides(linetype = guide_legend(
    override.aes = list(linewidth = 0.6,  # thickness of the lines in legend
                        keywidth = 10)   # width of the legend keys
  )) +
  theme_light() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

pi_pi0

ggsave("Fig4.png", plot = pi_pi0, width = 10, height = 6) 

