#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');
if (!require('sjPlot', quietly = T)) install.packages('sjPlot');
if (!require('RColorBrewer', quietly = T)) install.packages('RColorBrewer');
if (!require('dplyr', quietly = T)) install.packages('dyplr');
if (!require('readr', quietly = T)) install.packages('readr');

library(ggplot2)
library(sjPlot)
library(RColorBrewer)
library(dplyr)
library(readr)

########################################################################
# Function that upload and format the data
rm(list=ls())

input_string = "ctrl"
data = read.csv(paste0("../../../script/", input_string, "/results/t_", input_string, ".csv"))

# Add the alpha column
data = data %>% tibble() %>%
  mutate(alpha = as.factor(1 / GR)) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)

data  

########################################################################

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(0, max(data$ymax))

# Create color palette 
num_conditions <- length(unique(data$alpha))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

########################################################################

t = ggplot(data, aes(x = alpha, y = mean, color = factor(alpha))) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1, linewidth = 1.3) +
  labs(x = expression(alpha), 
       y = expression("Time to fixation"))+ 
  theme_light() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)) + 
  scale_y_continuous(limits = y_axis_limits) + 
  guides(color = "none", fill = "none") + 
  scale_color_manual(values = colors)  # Applying the predefined colors

ggsave(paste0("FigS1.png"), plot = t, width = 10, height = 6) 
rm(list=ls())