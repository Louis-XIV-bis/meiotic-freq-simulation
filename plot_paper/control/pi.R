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

csv_to_tibble = function(GR_list, subfolder, windows){
  for (GR in GR_list){ 
    data =  read_csv(paste0(subfolder,"/pi_",GR,".csv"),col_names = F)
    
    mean = apply(data, 2, mean)
    sd = apply(data, 2, sd)
    ymin = mean-sd
    ymax = mean+sd
    
    # Create a variable name based on the loop index
    GR_table_name = GR
    
    # Assign a value to the dynamically generated variable
    assign(GR_table_name, data.frame(windows,mean,ymin,ymax) %>% add_column(submodel=GR, .before=1))
  }
  
  rm(data,mean,sd,ymin,ymax,windows)
  
  pi_table = bind_rows(EGR,`100GR`) %>%
    as_tibble() %>% 
    mutate(submodel = case_when(
      submodel == "EGR" ~ "1",
      submodel == "100GR" ~ "0.01",
      TRUE ~ submodel  # Default: no change for other values
    ))
}

########################################################################
# Create tibble from csv data 

windows = seq(1,1000000,length.out=499)
GR_list = c("EGR", "100GR")

pi_table_sln = csv_to_tibble(GR_list, "sln", windows)
pi_table_sln

windows = seq(1,1000000,length.out=99)
pi_table_neutral = csv_to_tibble(GR_list, "nosln", windows)
pi_table_neutral

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(pi_table_sln$ymin, pi_table_neutral$ymin)),
                   max(c(pi_table_sln$ymax, pi_table_neutral$ymax)))

# Create color palette 
colors <- c("#1B9E77","#E7298A")

# Using the extracted color palette for both geom_line and geom_ribbon
pi_neutral = ggplot(pi_table_neutral, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "(branch length)"))+ 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
pi_neutral

pi_sln = ggplot(pi_table_sln, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel), alpha = 0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ "(branch length)"),
       color = "Meiotic frequency",
       fill = "Meiotic frequency",
       linetype = "Meiotic frequency") + 
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
pi_sln

# Arrange the plots side by side using plot_grid
combined_plot <- plot_grid(pi_neutral, pi_sln + theme(plot.margin = margin(l = 23)),
  labels = "AUTO",
  label_size = 20,
  ncol = 2,
  align = "h",  # Horizontal alignment
  rel_widths = c(1, 1.3)  # Adjust relative widths to account for the legend
)
combined_plot

ggsave("pi_control.png", plot = combined_plot, width = 16, height = 6, units = "in")

rm(list=ls())
