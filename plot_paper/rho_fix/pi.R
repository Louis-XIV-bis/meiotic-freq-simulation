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
  
  pi_table = bind_rows(EGR,`10GR`,`50GR`,`100GR`) %>%
    as_tibble() %>% 
    mutate(submodel = case_when(
      submodel == "EGR" ~ "1",
      submodel == "10GR" ~ "0.1",
      submodel == "50GR" ~ "0.05",
      submodel == "100GR" ~ "0.01",
      TRUE ~ submodel  # Default: no change for other values
    ))
}
########################################################################
# Create tibble from csv data 

windows = seq(1,2000000,length.out=499)

GR_list = c("EGR", "10GR", "50GR", "100GR")
pi_table_fix = csv_to_tibble(GR_list, "rhofixe", windows) %>% filter(windows < 1000000)
pi_table_fix

pi_table_var = csv_to_tibble(GR_list, "ctrl", windows) %>% filter(windows < 1000000)
pi_table_var

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(pi_table_fix$ymin, pi_table_var$ymin)),
                   max(c(pi_table_fix$ymax, pi_table_var$ymax)))

# Create color palette 
num_conditions <- length(unique(pi_table_var$submodel))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
pi_var = ggplot(pi_table_var, aes(x = windows/1000000, y = mean, group = submodel)) +
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
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  ggtitle(expression(rho ~ "= 5e-8 (control)")) +
  scale_y_continuous(limits = y_axis_limits)
pi_var

pi_fix = ggplot(pi_table_fix, aes(x = windows/1000000, y = mean, group = submodel)) +
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  ggtitle(expression(paste(rho[m]," = 5e-8"))) +
  scale_y_continuous(limits = y_axis_limits)
pi_fix

# Arrange the plots side by side using plot_grid
combined_plot_chr1 <- plot_grid(pi_var, pi_fix + theme(plot.margin = margin(l = 23)),
                           labels = "AUTO",
                           label_size = 20,
                           ncol = 2,
                           align = "h",  # Horizontal alignment
                           rel_widths = c(1, 1.3)  # Adjust relative widths to account for the legend
)
combined_plot_chr1

ggsave("pi_rhofixvar.png", plot = combined_plot_chr1, width = 16, height = 6, units = "in")

pi_chr2_var = csv_to_tibble(GR_list, "rhovar", windows) %>% filter(windows > 1000000)
pi_chr2_var
pi_chr2_fix = csv_to_tibble(GR_list, "rhofixe", windows) %>% filter(windows > 1000000)
pi_chr2_fix

pi_var_chr2 = ggplot(pi_chr2_var, aes(x = submodel, y = (mean - 2000), group = submodel, color = factor(submodel), fill = factor(submodel))) +
  geom_boxplot() +
  labs(x = "Meiotic frequency (m)",
       y = "Fixation time (generations)",
       title = expression(paste(rho," = 5e-8 (chromosome 2)")),
       color = "Meiotic frequency (m)",
       fill = "Meiotic frequency (m)") +
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  guides(color = "none", fill = "none") + 
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")
pi_var_chr2

pi_fix_chr2 = ggplot(pi_chr2_fix, aes(x = submodel, y = (mean - 2000), group = submodel, color = factor(submodel), fill = factor(submodel))) +
  geom_boxplot() +
  labs(x = "Meiotic frequency (m)",
       y = "Fixation time (generations)",
       title = expression(paste(rho[m]," = 5e-8 (chromosome 2)")),
       color = "Meiotic frequency (m)",
       fill = "Meiotic frequency (m)") +
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")
pi_fix_chr2

# Arrange the plots side by side using plot_grid
combined_plot_chr2 <- plot_grid(pi_var_chr2, pi_fix_chr2 + theme(plot.margin = margin(l = 23)),
                                labels = "AUTO",
                                label_size = 20,
                                ncol = 2,
                                align = "h",  # Horizontal alignment
                                rel_widths = c(1, 1.3)  # Adjust relative widths to account for the legend
)
combined_plot_chr2

ggsave("pi_rhofixvar_chr2.png", plot = combined_plot_chr2, width = 16, height = 6, units = "in")

rm(list=ls())
