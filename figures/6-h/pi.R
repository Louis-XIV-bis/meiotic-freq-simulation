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

csv_to_tibble = function(GR_list, h, windows){
  for (GR in GR_list){ 
    data =  read_csv(paste0("pi_",GR,"_",h,".csv"),col_names = F)
    
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
  
  pi_table = bind_rows(EGR, `100GR`) %>%
    as_tibble() %>% 
    mutate(submodel = case_when(
      submodel == "EGR" ~ "1",
      submodel == "100GR" ~ "0.01",
      TRUE ~ submodel  # Default: no change for other values
    ))
}

########################################################################
# Create tibble from csv data 

windows = seq(1,1000000,length.out=99)
GR_list = c("EGR", "100GR")

pi_table_h0p2 = csv_to_tibble(GR_list, "h0p2", windows)
pi_table_h0p2

pi_table_h0p4 = csv_to_tibble(GR_list, "h0p4", windows)
pi_table_h0p4

pi_table_h0p5 = csv_to_tibble(GR_list, "h0p5", windows)
pi_table_h0p5

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(pi_table_h0p2$ymin, pi_table_h0p4$ymin, pi_table_h0p5$ymin)),
                   max(c(pi_table_h0p2$ymax, pi_table_h0p4$ymax, pi_table_h0p5$ymax)))

# Create color palette 
colors <- c("#1B9E77","#E7298A")

# Using the extracted color palette for both geom_line and geom_ribbon
pi_h0p2 = ggplot(pi_table_h0p2, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi)) + 
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
pi_h0p2

pi_h0p4 = ggplot(pi_table_h0p4, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi ~ " (branch length)"),
       color = "m (meiotic frequency)",
       fill = "m (meiotic frequency)",
       linetype = "m (meiotic frequency)") + 
  theme_light() + 
  ggtitle("h = 0.4") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "top"
  ) + 
  scale_y_continuous(limits = y_axis_limits)
pi_h0p4

pi_h0p5 = ggplot(pi_table_h0p5, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)") +
  theme_light() + 
  ggtitle("h = 0.5") +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
pi_h0p5

# Extract the legend from the main plot
legend_only <- get_legend(pi_h0p4)

# Create the grid of plots using plot_grid
combined_plot <- plot_grid(
  pi_h0p2 + theme(plot.margin = margin(l = 23)),
  pi_h0p4 + theme(legend.position = "none", plot.margin = margin(l = 23)),
  pi_h0p5 + theme(plot.margin = margin(l = 23)), 
  NULL, legend_only, NULL,
  rel_heights = c(3,0.5),
  rel_widths = c(1.2,1,1),
  ncol = 3, nrow = 2, labels = c("A","B","C"), label_size = 20
)
combined_plot

png('pi_h.png',width=16,height=8,units="in",bg = "white", res=600)
combined_plot
dev.off()
