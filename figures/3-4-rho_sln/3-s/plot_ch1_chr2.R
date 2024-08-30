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
windows = seq(1,2000000,length.out=99)
GR_list = c("EGR", "10GR", "50GR", "100GR")

pi_table_ctrl = csv_to_tibble(GR_list, "ctrl", windows)
pi_ctrl_chr1 = pi_table_ctrl %>% filter(windows < 1000000)

pi_table_fav = csv_to_tibble(GR_list, "high_s", windows)
pi_fav_chr1 = pi_table_fav %>% filter(windows < 1000000)

pi_table_unfav = csv_to_tibble(GR_list, "low_s", windows)
pi_unfav_chr1 = pi_table_unfav %>% filter(windows < 1000000)

pi_table_ctrl = csv_to_tibble(GR_list, "ctrl", windows)
pi_ctrl_chr2 = pi_table_ctrl %>% filter(windows > 1000000)

pi_table_fav = csv_to_tibble(GR_list, "high_s", windows)
pi_fav_chr2 = pi_table_fav %>% filter(windows > 1000000)

pi_table_unfav = csv_to_tibble(GR_list, "low_s", windows)
pi_unfav_chr2 = pi_table_unfav %>% filter(windows > 1000000)

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(pi_table_ctrl$ymin, pi_table_fav$ymin, pi_table_unfav$ymin)),
                   max(c(pi_table_ctrl$ymax, pi_table_fav$ymax, pi_table_unfav$ymax)))

# Create color palette 
num_conditions <- length(unique(pi_table_ctrl$submodel))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
plot_unfav_chr1 = ggplot(pi_unfav_chr1, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi)) + 
  theme_light() + 
  ggtitle("s = 0.02") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_unfav_chr1

plot_ctrl_chr1 = ggplot(pi_ctrl_chr1, aes(x = windows/1000000, y = mean, group = submodel)) +
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
  ggtitle("s = 0.1") + 
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
plot_ctrl_chr1

plot_fav_chr1 = ggplot(pi_fav_chr1, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)") +
  theme_light() + 
  ggtitle("s = 0.5") +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_fav_chr1

# Using the extracted color palette for both geom_line and geom_ribbon
plot_unfav_chr2 = ggplot(pi_unfav_chr2, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "m (meiotic frequency)",
       y = expression(pi)) + 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18),
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_unfav_chr2

plot_ctrl_chr2 = ggplot(pi_ctrl_chr2, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "m (meiotic frequency)",
       y = expression(pi ~ " (branch length)"),
       color = "m (meiotic frequency)",
       fill = "m (meiotic frequency)",
       linetype = "m (meiotic frequency)") + 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
  ) +   
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_ctrl_chr2

plot_fav_chr2 = ggplot(pi_fav_chr2, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "m (meiotic frequency)") +
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(), 
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_fav_chr2

# Extract the legend from the main plot
legend_only <- get_legend(plot_ctrl_chr1)

# Create the grid of plots using plot_grid
combined_plot <- plot_grid(
  plot_unfav_chr1 + theme(plot.margin = margin(l = 23)),
  plot_ctrl_chr1 + theme(legend.position = "none", plot.margin = margin(l = 23)),
  plot_fav_chr1 + theme(plot.margin = margin(l = 23)), 
  plot_unfav_chr2 + theme(plot.margin = margin(l = 23, t = 23)),
  plot_ctrl_chr2 + theme(plot.margin = margin(l = 23, t = 23)),
  plot_fav_chr2 + theme(plot.margin = margin(l = 23, t = 23)), 
  NULL, legend_only, NULL,
  rel_heights = c(3,3,0.5),
  rel_widths = c(1.2,1,1),
  ncol = 3, nrow = 3, labels = c("A","B","C","D","E","F"), label_size = 20
)
combined_plot
png('pi_sln_s_wg.png',width=16,height=12,units="in",bg = "white", res=300)
combined_plot
dev.off()

