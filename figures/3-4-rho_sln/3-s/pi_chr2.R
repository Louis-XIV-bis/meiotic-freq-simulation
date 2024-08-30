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
if (!require('dunn.test', quietly = T)) install.packages('cowplot');

library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(sjPlot)
library(RColorBrewer)
library(cowplot)
library(dunn.test)

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
pi_ctrl_chr2 = pi_table_ctrl %>% filter(windows > 1000000)

pi_table_fav = csv_to_tibble(GR_list, "high_s", windows)
pi_fav_chr2 = pi_table_fav %>% filter(windows > 1000000)

pi_table_unfav = csv_to_tibble(GR_list, "low_s", windows)
pi_unfav_chr2 = pi_table_unfav %>% filter(windows > 1000000)

########################################################################
# Statistical tests 
# S = 0.02 # 
unfav_10GR = pi_unfav_chr2 %>% filter(submodel == 0.1) 
shapiro.test(unfav_10GR$mean) 
# p-value = 0.8033 => gaussian distribution

unfav_50GR = pi_unfav_chr2 %>% filter(submodel == 0.05)
shapiro.test(unfav_50GR$mean) 
# p-value = 0.5196 => gaussian distribution

unfav_100GR = pi_unfav_chr2 %>% filter(submodel == 0.01) 
shapiro.test(unfav_100GR$mean) 
# p-value = 0.548 => gaussian distribution

unfav_EGR = pi_unfav_chr2 %>% filter(submodel == 1) 
shapiro.test(unfav_EGR$mean) 
# p-value = 0.2225 => gaussian distribution

pi_table_chr2_unfav = rbind(unfav_EGR, unfav_100GR, unfav_50GR, unfav_10GR)
kruskal.test(mean ~ factor(submodel), data = pi_table_chr2_unfav)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 

dunn.test(pi_table_chr2_unfav$mean, factor(pi_table_chr2_unfav$submodel), method = "bonferroni")

mean(unfav_100GR$mean)
mean(unfav_10GR$mean)
mean(unfav_50GR$mean)
mean(unfav_EGR$mean)

# S = 0.1 # 
ctrl_10GR = pi_ctrl_chr2 %>% filter(submodel == 0.1) 
shapiro.test(ctrl_10GR$mean) 
# p-value = 0.2155 => gaussian distribution

ctrl_50GR = pi_ctrl_chr2 %>% filter(submodel == 0.05)
shapiro.test(ctrl_50GR$mean) 
# p-value = 0.2173 => gaussian distribution

ctrl_100GR = pi_ctrl_chr2 %>% filter(submodel == 0.01) 
shapiro.test(ctrl_100GR$mean) 
# p-value = 0.0376 => NOT a gaussian distribution

ctrl_EGR = pi_ctrl_chr2 %>% filter(submodel == 1) 
shapiro.test(ctrl_EGR$mean) 
# p-value = 0.03252 => NOT a gaussian distribution

pi_table_chr2_ctrl = rbind(ctrl_EGR, ctrl_100GR, ctrl_50GR, ctrl_10GR)
kruskal.test(mean ~ factor(submodel), data = pi_table_chr2_ctrl)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 

dunn.test(pi_table_chr2_ctrl$mean, factor(pi_table_chr2_ctrl$submodel), method = "bonferroni")
mean(ctrl_50GR$mean)
mean(ctrl_10GR$mean)

# S = 0.5 # 
fav_10GR = pi_fav_chr2 %>% filter(submodel == 0.1) 
shapiro.test(fav_10GR$mean) 
# p-value = 0.8033 => gaussian distribution

fav_50GR = pi_fav_chr2 %>% filter(submodel == 0.05)
shapiro.test(fav_50GR$mean) 
# p-value = 0.5196 => gaussian distribution

fav_100GR = pi_fav_chr2 %>% filter(submodel == 0.01) 
shapiro.test(fav_100GR$mean) 
# p-value = 0.548 => gaussian distribution

fav_EGR = pi_fav_chr2 %>% filter(submodel == 1) 
shapiro.test(fav_EGR$mean) 
# p-value = 0.2225 => gaussian distribution

pi_table_chr2_fav = rbind(fav_EGR, fav_100GR, fav_50GR, fav_10GR)
kruskal.test(mean ~ factor(submodel), data = pi_table_chr2_fav)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 

dunn.test(pi_table_chr2_fav$mean, factor(pi_table_chr2_fav$submodel), method = "bonferroni")

fav_EGR$s = 0.5
unfav_EGR$s = 0.02
ctrl_EGR$s = 0.1

shapiro.test(fav_EGR$mean) # gaussian
shapiro.test(unfav_EGR$mean) # gaussian 
shapiro.test(ctrl_EGR$mean) # not gaussian

data_m1 = rbind(fav_EGR, unfav_EGR, ctrl_EGR)
data_m1

ggplot(data_m1, aes(x = factor(s), y = mean, fill = factor(s), group = factor(s))) +
  geom_boxplot() 

kruskal.test(mean ~ factor(s), data = data_m1)
dunn.test(data_m1$mean, factor(data_m1$s), method = "bonferroni")

print(mean(fav_EGR$mean))
print(mean(unfav_EGR$mean))
print(mean(ctrl_EGR$mean))

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(pi_ctrl_chr2$ymin, pi_fav_chr2$ymin, pi_unfav_chr2$ymin)),
                   max(c(pi_ctrl_chr2$ymax, pi_fav_chr2$ymax, pi_unfav_chr2$ymax)))

# Create color palette 
num_conditions <- length(unique(pi_ctrl_chr2$submodel))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
plot_unfav_chr2 = ggplot(pi_unfav_chr2, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Meiotic frequency (m)",
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
plot_unfav_chr2

print(min(pi_fav_chr2$mean))

plot_ctrl_chr2 = ggplot(pi_ctrl_chr2, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Meiotic frequency (m)",
       y = expression(pi ~ " (branch length)"),
       color = "Meiotic frequency",
       fill = "Meiotic frequency",
       linetype = "Meiotic frequency") + 
  theme_light() + 
  ggtitle("s = 0.1 (control)") + 
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
plot_ctrl_chr2

plot_fav_chr2 = ggplot(pi_fav_chr2, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Meiotic frequency (m)") +
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
plot_fav_chr2

# Extract the legend from the main plot
legend_only <- get_legend(plot_ctrl_chr2)

# Create the grid of plots using plot_grid
combined_plot <- plot_grid(
  plot_unfav_chr2 + theme(plot.margin = margin(l = 23)),
  plot_ctrl_chr2 + theme(legend.position = "none", plot.margin = margin(l = 23)),
  plot_fav_chr2 + theme(plot.margin = margin(l = 23)), 
  NULL, legend_only, NULL,
  rel_heights = c(3,0.5),
  rel_widths = c(1.2,1,1),
  ncol = 3, nrow = 2, labels = c("A","B","C"), label_size = 20
)
combined_plot
png('pi_sln_s_chr2.png',width=16,height=8,units="in",bg = "white", res=600)
combined_plot
dev.off()
