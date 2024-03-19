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
windows = seq(1,2000000,length.out=499)
GR_list = c("EGR", "10GR", "50GR", "100GR")

pi_table_ctrl = csv_to_tibble(GR_list, "ctrl", windows)
pi_ctrl_chr2 = pi_table_ctrl %>% filter(windows > 1000000)
pi_ctrl_chr2$rho = 5e-8

pi_table_low = csv_to_tibble(GR_list, "fav", windows)
pi_low_chr2 = pi_table_low %>% filter(windows > 1000000)
pi_low_chr2$rho = 1e-8 

pi_table_high = csv_to_tibble(GR_list, "unfav", windows)
pi_high_chr2 = pi_table_high %>% filter(windows > 1000000)
pi_high_chr2$rho = 1e-7

merged_df = rbind(pi_ctrl_chr2, pi_low_chr2, pi_high_chr2)

m_1 = merged_df %>% filter(submodel == 1)
m_0p1 = merged_df %>% filter(submodel == 0.1)
m_0p05 = merged_df %>% filter(submodel == 0.05)
m_0p01 = merged_df %>% filter(submodel == 0.01)

########################################################################

# Statistical tests 
# m = 1 # 
low_m1 = m_1 %>% filter(rho == 1e-8)
shapiro.test(low_m1$mean) 
# p-value = 0.000386 => non gaussian distribution

mid_m1 = m_1 %>% filter(rho == 5e-8) 
shapiro.test(mid_m1$mean)
# p-value = 0.04208 => non gaussian distribution

high_m1 = m_1 %>% filter(rho == 1e-7) 
shapiro.test(high_m1$mean)
# p-value = 0.7007 => gaussian distribution

pi_m1= rbind(low_m1, mid_m1, high_m1)

kruskal.test(mean ~ factor(rho), data = pi_m1)
# p-value < 0.0825 ==> No significiant != between some groups 

# m = 0.1 # 
low_m0p1 = m_0p1 %>% filter(rho == 1e-8)
shapiro.test(low_m0p1$mean) 
# p-value = 0.1319 => gaussian distribution

mid_m0p1 = m_0p1 %>% filter(rho == 5e-8) 
shapiro.test(mid_m0p1$mean)
# p-value = 0.6104 => gaussian distribution

high_m0p1 = m_0p1 %>% filter(rho == 1e-7) 
shapiro.test(high_m0p1$mean)
# p-value = 0.6104 => gaussian distribution

pi_m0p1 = rbind(low_m0p1, mid_m0p1, high_m0p1)

kruskal.test(mean ~ factor(rho), data = pi_m0p1) # run krustal to have the same tests everywhere
# p-value < 0.8979 ==> No significiant != between some groups

# m = 0.05 # 
low_m0p05 = m_0p05 %>% filter(rho == 1e-8)
shapiro.test(low_m0p05$mean) 
# p-value = 0.9431 => gaussian distribution

mid_m0p05 = m_0p05 %>% filter(rho == 5e-8) 
shapiro.test(mid_m0p05$mean)
# p-value = 0.01173 => non gaussian distribution

high_m0p05 = m_0p05 %>% filter(rho == 1e-7) 
shapiro.test(high_m0p05$mean)
# p-value = 4.669e-05 => non gaussian distribution

pi_m0p05 = rbind(low_m0p05, mid_m0p05, high_m0p05)

kruskal.test(mean ~ factor(rho), data = pi_m0p05) # run krustal to have the same tests everywhere
# p-value < 3.56e-05 ==> No significiant != between some groups
dunn.test(pi_m0p05$mean, factor(pi_m0p05$rho), method = "bonferroni")

# m = 0.01 # 
low_m0p01 = m_0p01 %>% filter(rho == 1e-8)
shapiro.test(low_m0p01$mean) 
# p-value = 0.6631 => gaussian distribution

mid_m0p01 = m_0p01 %>% filter(rho == 5e-8) 
shapiro.test(mid_m0p01$mean)
# p-value = 0.6383 => gaussian distribution

high_m0p01 = m_0p01 %>% filter(rho == 1e-7) 
shapiro.test(high_m0p01$mean)
# p-value = 0.272 => gaussian distribution

pi_m0p01 = rbind(low_m0p01, mid_m0p01, high_m0p01)

kruskal.test(mean ~ factor(rho), data = pi_m0p01) # run krustal to have the same tests everywhere
# p-value < 2.2e-16 ==> No significiant != between some groups
dunn.test(pi_m0p01$mean, factor(pi_m0p01$rho), method = "bonferroni")

########################################################################
# Create plots individually #

# Explicitly set the y-axis limits to be the same for both plots
y_axis_limits <- c(min(c(pi_ctrl_chr2$ymin, pi_low_chr2$ymin, pi_high_chr2$ymin)),
                   max(c(pi_ctrl_chr2$ymax, pi_low_chr2$ymax, pi_high_chr2$ymax)))

# Create color palette 
num_conditions <- length(unique(pi_ctrl_chr2$submodel))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

plot_m1 = ggplot(m_1, aes(x = factor(rho), y = mean, fill = factor(rho))) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = expression(rho),
       y = expression(pi ~ "(branch length)")) + 
  theme_light() + 
  ggtitle("m = 1") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_m1

plot_m0p1 = ggplot(m_0p1, aes(x = factor(rho), y = mean, fill = factor(rho))) +
  geom_boxplot() +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = expression(rho),
       y = expression(pi ~ "(branch length)")) + 
  theme_light() + 
  ggtitle("m = 0.1") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
  ) + 
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_m0p1

plot_m0p05 = ggplot(m_0p05, aes(x = factor(rho), y = mean, fill = factor(rho))) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = expression(rho),
       y = expression(pi ~ "(branch length)")) + 
  theme_light() + 
  ggtitle("m = 0.05") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_m0p05

plot_m0p01 = ggplot(m_0p01, aes(x = factor(rho), y = mean, fill = factor(rho))) +
  geom_boxplot() +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = expression(rho),
       y = expression(pi ~ "(branch length)")) + 
  theme_light() + 
  ggtitle("m = 0.01") + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
  ) + 
  guides(color = "none", fill = "none") + 
  scale_y_continuous(limits = y_axis_limits)
plot_m0p01

# Create the grid of plots using plot_grid
combined_plot <- plot_grid(
  plot_m1 + theme(plot.margin = margin(l = 23)),
  plot_m0p1 + theme(plot.margin = margin(l = 23)),
  plot_m0p05 + theme(plot.margin = margin(l = 23)), 
  plot_m0p01 + theme(plot.margin = margin(l = 23)), 
  ncol = 2, nrow = 2, labels = c("A","B","C", "D"), label_size = 20,
  rel_widths = c(1.2,1,1))
combined_plot
png('rhovar_chr2_m.png',width=16,height=8,units="in",bg = "white", res=600)
combined_plot
dev.off()
