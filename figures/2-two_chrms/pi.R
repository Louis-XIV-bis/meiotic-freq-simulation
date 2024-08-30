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

rm(list=ls())

windows = seq(1,2000000,length.out=99)

########################################################################
# Upload, process and store csv data to tibble

GR_list = c("EGR", "100GR")

for (GR in GR_list){ 
  data =  read_csv(paste0("pi_",GR,".csv"),col_names = F)
  
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
pi_table

rm(EGR,`100GR`)

########################################################################
# Create plot #

# Create color palette 
colors <- c("#1B9E77","#E7298A")

# Using the extracted color palette for both geom_line and geom_ribbon
pi = ggplot(pi_table, aes(x = windows/1000000, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (Mbp)",
       y = expression(pi),
       color="Meiotic frequency",
       fill="Meiotic frequency",
       linetype="Meiotic frequency") + 
  theme_light() + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "black") + 
  ylim(-1e-5,8e-5)
pi 

# Adding labels 
pi = ggdraw() +
  draw_plot(pi, width = 1, height = 1) +
  draw_text("Chromosome 1", x = 0.27, y = 0.92, size = 18, color = "black") +
  draw_text("Chromosome 2", x = 0.61, y = 0.92, size = 18, color = "black")
pi

png('pi_twochrms.png',width=14,height=6,units="in",bg = "white", res=600)
pi
dev.off()

# Statistical tests on the second chromosome to see if the delta pi is significant 
# supposed that it's constant for the whole chromosome

chr2_data_EGR = pi_table %>%  # m = 1
  filter(submodel == 1 & windows > 1000000)
chr2_pi_EGR <- as.numeric(chr2_data_EGR$mean)
hist(chr2_pi_EGR, probability = TRUE)
shapiro.test(chr2_pi_EGR)
# p-value = 0.7131 => gaussian 

chr2_data_100GR = pi_table %>%  # m = 0.01
  filter(submodel == 0.01 & windows > 1000000)
chr2_pi_100GR <- as.numeric(chr2_data_100GR$mean)
hist(chr2_pi_100GR, probability = TRUE)
shapiro.test(chr2_pi_100GR) 
# p-value = 0.8647 => gaussian 

t.test(chr2_pi_EGR, chr2_pi_100GR, alternative="greater")
# p-value < 2.2e-16
