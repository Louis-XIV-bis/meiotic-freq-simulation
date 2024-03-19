#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

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

rm(list=ls())

windows = seq(1,1000000,length.out=499)

########################################################################

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
    # Add more conditions as needed
    TRUE ~ submodel  # Default: no change for other values
  ))
pi_table

rm(EGR,`100GR`)

########################################################################
num_conditions <- length(unique(pi_table$submodel))
palette_name <- "Dark2"
colors <- brewer.pal(n = num_conditions, name = palette_name)

# Using the extracted color palette for both geom_line and geom_ribbon
pi = ggplot(pi_table, aes(x = windows, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (bp)",
       y = expression(pi ~ "(branch length)"),
       color="Meiotic frequency",
       fill="Meiotic frequency",
       linetype="Meiotic frequency") + 
  theme_light() + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )
pi
ggsave("pi.pdf") 
save_plot("pi.svg", fig = pi, width=30, height=20)


# Same plot without legend 
pi = ggplot(pi_table, aes(x = windows, y = mean, group = submodel)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = submodel),alpha=0.2) +
  geom_line(aes(color = submodel)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "Sequence position (bp)",
       y = expression(pi ~ "(branch length)"))+ 
  theme_light() + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  guides(color = "none", fill = "none")
pi
ggsave("pi_noleg.pdf") 
save_plot("pi_noleg.svg", fig = pi, width=30, height=20)

rm(list=ls())