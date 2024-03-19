#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');
if (!require('sjPlot', quietly = T)) install.packages('sjPlot');
if (!require('RColorBrewer', quietly = T)) install.packages('RColorBrewer');

library(ggplot2)
library(sjPlot)
library(RColorBrewer)

########################################################################

rm(list=ls())

GR_list = c("EGR", "10GR", "50GR", "100GR")

for (GR in GR_list){ 
  # Create a variable name based on the loop index
  GR_table_name = GR
  
  # Assign a value to the dynamically generated variable
  assign(GR_table_name, read.csv(paste0("fix_",GR,".txt"),header = F))

}

########################################################################

merged_data <- data.frame(
  Value = c(EGR$V1, `10GR`$V1, `50GR`$V1, `100GR`$V1),
  Group = rep(c("1/1", "1/10", "1/50", "1/100"), each = 200)
)

########################################################################

p = ggplot(merged_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(x = "Meiosis to mitosis ratio",
       y = "Generation of fixation",
       color="Meiosis to mitosis ratio",
       fill="Meiosis to mitosis ratio",
       linetype="Meiosis to mitosis ratio") +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )
save_plot("gen_fix.svg", fig = p, width=30, height=20)
ggsave('gen_fix.pdf')
