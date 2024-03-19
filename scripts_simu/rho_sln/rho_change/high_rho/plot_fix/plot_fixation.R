#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');

library(ggplot2)

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
  Group = rep(c("EGR", "10GR", "50GR", "100GR"), each = 200)
)

########################################################################

ggplot(merged_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(title = "Génération de finxation pour chaque GR",
       x = "GR",
       y = "Génération")
ggsave('gen_fix.pdf')
