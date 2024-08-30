#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('ggplot2', quietly = T)) install.packages('ggplot2');
if (!require('sjPlot', quietly = T)) install.packages('sjPlot');
if (!require('RColorBrewer', quietly = T)) install.packages('RColorBrewer');
if (!require('dplyr', quietly = T)) install.packages('dyplr');
if (!require('readr', quietly = T)) install.packages('readr');

library(ggplot2)
library(sjPlot)
library(RColorBrewer)
library(dplyr)
library(readr)

########################################################################
# Function that upload and format the data
rm(list=ls())

csv_to_tibble = function(GR_list){
  for (GR in GR_list){ 
    # Create a variable name based on the loop index
    GR_table_name = GR
    
    # Assign a value to the dynamically generated variable
    assign(GR_table_name, read.csv(paste0("fix_",GR,".txt"),header = F))
    
  }
  
  merged_data <- data.frame(
    Value = c(EGR$V1, `100GR`$V1),
    Group = rep(c("1", "0.01"), each = 200)
  )
}

########################################################################
# Create tibble from csv data 

GR_list = c("EGR", "100GR")
time_table = csv_to_tibble(GR_list)
colnames(time_table) = c("time", "submodel")

########################################################################
# Statistical tests #  

data_EGR = time_table %>% filter(submodel == 1)
hist(as.vector(time_table$time), probability = TRUE)
shapiro.test(data_EGR$time) 
# p-value = 6.494e-08 => not Gaussian 

data_100GR = time_table %>% filter(submodel == 0.01)
hist(as.vector(data_100GR$time), probability = TRUE)
shapiro.test(data_100GR$time) 
# p-value = 2.318e-09 => not Gaussian 

# Both not Gaussian : t.test 
wilcox.test(data_EGR$time,data_100GR$time)
# p-value = 0.6383 ==> pas significativement != 

########################################################################
colors <- c("#1B9E77","#E7298A")

p = ggplot(time_table, aes(x = submodel, y = time - 2000, fill = submodel)) +
  geom_boxplot() +
  labs(x = "Meiotic frequency",
       y = "Fixation time (generations)",
       color="Meiotic frequency",
       fill="Meiotic frequency",
       linetype="Meiosis frequency") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_light() + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )
p
ggsave('ctrl_genfix.png')
