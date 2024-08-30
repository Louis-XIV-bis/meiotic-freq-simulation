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
  
  pi_table = bind_rows(EGR,`100GR`) %>%
    as_tibble() %>% 
    mutate(submodel = case_when(
      submodel == "EGR" ~ "1",
      submodel == "100GR" ~ "0.01",
      TRUE ~ submodel  # Default: no change for other values
    ))
}

########################################################################
# Create tibble from csv data 

GR_list = c("EGR", "100GR")
windows = seq(1,2000000,length.out=99)
pi_table_sln = csv_to_tibble(GR_list, "sln", windows) %>% filter(windows > 1000000) # <=> 2nd chromosome
pi_table_sln

########################################################################
# Statistical tests #  

data_EGR = pi_table_sln %>% filter(submodel == 1)
hist(as.vector(data_EGR$mean), probability = TRUE)
shapiro.test(data_EGR$mean) 
# p-value = 0.7131 => Gaussian 

data_100GR = pi_table_sln %>% filter(submodel == 0.01)
hist(as.vector(data_100GR$mean), probability = TRUE)
shapiro.test(data_100GR$mean) 
# p-value = 0.8647 => Gaussian 

# Both gaussian : t.test 

t.test(data_EGR$mean,data_100GR$mean, altenative = "greater")
# t = 35.7, df = 85.301, p-value < 2.2e-16 ==> significatively different
summary(pi_table_sln)

print(aggregate(mean ~ submodel, data = pi_table_sln, FUN = mean))
# submodel         mean
# 1     0.01 2.692781e-05
# 2        1 4.017834e-05

########################################################################
# Plot #  

colors <- c("#1B9E77","#E7298A")

p = ggplot(pi_table_sln, aes(x = submodel, y = mean, fill = submodel)) +
  geom_boxplot() +
  labs(x = "Meiotic frequency",
       y = expression(pi),
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
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  ggtitle("With selection")
p
ggsave('ctrl_pi_sln.png')
