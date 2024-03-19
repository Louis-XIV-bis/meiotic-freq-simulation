#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Data upload and formatting ----------------------------------------

rm(list=ls()) 

freq = c(1/40,2/40,3/40,4/40,5/40,6/40,7/40,8/40,9/40,10/40,11/40,12/40,13/40,14/40,15/40,16/40,17/40,18/40,19/40,20/40,21/40,22/40,23/40,24/40,25/40,26/40,27/40,28/40,29/40,30/40,31/40,32/40,33/40,34/40,35/40,36/40,37/40,38/40,39/40)

data_EGR_h0p4 = read_csv("normSFS_EGR_h0p4.csv",col_names = F)

mean = apply(data_EGR_h0p4, 2, mean)
sd = apply(data_EGR_h0p4, 2, sd)
ymin = mean-sd
ymax = mean+sd

EGR_h0p4 = data.frame(freq,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
EGR_h0p4

rm(data_EGR_h0p4)

data_100GR_h0p4 = read_csv("normSFS_100GR_h0p4.csv",col_names = F)

mean = apply(data_100GR_h0p4, 2, mean)
sd = apply(data_100GR_h0p4, 2, sd)
ymin = mean-sd
ymax = mean+sd

r100GR_h0p4 = data.frame(freq,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
r100GR_h0p4

rm(data_100GR_h0p4)
rm(mean,sd,ymin,ymax,freq)

norm_SFS_data = bind_rows(EGR_h0p4,r100GR_h0p4) %>%
  as_tibble()
norm_SFS_data

rm(EGR_h0p4,r100GR_h0p4)

norm_SFS_data$submodel <- factor(norm_SFS_data$submodel , levels=c("EGR","100GR"))

###### Plot : normalized SFS for each frequency for each submodel ----------------------------------------

SFSnorm = ggplot(norm_SFS_data, aes(x = freq, y = mean, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel) )) +
  geom_hline(yintercept=1, color="blue",size=0.2) +
  labs(title="Mesure du SFS normalisé",
       x = "Derived allele frequency",
       y = "normalized SFS",
       color="Recombination frequency",
       fill="Recombination frequency",
       linetype="Recombination frequency") + 
  theme_light() +
  scale_fill_manual(values=c("#35C0CA","#F8766D")) + 
  scale_color_manual(values=c("#35C0CA","#F8766D"))


ggsave("SFSnorm.pdf") 
SFSnorm

rm(list=ls())

