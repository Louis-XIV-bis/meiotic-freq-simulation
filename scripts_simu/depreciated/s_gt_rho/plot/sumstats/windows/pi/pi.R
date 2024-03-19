#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Plot 1 : pi --------------------------------------------------------------------------------

rm(list=ls())

windows = seq(1,1000000,length.out=499)

## DATA EGR ######################################################################
data_pi_EGR_h0p4 = read_csv("pi_EGR_h0p4.csv",col_names = F)

mean = apply(data_pi_EGR_h0p4, 2, mean)
sd = apply(data_pi_EGR_h0p4, 2, sd)
ymin = mean-sd
ymax = mean+sd

pi_EGR_h0p4 = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
pi_EGR_h0p4

rm(data_pi_EGR_h0p4)

## DATA 100GR ######################################################################
data_pi_100GR_h0p4 = read_csv("pi_100GR_h0p4.csv",col_names = F) 

mean = apply(data_pi_100GR_h0p4, 2, mean)
sd = apply(data_pi_100GR_h0p4, 2, sd)
ymin = mean-sd
ymax = mean+sd

pi_100GR_h0p4 = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
pi_100GR_h0p4

rm(data_pi_100GR_h0p4)

#####################################################################################

rm(mean,sd,ymin,ymax,windows)

pi_h0p4 = bind_rows(pi_EGR_h0p4,pi_100GR_h0p4) %>%
  as_tibble()
pi_h0p4

rm(pi_EGR_h0p4,pi_100GR_h0p4)

###### Plot pi ######################################################################

pi = ggplot(pi_h0p4, aes(x = windows, y = mean, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel))) +
  labs(title="Mesure du pi le long du génome", #
       x = "Position (base)",
       y = expression(pi),
       color="Sous-modèles étudiés",
       fill="Sous-modèles étudiés",
       linetype="Sous-modèles étudiés") + 
  theme_light() +   
  theme(
    plot.title = element_text(size = 0),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )
ggsave("pi.pdf") 
pi

rm(list=ls())

