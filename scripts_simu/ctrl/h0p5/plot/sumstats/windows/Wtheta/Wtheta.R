#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Plot 1 : Wtheta --------------------------------------------------------------------------------

rm(list=ls())

windows = seq(1,1000000,length.out=499)

## DATA EGR ######################################################################
data_Wtheta_EGR_h0p5 = read_csv("Wtheta_EGR_h0p5.csv",col_names = F)

mean = apply(data_Wtheta_EGR_h0p5, 2, mean)
sd = apply(data_Wtheta_EGR_h0p5, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_EGR_h0p5 = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
Wtheta_EGR_h0p5

rm(data_Wtheta_EGR_h0p5)

## DATA 100GR ######################################################################
data_Wtheta_100GR_h0p5 = read_csv("Wtheta_100GR_h0p5.csv",col_names = F) 

mean = apply(data_Wtheta_100GR_h0p5, 2, mean)
sd = apply(data_Wtheta_100GR_h0p5, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_100GR_h0p5 = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
Wtheta_100GR_h0p5

rm(data_Wtheta_100GR_h0p5)

#####################################################################################

rm(mean,sd,ymin,ymax,windows)

Wtheta_h0p5 = bind_rows(Wtheta_EGR_h0p5,Wtheta_100GR_h0p5) %>%
  as_tibble()
Wtheta_h0p5

rm(Wtheta_EGR_h0p5,Wtheta_100GR_h0p5)

###### Plot Wtheta ######################################################################

Wtheta = ggplot(Wtheta_h0p5, aes(x = windows, y = mean, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel))) +
  labs(title="Mesure du Wtheta le long du génome", #
       x = "Position (base)",
       y = "Wtheta",
       color="Sous-modèles étudiés",
       fill="Sous-modèles étudiés",
       linetype="Sous-modèles étudiés") + 
  theme_light()
ggsave("Wtheta.pdf") 
Wtheta

rm(list=ls())
