#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Plot 1 : Wtheta --------------------------------------------------------------------------------

rm(list=ls())

windows = seq(1,1000000,length.out=99)

## DATA EGR ######################################################################
data_Wtheta_EGR_neutral = read_csv("Wtheta_EGR_neutral.csv",col_names = F)

mean = apply(data_Wtheta_EGR_neutral, 2, mean)
sd = apply(data_Wtheta_EGR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_EGR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
Wtheta_EGR_neutral

rm(data_Wtheta_EGR_neutral)

## dATA 10GR ######################################################################

data_Wtheta_10GR_neutral = read_csv("Wtheta_10GR_neutral.csv",col_names = F) 

mean = apply(data_Wtheta_10GR_neutral, 2, mean)
sd = apply(data_Wtheta_10GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_10GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="10GR", .before=1)
Wtheta_10GR_neutral

rm(data_Wtheta_10GR_neutral)

## dATA 100GR ######################################################################
data_Wtheta_100GR_neutral = read_csv("Wtheta_100GR_neutral.csv",col_names = F) 

mean = apply(data_Wtheta_100GR_neutral, 2, mean)
sd = apply(data_Wtheta_100GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_100GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
Wtheta_100GR_neutral

rm(data_Wtheta_100GR_neutral)

#####################################################################################

rm(mean,sd,ymin,ymax,windows)

Wtheta_neutral = bind_rows(Wtheta_EGR_neutral,Wtheta_10GR_neutral,Wtheta_100GR_neutral) %>%
  as_tibble()
Wtheta_neutral

rm(Wtheta_EGR_neutral,Wtheta_10GR_neutral,Wtheta_100GR_neutral)

###### Plot Wtheta ######################################################################

Wtheta = ggplot(Wtheta_neutral, aes(x = windows, y = mean, group=submodel, fill=submodel)) +
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

windows = seq(1,1000000,length.out=99)

## WthetaATA EGR ######################################################################
data_Wtheta_EGR_neutral = read_csv("Wtheta_EGR_neutral.csv",col_names = F)

mean = apply(data_Wtheta_EGR_neutral, 2, mean)
sd = apply(data_Wtheta_EGR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_EGR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
Wtheta_EGR_neutral

rm(data_Wtheta_EGR_neutral)

## WthetaATA 10GR ######################################################################

data_Wtheta_10GR_neutral = read_csv("Wtheta_10GR_neutral.csv",col_names = F) 

mean = apply(data_Wtheta_10GR_neutral, 2, mean)
sd = apply(data_Wtheta_10GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_10GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="10GR", .before=1)
Wtheta_10GR_neutral

rm(data_Wtheta_10GR_neutral)

## WthetaATA 100GR ######################################################################
data_Wtheta_100GR_neutral = read_csv("Wtheta_100GR_neutral.csv",col_names = F) 

mean = apply(data_Wtheta_100GR_neutral, 2, mean)
sd = apply(data_Wtheta_100GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

Wtheta_100GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
Wtheta_100GR_neutral

rm(data_Wtheta_100GR_neutral)

#####################################################################################

rm(mean,sd,ymin,ymax,windows)

Wtheta_neutral = bind_rows(Wtheta_EGR_neutral,Wtheta_10GR_neutral,Wtheta_100GR_neutral) %>%
  as_tibble()
Wtheta_neutral

rm(Wtheta_EGR_neutral,Wtheta_10GR_neutral,Wtheta_100GR_neutral)

###### Plot Wtheta ######################################################################

Wtheta = ggplot(Wtheta_neutral, aes(x = windows, y = mean, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel))) +
  labs(title="Mesure du Wtheta de Tajima le long du génome", #
       x = "Position (base)",
       y = "Wtheta",
       color="Sous-modèles étudiés",
       fill="Sous-modèles étudiés",
       linetype="Sous-modèles étudiés") + 
  theme_light()
ggsave("Wtheta.pdf") 
Wtheta

rm(list=ls())

