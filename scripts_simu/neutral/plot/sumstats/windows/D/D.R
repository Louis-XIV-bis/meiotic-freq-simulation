#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Plot 1 : D --------------------------------------------------------------------------------

rm(list=ls())

windows = seq(1,1000000,length.out=99)

## DATA EGR ######################################################################
data_D_EGR_neutral = read_csv("D_EGR_neutral.csv",col_names = F)

mean = apply(data_D_EGR_neutral, 2, mean)
sd = apply(data_D_EGR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

D_EGR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
D_EGR_neutral

rm(data_D_EGR_neutral)

## DATA 10GR ######################################################################

data_D_10GR_neutral = read_csv("D_10GR_neutral.csv",col_names = F) 

mean = apply(data_D_10GR_neutral, 2, mean)
sd = apply(data_D_10GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

D_10GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="10GR", .before=1)
D_10GR_neutral

rm(data_D_10GR_neutral)

## DATA 100GR ######################################################################
data_D_100GR_neutral = read_csv("D_100GR_neutral.csv",col_names = F) 

mean = apply(data_D_100GR_neutral, 2, mean)
sd = apply(data_D_100GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

D_100GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
D_100GR_neutral

rm(data_D_100GR_neutral)

#####################################################################################

rm(mean,sd,ymin,ymax,windows)

D_neutral = bind_rows(D_EGR_neutral,D_10GR_neutral,D_100GR_neutral) %>%
  as_tibble()
D_neutral

rm(D_EGR_neutral,D_10GR_neutral,D_100GR_neutral)

###### Plot D ######################################################################

D = ggplot(D_neutral, aes(x = windows, y = mean, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel))) +
  labs(title="Mesure du D de Tajima le long du génome", #
       x = "Position (base)",
       y = "D",
       color="Sous-modèles étudiés",
       fill="Sous-modèles étudiés",
       linetype="Sous-modèles étudiés") + 
  theme_light()
ggsave("D.pdf") 
D

rm(list=ls())

windows = seq(1,1000000,length.out=99)

## DATA EGR ######################################################################
data_D_EGR_neutral = read_csv("D_EGR_neutral.csv",col_names = F)

mean = apply(data_D_EGR_neutral, 2, mean)
sd = apply(data_D_EGR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

D_EGR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="EGR", .before=1)
D_EGR_neutral

rm(data_D_EGR_neutral)

## DATA 10GR ######################################################################

data_D_10GR_neutral = read_csv("D_10GR_neutral.csv",col_names = F) 

mean = apply(data_D_10GR_neutral, 2, mean)
sd = apply(data_D_10GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

D_10GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="10GR", .before=1)
D_10GR_neutral

rm(data_D_10GR_neutral)

## DATA 100GR ######################################################################
data_D_100GR_neutral = read_csv("D_100GR_neutral.csv",col_names = F) 

mean = apply(data_D_100GR_neutral, 2, mean)
sd = apply(data_D_100GR_neutral, 2, sd)
ymin = mean-sd
ymax = mean+sd

D_100GR_neutral = data.frame(windows,mean,ymin,ymax) %>%
  add_column(submodel="100GR", .before=1)
D_100GR_neutral

rm(data_D_100GR_neutral)

#####################################################################################

rm(mean,sd,ymin,ymax,windows)

D_neutral = bind_rows(D_EGR_neutral,D_10GR_neutral,D_100GR_neutral) %>%
  as_tibble()
D_neutral

rm(D_EGR_neutral,D_10GR_neutral,D_100GR_neutral)

###### Plot D ######################################################################

D = ggplot(D_neutral, aes(x = windows, y = mean, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel))) +
  labs(title="Mesure du D de Tajima le long du génome", #
       x = "Position (base)",
       y = "D",
       color="Sous-modèles étudiés",
       fill="Sous-modèles étudiés",
       linetype="Sous-modèles étudiés") + 
  theme_light()
ggsave("D.pdf") 
D

rm(list=ls())

