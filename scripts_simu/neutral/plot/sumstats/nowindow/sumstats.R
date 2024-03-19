#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Data upload and formatting ----------------------------------------

rm(list=ls())

EGR_exponential = read_csv("sumstats_EGR_expo.csv",col_names = F) %>%
  add_column(submodel="EGR", .before=1)
EGR_exponential

r1000GR_exponential = read_csv("sumstats_1000GR_expo.csv",col_names = F) %>%
  add_column(submodel="1000GR", .before=1)
r1000GR_exponential

sumstats = bind_rows(EGR_exponential,r1000GR_exponential)
colnames(sumstats) = c("submodel","pi","D","Wtheta","DaF")
sumstats 

rm(EGR_exponential,r1000GR_exponential)

###### Plot 1 : pi ----------------------------------------

pi = sumstats %>%
  ggplot(aes(x=factor(submodel), y=pi, fill=submodel)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "Valeurs de pi pour les différents sous-modèles", #
       x = "Sous-modèles",
       y = "pi") +  
  guides(fill=guide_legend(title="Frequency of recombination events")) + 
  theme_light()
ggsave("pi.pdf") #
pi

###### Plot 2 : D ----------------------------------------

D = sumstats %>%
  ggplot(aes(x=factor(submodel), y=D, fill=submodel)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "Valeurs de pi pour les différents sous-modèles", #
       x = "Sous-modèles",
       y = "D") +  
  guides(fill=guide_legend(title="Frequency of recombination events")) + 
  theme_light()
ggsave("D.pdf") #
D

###### Plot 3 : Wtheta ----------------------------------------

Wtheta = sumstats %>%
  ggplot(aes(x=factor(submodel), y=Wtheta, fill=submodel)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "Valeurs de Wtheta pour les différents sous-modèles", #
       x = "Sous-modèles",
       y = "Wtheta") +  
  guides(fill=guide_legend(title="Frequency of recombination events")) + 
  theme_light()
ggsave("Wtheta.pdf") #
Wtheta

###### Plot 4 : DaF ----------------------------------------

DaF = sumstats %>%
  ggplot(aes(x=factor(submodel), y=DaF, fill=submodel)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "Valeurs de DaF pour les différents sous-modèles", #
       x = "Sous-modèles",
       y = "DaF") +  
  guides(fill=guide_legend(title="Frequency of recombination events")) + 
  theme_light()
ggsave("DaF.pdf") #
DaF

library(cowplot)
x = plot_grid(pi,D,Wtheta,DaF)
x
save_plot("multiplot.pdf",x)

rm(list=ls())
