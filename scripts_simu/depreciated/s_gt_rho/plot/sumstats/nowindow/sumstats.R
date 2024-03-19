#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Data upload and formatting ----------------------------------------

rm(list=ls())

EGR_h0p4nential = read_csv("sumstats_EGR_h0p4.csv",col_names = F) %>%
  add_column(submodel="EGR", .before=1)
EGR_h0p4nential

r100GR_h0p4nential = read_csv("sumstats_100GR_h0p4.csv",col_names = F) %>%
  add_column(submodel="100GR", .before=1)
r100GR_h0p4nential

sumstats = bind_rows(EGR_h0p4nential,r100GR_h0p4nential)
colnames(sumstats) = c("submodel","pi","D","Wtheta","DaF")
sumstats 

rm(EGR_h0p4nential,r100GR_h0p4nential)

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

