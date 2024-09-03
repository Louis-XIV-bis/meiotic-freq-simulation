#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@universite-paris-saclay.fr
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('readr', quietly = T)) install.packages('readr');
if (!require('tibble', quietly = T)) install.packages('tibble');
if (!require('dplyr', quietly = T)) install.packages('dplyr');

library(readr)
library(tibble)
library(dplyr)

########################################################################
# Function that upload and format the data

rm(list=ls())

# Load all full datasets 
data_pi = read_csv('../data/pi_merged.csv') %>%
  mutate(alpha = as.factor(1 / GR))
data_pi

data_pi_chr2 = read_csv('../data/pi_full_chr2_ctrl.csv') %>%
  mutate(alpha = as.factor(1 / GR))
data_pi_chr2

data_t = read_csv('../data/t_merged.csv') %>%
  mutate(alpha = as.factor(1 / GR)) 
data_t

########################################################################
# Statistical tests 

###### Fig 1 ###########################################################

### pi chr 2 alpha = 0.01 vs alpha == 1 ###
chr2_pi_1 <- data_pi_chr2 %>%
  filter(alpha == 1) %>% 
  pull(mean)

hist(chr2_pi_1)
shapiro.test(chr2_pi_1)
# p-value = 0.1855 => gaussian 

chr2_pi_0p01 <- data_pi_chr2 %>%
  filter(alpha == 0.01) %>% 
  pull(mean)

hist(chr2_pi_0p01)
shapiro.test(chr2_pi_0p01)
# p-value = 0.5167 => gaussian 

t.test(chr2_pi_1, chr2_pi_0p01, alternative="greater")
# p-value < 2.2e-16 ==> pi_alpha=1 > pi_alpha=0.01
rm(chr2_pi_1,chr2_pi_0p01)

### T chr 2 alpha = 0.01 vs alpha == 1 ### PB distrib 
chr2_t_1 <- data_t %>%
  filter(alpha == 1 & s == 0.05 & h == 0.5 & rho == '5e-08') %>% 
  pull(mean)

hist(chr2_t_1)
shapiro.test(chr2_t_1)
# p-value = 0.02274 => not gaussian 

chr2_t_0p01 <- data_t %>%
  filter(alpha == 0.01) %>% 
  pull(mean)

hist(chr2_t_0p01)
shapiro.test(chr2_t_0p01)
# p-value = 0.5167 => gaussian 

t.test(chr2_t_1, chr2_t_0p01)
# p-value < 2.2e-16 ==> pi_alpha=1 > pi_alpha=0.01
rm(chr2_pi_1,chr2_pi_0p01)