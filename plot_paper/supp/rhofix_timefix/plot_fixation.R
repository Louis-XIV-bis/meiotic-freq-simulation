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
if (!require('dunn.test', quietly = T)) install.packages('cowplot');

library(ggplot2)
library(sjPlot)
library(RColorBrewer)
library(dplyr)
library(readr)
library(cowplot)
library(dunn.test)

########################################################################
rm(list=ls())

# Done only for neutral (not fav / unfav)
t_100GR <- read.table("fix_100GR.txt", header = FALSE, col.names = "time")
t_100GR$m = "0.01"

t_50GR <- read.table("fix_50GR.txt", header = FALSE, col.names = "time")
t_50GR$m = "0.05"

t_10GR <- read.table("fix_10GR.txt", header = FALSE, col.names = "time")
t_10GR$m = "0.1"

t_EGR <- read.table("fix_EGR.txt", header = FALSE, col.names = "time")
t_EGR$m = "1"

merged_t <- rbind(t_100GR, t_50GR, t_10GR, t_EGR)

########################################################################

t_rhofix = ggplot(merged_t, aes(x = m, y = time - 2000, fill = m)) +
  geom_boxplot() +
  labs(x = expression(s),
       y = "Fixation time (generations)",
       color=expression(m),
       fill=expression(m),
       linetype=expression(m)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  ggtitle(expression(paste(rho[m]," = 5e-8 (chromosome 2)"))) + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size=16), 
    plot.title = element_text(size = 20, hjust = 0.5)
  )
t_rhofix

x = merged_t %>% filter(m == 0.01)
mean(x$time -2000)

ggsave("rhofix_timefix.png", plot = t_rhofix)


# Statistical tests 
shapiro.test(low_100GR$time-2000) 
# p-value = 2.2e-16 => not gaussian 
hist(as.vector(low_100GR$time-2000), probability = TRUE)

shapiro.test(mid_100GR$time-2000) 
# p-value = 6.645e-08 => not gaussian 
hist(as.vector(mid_100GR$time-2000), probability = TRUE)

shapiro.test(high_100GR$time-2000) 
# p-value = 1.645e-07 => not gaussian 
hist(as.vector(high_100GR$time-2000), probability = TRUE)

kruskal.test(time -2000 ~ factor(s), data = merged_100GR)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 
dunn.test(merged_100GR$time - 2000, factor(merged_100GR$s), method = "bonferroni")

# Statistical tests 
shapiro.test(low_EGR$time-2000) 
# p-value = 1.393e-14 => not gaussian 
hist(as.vector(low_EGR$time-2000), probability = TRUE)

shapiro.test(mid_EGR$time-2000) 
# p-value = 1.876e-05 => not gaussian 
hist(as.vector(mid_EGR$time-2000), probability = TRUE)

shapiro.test(high_EGR$time-2000) 
# p-value = 9.185e-06 => not gaussian 
hist(as.vector(high_EGR$time-2000), probability = TRUE)

kruskal.test(time -2000 ~ factor(s), data = merged_EGR)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 
dunn.test(merged_EGR$time - 2000, factor(merged_EGR$s), method = "bonferroni")

# Test pour s = 0.5 et != m 

rm(list=ls())
shapiro.test(mid_EGR$time-2000) # not gaussian
shapiro.test(high_EGR$time-2000)# not gaussian 
mean(mid_EGR$time-2000)
mean(high_EGR$time-2000)

wilcox.test(mid_EGR$time - 2000, high_EGR$time - 2000)
