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

result <- merged_t %>%
  group_by(m) %>%
  summarize(mean_time = mean(time, na.rm = TRUE) -2000)

# Print the result
print(result)
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
shapiro.test(t_100GR$time-2000) 
# p-value = 4.559e-05 => not gaussian 
hist(as.vector(t_100GR$time-2000), probability = TRUE)

shapiro.test(t_100GR$time-2000) 
# p-value = 4.559e-05 => not gaussian 
hist(as.vector(t_100GR$time-2000), probability = TRUE)

shapiro.test(t_100GR$time-2000) 
# p-value = 4.559e-05 => not gaussian 
hist(as.vector(t_100GR$time-2000), probability = TRUE)

kruskal.test(time - 2000 ~ factor(m), data = merged_t)
# p-value = 0.5495 ==> There's NO significiant != between some groups 
dunn.test(merged_t$time - 2000, factor(merged_t$m), method = "bonferroni")
