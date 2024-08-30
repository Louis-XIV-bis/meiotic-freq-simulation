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

low_100GR <- read.table("fix_100GR_0p02.txt", header = FALSE, col.names = "time")
low_100GR$s = "0.02"

mid_100GR <- read.table("fix_100GR_0p1.txt", header = FALSE, col.names = "time")
mid_100GR$s = "0.1"

high_100GR <- read.table("fix_100GR_0p5.txt", header = FALSE, col.names = "time")
high_100GR$s = "0.5"

merged_100GR <- rbind(low_100GR, mid_100GR, high_100GR)

low_EGR <- read.table("fix_EGR_0p02.txt", header = FALSE, col.names = "time")
low_EGR$s = "0.02"

mean(low_EGR$time -2000)

mid_EGR <- read.table("fix_EGR_0p1.txt", header = FALSE, col.names = "time")
mid_EGR$s = "0.1"

high_EGR <- read.table("fix_EGR_0p5.txt", header = FALSE, col.names = "time")
high_EGR$s = "0.5"

merged_EGR <- rbind(low_EGR, mid_EGR, high_EGR)

result <- merged_100GR %>%
  group_by(s) %>%
  summarize(mean_time = mean(time, na.rm = TRUE) -2000)
print(result)

########################################################################

p_EGR = ggplot(merged_EGR, aes(x = s, y = time - 2000, fill = s)) +
  geom_boxplot() +
  labs(x = expression(s),
       y = "Fixation time (generations)",
       color=expression(s),
       fill=expression(s),
       linetype=expression(s)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  ggtitle("m = 0.01") + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16), 
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  guides(color = "none", fill = "none") 
  
p_EGR

p_100GR = ggplot(merged_100GR, aes(x = s, y = time - 2000, fill = s)) +
  geom_boxplot() +
  labs(x = expression(s),
       y = "Fixation time (generations)",
       color=expression(s),
       fill=expression(s),
       linetype=expression(s)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light() + 
  ggtitle("m = 1") + 
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5)
  )
p_100GR

# Arrange the plots side by side using plot_grid
combined_plot <- plot_grid(p_EGR, p_100GR + theme(plot.margin = margin(l = 23)),
                           labels = "AUTO",
                           label_size = 20,
                           ncol = 2,
                           align = "h"  # Horizontal alignment
)
combined_plot

ggsave("s_timefix.png", plot = combined_plot, width = 16, height = 6, units = "in")

# Statistical tests 
shapiro.test(low_100GR$time-2000) 
# p-value = 3.638e-16 => not gaussian 
hist(as.vector(low_100GR$time-2000), probability = TRUE)

shapiro.test(mid_100GR$time-2000) 
# p-value = 1.326e-08 => not gaussian 
hist(as.vector(mid_100GR$time-2000), probability = TRUE)

shapiro.test(high_100GR$time-2000) 
# p-value = 5.069e-09=> not gaussian 
hist(as.vector(high_100GR$time-2000), probability = TRUE)

kruskal.test(time -2000 ~ factor(s), data = merged_100GR)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 
dunn.test(merged_100GR$time - 2000, factor(merged_100GR$s), method = "bonferroni")
# All 100GR groups are different 

# Statistical tests 
shapiro.test(low_EGR$time-2000) 
# p-value = 3.773e-15 => not gaussian 
hist(as.vector(low_EGR$time-2000), probability = TRUE)

shapiro.test(mid_EGR$time-2000) 
# p-value = 0.001881 => not gaussian 
hist(as.vector(mid_EGR$time-2000), probability = TRUE)

shapiro.test(high_EGR$time-2000) 
# p-value = 1.574e-06 => not gaussian 
hist(as.vector(high_EGR$time-2000), probability = TRUE)

kruskal.test(time -2000 ~ factor(s), data = merged_EGR)
# p-value < 2.2e-16 ==> There's a significiant != between some groups 
dunn.test(merged_EGR$time - 2000, factor(merged_EGR$s), method = "bonferroni")
# All EGR groups are different 

# Test pour s = 0.5 et != m 

shapiro.test(mid_EGR$time-2000) # not gaussian
shapiro.test(high_EGR$time-2000)# not gaussian 
mean(mid_EGR$time-2000)
mean(high_EGR$time-2000)

wilcox.test(mid_EGR$time - 2000, high_EGR$time - 2000)
# p-value < 2.2e-16 : for s = 0.5 and values of m (0.01 and 1) are different