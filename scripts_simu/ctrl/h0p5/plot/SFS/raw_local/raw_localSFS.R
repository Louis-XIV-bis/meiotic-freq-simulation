#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr 
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('tidyverse', quietly = T)) install.packages('tidyverse');
library(tidyverse)

###### Data upload and formatting ----------------------------------------

rm(list=ls()) 

# Define interval of windows we want (500 windows in my analysis)
start_windows = 200
stop_windows = 300
gap_windows = stop_windows - start_windows

## DATA EGR ## 

data_EGR_h0p5 = read_csv("rawSFS_local_EGR_h0p5.csv",col_names = F)
data_EGR_h0p5 = data_EGR_h0p5[,(1+41*(start_windows-1)):(1+41*(start_windows-1)+41*gap_windows-1)]

# Conversion to absolute values to frequency
for (i in c(1:100)){
  for (windows in c(0:(gap_windows-1))){
    list_col_windows = c((41*windows+1):(41*windows+41))
    data_EGR_h0p5[i,list_col_windows] = data_EGR_h0p5[i,list_col_windows]/sum(data_EGR_h0p5[i,list_col_windows])
  }
}

rm(i, list_col_windows, windows)

# The output give frequency from 0/40 and 40/40 which are not used for SFS, we need to delete them
to_delet_freq = sort(c(seq(41,(41*gap_windows),by=41),seq(1,41*gap_windows-(41-1), by=41)))
to_delet_freq = seq(1,41*gap_windows,by=1)[-c(to_delet_freq)]
to_delet_freq

rawSFS_local = apply(data_EGR_h0p5, 2, mean) %>%
  as_tibble() %>% 
  slice(to_delet_freq) 

rm(to_delet_freq)

# DAF for each windows (mean(SFSi*freq(i))
freq_SFS = seq(1/40,39/40, by=1/40)
freq_SFS

DAF = rep(-1,gap_windows-1) # init
sd = rep(-1,gap_windows-1)

for(i in c(0:(gap_windows-1))){ # gap _ windows -1
  start = i*39+1
  stop = (i+1)*39
  
  DAF_data = rawSFS_local[start:stop,]*freq_SFS
  DAF[i+1] = mean(DAF_data$value)
  sd[i+1] = sd(DAF_data$value)
}

rm(DAF_data, start, stop, rawSFS_local)

ymin = DAF-sd
ymax = DAF+sd
rm(sd)

# From windows number to position (1Mb genome) : a windows for each DAF 
win_size = 2000 # 1 windows = 2kb (1Mb / 500 windows)
position = seq((start_windows * win_size), ((stop_windows - 1) * win_size), by = win_size)
position

SFSraw_local_EGR = data.frame(DAF,ymin,ymax,position) %>% 
  as_tibble() %>% 
  add_column(submodel="EGR", .before=1)
SFSraw_local_EGR

rm(DAF,ymin,ymax,position,i)
rm(data_EGR_h0p5)

## DATA 100GR ## 

data_100GR_h0p5 = read_csv("rawSFS_local_100GR_h0p5.csv",col_names = F)
data_100GR_h0p5 = data_100GR_h0p5[,(1+41*(start_windows-1)):(1+41*(start_windows-1)+41*gap_windows-1)]

# Conversion to absolute values to frequency
for (i in c(1:100)){
  for (windows in c(0:(gap_windows-1))){
    list_col_windows = c((41*windows+1):(41*windows+41))
    data_100GR_h0p5[i,list_col_windows] = data_100GR_h0p5[i,list_col_windows]/sum(data_100GR_h0p5[i,list_col_windows])
  }
}

rm(i, list_col_windows, windows)

# The output give frequency from 0/40 and 40/40 which are not used for SFS, we need to delete them
to_delet_freq = sort(c(seq(41,(41*gap_windows),by=41),seq(1,41*gap_windows-(41-1), by=41)))
to_delet_freq = seq(1,41*gap_windows,by=1)[-c(to_delet_freq)]
to_delet_freq

rawSFS_local = apply(data_100GR_h0p5, 2, mean) %>%
  as_tibble() %>% 
  slice(to_delet_freq) 

rm(to_delet_freq)

# DAF for each windows (mean(SFSi*freq(i))

freq_SFS = seq(1/40,39/40, by=1/40)
freq_SFS

DAF = rep(-1,gap_windows-1) # init
sd = rep(-1,gap_windows-1)

for(i in c(0:(gap_windows-1))){ # gap _ windows -1
  start = i*39+1
  stop = (i+1)*39
  
  DAF_data = rawSFS_local[start:stop,]*freq_SFS
  DAF[i+1] = mean(DAF_data$value)
  sd[i+1] = sd(DAF_data$value)
}

rm(DAF_data, start, stop, rawSFS_local)

ymin = DAF-sd
ymax = DAF+sd
rm(sd)

# From windows number to position (1Mb genome) : a windows for each DAF 
win_size = 2000 # 1 windows = 2kb (1Mb / 500 windows)
position = seq((start_windows * win_size), ((stop_windows - 1) * win_size), by = win_size)
position

SFSraw_local_100GR = data.frame(DAF,ymin,ymax,position) %>% 
  as_tibble() %>% 
  add_column(submodel="100GR", .before=1)
SFSraw_local_100GR

rm(DAF,ymin,ymax,position,i)
rm(data_100GR_h0p5)

### 

DAF_local = bind_rows(SFSraw_local_EGR,SFSraw_local_100GR) 

DAF_local$submodel <- factor(DAF_local$submodel , levels=c("EGR","100GR"))

DAF_local
rm(SFSraw_local_EGR,SFSraw_local_100GR)
rm(start_windows, stop_windows, win_size, gap_windows,freq_SFS)

###### Plot : rawalized SFS for each frequency for each submodel ----------------------------------------

DAF = ggplot(DAF_local, aes(x=position, y=DAF, group=submodel, fill=submodel)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax),alpha=0.2) +
  geom_line(aes(color = factor(submodel))) +
  labs(title="Mesure du DAF le long du génome",
       x = "Derived allele frequency",
       y = "Observed DAF",
       color="Recombination frequency",
       fill="Recombination frequency",
       linetype="Recombination frequency") + 
  theme_light() +
  scale_fill_manual(values=c("#35C0CA","#F8766D")) + 
  scale_color_manual(values=c("#35C0CA","#F8766D"))

ggsave("DAFlocal.pdf") 
DAF

rm(list=ls())

