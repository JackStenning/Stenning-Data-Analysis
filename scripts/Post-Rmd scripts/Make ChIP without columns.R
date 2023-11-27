############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#Load software
library(ggplot2)
library(ggthemes)

#Set working directory
#setwd("//userfs/jps558/w2k/Desktop/Local Calling Cards Analysis/Data/JPS/20231117RmdOutputs") 


#Set column names
col_names = c("Chromosome", "Position Start", "Position End", "Peak name", "Score")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_EK_HyPB_A90 = read.table("GSM2970404_A90_peakss.bed", header = TRUE, sep = "\t", col.names = col_names)
