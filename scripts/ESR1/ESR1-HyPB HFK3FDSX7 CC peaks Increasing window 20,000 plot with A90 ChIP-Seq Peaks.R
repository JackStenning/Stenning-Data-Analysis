############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#Load software
library(ggplot2)
library(ggthemes)

#Set working directory

## SEARCH "C:/<Your Directory>" and replace it with the correct path

#To test
setwd("C:/Users/JackP/Documents/GitHub/Stenning_data_analysis/Data/20231017outputs/increasingwindow/") 
#For markdown
#setwd("C:/<Your Directory>/Stenning_data_analysis/Data") 

#Set column names
col_names = c("Chromosome", "Position Start", "Position End", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
##ESR1
ESR1_HyPB_HF5_A90_1000bp = read.table("A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed/1000bp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

ESR1_HyPB_HF5_A90_5000bp = read.table("A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed/5000bp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

ESR1_HyPB_HF5_A90_10000bp = read.table("A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed/10000bp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

ESR1_HyPB_HF5_A90_15000bp = read.table("A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed/15000bp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

ESR1_HyPB_HF5_A90_20000bp = read.table("A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed/20000bp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

##Undirected HyPB
HyPB_HF5_A90_1000bp = read.table("A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed/1000bp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_HF5_A90_5000bp = read.table("A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed/5000bp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_HF5_A90_10000bp = read.table("A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed/10000bp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_HF5_A90_15000bp = read.table("A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed/15000bp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_HF5_A90_20000bp = read.table("A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed/20000bp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

#Genomic TTAA
#Not run yet to save space

#Get total peaks
ESR1_HF5_total=nrow(ESR1_HyPB_HF5_A90_1000bp)
HyPB_HF5_total=nrow(HyPB_HF5_A90_1000bp)
#TTAA_HF5_total=nrow(TTAA_A90)

#Get overlapping peaks
#ESR1
ESR1_HF5_overlap_1000bp = as.numeric(length(which(ESR1_HyPB_HF5_A90_1000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_HF5_overlap_5000bp = as.numeric(length(which(ESR1_HyPB_HF5_A90_5000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_HF5_overlap_10000bp = as.numeric(length(which(ESR1_HyPB_HF5_A90_10000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_HF5_overlap_15000bp = as.numeric(length(which(ESR1_HyPB_HF5_A90_15000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_HF5_overlap_20000bp = as.numeric(length(which(ESR1_HyPB_HF5_A90_20000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#HyPB
HyPB_HF5_overlap_1000bp = as.numeric(length(which(HyPB_HF5_A90_1000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_5000bp = as.numeric(length(which(HyPB_HF5_A90_5000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_10000bp = as.numeric(length(which(HyPB_HF5_A90_10000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_15000bp = as.numeric(length(which(HyPB_HF5_A90_15000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_20000bp = as.numeric(length(which(HyPB_HF5_A90_20000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
#ESR1
Percentage_ESR1_HF5_1000bp = (1-(ESR1_HF5_overlap_1000bp/ESR1_HF5_total))*100
Percentage_ESR1_HF5_5000bp = (1-(ESR1_HF5_overlap_5000bp/ESR1_HF5_total))*100
Percentage_ESR1_HF5_10000bp = (1-(ESR1_HF5_overlap_10000bp/ESR1_HF5_total))*100
Percentage_ESR1_HF5_15000bp = (1-(ESR1_HF5_overlap_15000bp/ESR1_HF5_total))*100
Percentage_ESR1_HF5_20000bp = (1-(ESR1_HF5_overlap_20000bp/ESR1_HF5_total))*100
#HyPB
Percentage_HyPB_HF5_1000bp = (1-(HyPB_HF5_overlap_1000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_5000bp = (1-(HyPB_HF5_overlap_5000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_10000bp = (1-(HyPB_HF5_overlap_10000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_15000bp = (1-(HyPB_HF5_overlap_15000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_20000bp = (1-(HyPB_HF5_overlap_20000bp/HyPB_HF5_total))*100

#Percetage table
graphorder= c('ER - 1', 
              'ER - 5',
              'ER - 10',
              'ER - 15',
              'ER - 20', 
              'HyPB 1', 
              'HyPB 5',
              'HyPB 10',
              'HyPB 15',
              'HyPB 20')
HF5_A90_Incrw_Finalpercentage_table = data.frame(Name=graphorder,
                                       Total=rep(c(ESR1_total, 
                                                   ESR1_total, 
                                                   ESR1_total, 
                                                   ESR1_total, 
                                                   ESR1_total, 
                                                   HyPB_total, 
                                                   HyPB_total, 
                                                   HyPB_total, 
                                                   HyPB_total, 
                                                   HyPB_total)),
                                       Overlap_withChIP= c(ESR1_HF5_overlap_1000bp,
                                                           ESR1_HF5_overlap_5000bp,
                                                           ESR1_HF5_overlap_10000bp,
                                                           ESR1_HF5_overlap_15000bp,
                                                           ESR1_HF5_overlap_20000bp,
                                                           HyPB_HF5_overlap_1000bp,
                                                           HyPB_HF5_overlap_5000bp,
                                                           HyPB_HF5_overlap_10000bp,
                                                           HyPB_HF5_overlap_15000bp,
                                                           HyPB_HF5_overlap_20000bp),
                                       PercentageOverlapping= c(
                                         Percentage_ESR1_HF5_1000bp,
                                         Percentage_ESR1_HF5_5000bp,
                                         Percentage_ESR1_HF5_10000bp,
                                         Percentage_ESR1_HF5_15000bp,
                                         Percentage_ESR1_HF5_20000bp,
                                         Percentage_HyPB_HF5_1000bp,
                                         Percentage_HyPB_HF5_5000bp,
                                         Percentage_HyPB_HF5_10000bp,
                                         Percentage_HyPB_HF5_15000bp,
                                         Percentage_HyPB_HF5_20000bp))

#WhenTTAA needed
#'Genomic TTAA 1kbp',
#'Genomic TTAA 5kbp',
#'Genomic TTAA 10kbp',
#'Genomic TTAA 15kbp',
#'Genomic TTAA 20kbp'
#TTAA_total,
#TTAA_total,
#TTAA_total,
#TTAA_total,
#TTAA_total,

#Plot
HF5_A90_Incrw_Figure <- ggplot(data = HF5_A90_Incrw_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak (kbp)", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "Calling Cards and ChIP-Seq Overlap",
                                                                                                    subtitle = "ESR1-HyPB and HyPB sample HFK3FDSX7 overlap with ESR1 ChIP-Seq \npeaks from sample A90, increasing the distance window from 1kbp to 20kbp") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed", linewidth = 1) 


HF5_A90_Incrw_Figure
########################################################################
