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
setwd("C:/Users/JackP/Documents/GitHub/Stenning_data_analysis/Data/20231018Mitra/") 
#For markdown
#setwd("C:/<Your Directory>/Stenning_data_analysis/Data") 

#Set column names
col_names = c("Chromosome", "Position Start", "Position End", "Number of Overlaps with ChIP-Seq Peaks")
Genocol_names = c("Chromosome", "Position Start", "Position End", 'Name', 'Score', 'Strand',"Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#SP1
SP1_HyPB_ChIP_ENCODE_1000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/1000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_5000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/5000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_10000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/10000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_20000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/20000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_30000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/30000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_40000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/40000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_50000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/50000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_60000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/60000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_70000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/70000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_80000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/80000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_90000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/90000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

SP1_HyPB_ChIP_ENCODE_100000bp = read.table("A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/100000bp_A_GSM4471639_HCT-116_SP1-HyPBase_filtered_fdr_bh_0.05BENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

##Undirected HyPB
HyPB_ChIP_ENCODE_1000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/1000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_5000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/5000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_10000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/10000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_20000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/20000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_30000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/30000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_40000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/40000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_50000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/50000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_60000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/60000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_70000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/70000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_80000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/80000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_90000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/90000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

HyPB_ChIP_ENCODE_100000bp = read.table("A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_SP1_ChIP_ENCODE_controlinput_peaks.narrowPeak/100000bp_A_GSM4471638_HCT-116_HyPBase_filtered_p9_peaks.bed_B_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)

#Genomic TTAA
#Not run yet to save space

#Get total peaks
SP1_HF5_total=nrow(SP1_HyPB_ChIP_ENCODE_1000bp)
HyPB_HF5_total=nrow(HyPB_ChIP_ENCODE_1000bp)
#TTAA_HF5_total=nrow(TTAA_A90)

#Get overlapping peaks
#SP1
SP1_HF5_overlap_1000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_1000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_5000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_5000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_10000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_10000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_20000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_20000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_30000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_30000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_40000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_40000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_50000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_50000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_60000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_60000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_70000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_70000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_80000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_80000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_90000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_90000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
SP1_HF5_overlap_100000bp = as.numeric(length(which(SP1_HyPB_ChIP_ENCODE_100000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#HyPB
HyPB_HF5_overlap_1000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_1000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_5000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_5000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_10000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_10000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_20000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_20000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_30000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_30000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_40000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_40000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_50000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_50000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_60000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_60000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_70000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_70000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_80000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_80000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_90000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_90000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap_100000bp = as.numeric(length(which(HyPB_ChIP_ENCODE_100000bp$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
#SP1
Percentage_SP1_HF5_1000bp = (1-(SP1_HF5_overlap_1000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_5000bp = (1-(SP1_HF5_overlap_5000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_10000bp = (1-(SP1_HF5_overlap_10000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_20000bp = (1-(SP1_HF5_overlap_20000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_30000bp = (1-(SP1_HF5_overlap_30000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_40000bp = (1-(SP1_HF5_overlap_40000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_50000bp = (1-(SP1_HF5_overlap_50000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_60000bp = (1-(SP1_HF5_overlap_60000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_70000bp = (1-(SP1_HF5_overlap_70000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_80000bp = (1-(SP1_HF5_overlap_80000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_90000bp = (1-(SP1_HF5_overlap_90000bp/SP1_HF5_total))*100
Percentage_SP1_HF5_100000bp = (1-(SP1_HF5_overlap_100000bp/SP1_HF5_total))*100

#HyPB
Percentage_HyPB_HF5_1000bp = (1-(HyPB_HF5_overlap_1000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_5000bp = (1-(HyPB_HF5_overlap_5000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_10000bp = (1-(HyPB_HF5_overlap_10000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_20000bp = (1-(HyPB_HF5_overlap_20000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_30000bp = (1-(HyPB_HF5_overlap_30000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_40000bp = (1-(HyPB_HF5_overlap_40000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_50000bp = (1-(HyPB_HF5_overlap_50000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_60000bp = (1-(HyPB_HF5_overlap_60000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_70000bp = (1-(HyPB_HF5_overlap_70000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_80000bp = (1-(HyPB_HF5_overlap_80000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_90000bp = (1-(HyPB_HF5_overlap_90000bp/HyPB_HF5_total))*100
Percentage_HyPB_HF5_100000bp = (1-(HyPB_HF5_overlap_100000bp/HyPB_HF5_total))*100

#Percetage table
graphorder= c('S1', 
              'S5',
              'S10',
              'S20',
              'S30',
              'S40',
              'S50',
              'S60',
              'S70',
              'S80',
              'S90',
              'S100',
              'H1', 
              'H5',
              'H10',
              'H20',
              'H30',
              'H40',
              'H50',
              'H60',
              'H70',
              'H80',
              'H90',
              'H100')
SP1_HyPB_ENCODE_Incrw_Finalpercentage_table = data.frame(Name=graphorder,
                                                 Total=rep(c(SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             SP1_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total,
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total,
                                                             HyPB_HF5_total, 
                                                             HyPB_HF5_total)),
                                                 Overlap_withChIP= c(SP1_HF5_overlap_1000bp,
                                                                     SP1_HF5_overlap_5000bp,
                                                                     SP1_HF5_overlap_10000bp,
                                                                     SP1_HF5_overlap_20000bp,
                                                                     SP1_HF5_overlap_30000bp,
                                                                     SP1_HF5_overlap_40000bp,
                                                                     SP1_HF5_overlap_50000bp,
                                                                     SP1_HF5_overlap_60000bp,
                                                                     SP1_HF5_overlap_70000bp,
                                                                     SP1_HF5_overlap_80000bp,
                                                                     SP1_HF5_overlap_90000bp,
                                                                     SP1_HF5_overlap_100000bp,
                                                                     HyPB_HF5_overlap_1000bp,
                                                                     HyPB_HF5_overlap_5000bp,
                                                                     HyPB_HF5_overlap_10000bp,
                                                                     HyPB_HF5_overlap_20000bp,
                                                                     HyPB_HF5_overlap_30000bp,
                                                                     HyPB_HF5_overlap_40000bp,
                                                                     HyPB_HF5_overlap_50000bp,
                                                                     HyPB_HF5_overlap_60000bp,
                                                                     HyPB_HF5_overlap_70000bp,
                                                                     HyPB_HF5_overlap_80000bp,
                                                                     HyPB_HF5_overlap_90000bp,
                                                                     HyPB_HF5_overlap_100000bp),
                                                 PercentageOverlapping= c(
                                                   Percentage_SP1_HF5_1000bp,
                                                   Percentage_SP1_HF5_5000bp,
                                                   Percentage_SP1_HF5_10000bp,
                                                   Percentage_SP1_HF5_20000bp,
                                                   Percentage_SP1_HF5_30000bp,
                                                   Percentage_SP1_HF5_40000bp,
                                                   Percentage_SP1_HF5_50000bp,
                                                   Percentage_SP1_HF5_60000bp,
                                                   Percentage_SP1_HF5_70000bp,
                                                   Percentage_SP1_HF5_80000bp,
                                                   Percentage_SP1_HF5_90000bp,
                                                   Percentage_SP1_HF5_100000bp,
                                                   Percentage_HyPB_HF5_1000bp,
                                                   Percentage_HyPB_HF5_5000bp,
                                                   Percentage_HyPB_HF5_10000bp,
                                                   Percentage_HyPB_HF5_20000bp,
                                                   Percentage_HyPB_HF5_30000bp,
                                                   Percentage_HyPB_HF5_40000bp,
                                                   Percentage_HyPB_HF5_50000bp,
                                                   Percentage_HyPB_HF5_60000bp,
                                                   Percentage_HyPB_HF5_70000bp,
                                                   Percentage_HyPB_HF5_80000bp,
                                                   Percentage_HyPB_HF5_90000bp,
                                                   Percentage_HyPB_HF5_100000bp))

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
SP1_HyPB_ENCODE_Incrw_Figure <- ggplot(data = SP1_HyPB_ENCODE_Incrw_Finalpercentage_table, 
                               aes(x= factor(Name, graphorder), 
                                   y = PercentageOverlapping,
                                   fill = Name, 
                                   stat("identity"))) +  geom_col(color = "black",
                                                                  show.legend = FALSE) +  ylim(0, 100) + labs(x = "Distance from ChIP-Seq Peak (kbp)", 
                                                                                                              y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                              title = "Calling Cards and ChIP-Seq Overlap",
                                                                                                              subtitle = "SP1-HyPB and HyPB sample produced in my analsysis overlap with SP1 ChIP-Seq \npeaks Using the ENCODE control, increasing the distance window from 1kbp to 100kbp") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed", linewidth = 1) 


SP1_HyPB_ENCODE_Incrw_Figure
########################################################################
