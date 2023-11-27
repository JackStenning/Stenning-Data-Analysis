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
col_names = c("Chromosome", "Position Start", "Position End", "Number of Overlaps with ChIP-Seq Peaks")
Genocol_names = c("Chromosome", "Position Start", "Position End", 'Name', 'Score', 'Strand',"Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_EK_HyPB_A90 = read.table("1kbp_A_ER_HyPB_EKDL230008858_fdr_bh_0.05_peaks_r.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
Undir_EK_HyPB_A90 = read.table("1kbp_A_WT_HyPB_EKDL230008858_p9_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
TTAA_A90 = read.table("1kbp_A_hg38_TTAA.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = Genocol_names)

ESR1_HF5_HyPB_A90 = read.table("1kbp_A_ER_HyPB_HF552DSX7_L2_1_fdr_bh_0.05_peaks_r.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
Undir_HF5_HyPB_A90 = read.table("1kbp_A_WT_HyPB_HF552DSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

ESR1_HF3_HyPB_A90 = read.table("1kbp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_peaks_r.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
Undir_HF3_HyPB_A90 = read.table("1kbp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

#Counting the total reads
ESR1_EK_total=nrow(ESR1_EK_HyPB_A90)
HyPB_EK_total=nrow(Undir_EK_HyPB_A90)
TTAA_total=nrow(TTAA_A90)

ESR1_HF5_total=nrow(ESR1_HF5_HyPB_A90)
HyPB_HF5_total=nrow(Undir_HF5_HyPB_A90)

ESR1_HF3_total=nrow(ESR1_HF3_HyPB_A90)
HyPB_HF3_total=nrow(Undir_HF3_HyPB_A90)
#Count the number of reads with no overlap
ESR1_EK_overlap = as.numeric(length(which(ESR1_EK_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_EK_overlap = as.numeric(length(which(Undir_EK_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
TTAA_overlap = as.numeric(length(which(TTAA_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

ESR1_HF5_overlap = as.numeric(length(which(ESR1_HF5_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF5_overlap = as.numeric(length(which(Undir_HF5_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

ESR1_HF3_overlap = as.numeric(length(which(ESR1_HF3_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_HF3_overlap = as.numeric(length(which(Undir_HF3_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
#Determine overlap percentage
Percentage_ESR1_EK = (1-(ESR1_EK_overlap/ESR1_EK_total))*100
Percentage_HyPB_EK = (1-(HyPB_EK_overlap/HyPB_EK_total))*100
Percentage_TTAA = (1-(TTAA_overlap/TTAA_total))*100

Percentage_ESR1_HF5 = (1-(ESR1_HF5_overlap/ESR1_HF5_total))*100
Percentage_HyPB_HF5 = (1-(HyPB_HF5_overlap/HyPB_HF5_total))*100

Percentage_ESR1_HF3 = (1-(ESR1_HF3_overlap/ESR1_HF3_total))*100
Percentage_HyPB_HF3 = (1-(HyPB_HF3_overlap/HyPB_HF3_total))*100
#Create table to plot
A90_Finalpercentage_table = data.frame(Name=c('ER_EK', 'HyPB_EK', 'TTAA_1', 'ER_HF5', 'HyPB_HF5', 'TTAA_2', 'ER_HF3', 'HyPB_HF3', 'TTAA_3'),
                                       Total=rep(c(ESR1_EK_total, HyPB_EK_total, TTAA_total, ESR1_HF5_total, HyPB_HF5_total, TTAA_total, ESR1_HF3_total, HyPB_HF3_total, TTAA_total)),
                                       Overlap_withChIP= c(ESR1_EK_overlap, HyPB_EK_overlap, TTAA_overlap, ESR1_HF5_overlap, HyPB_HF5_overlap, TTAA_overlap, ESR1_HF3_overlap, HyPB_HF3_overlap, TTAA_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_EK, Percentage_HyPB_EK, Percentage_TTAA, Percentage_ESR1_HF5, Percentage_HyPB_HF5, Percentage_TTAA, Percentage_ESR1_HF3, Percentage_HyPB_HF3, Percentage_TTAA))

#Then we plot
graphorder= c('ER_EK', 'HyPB_EK', 'TTAA_1', 'ER_HF5', 'HyPB_HF5', 'TTAA_2', 'ER_HF3', 'HyPB_HF3', 'TTAA_3')
A90_Figure <- ggplot(data = A90_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "Calling Cards and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from ESR1-HyPB and HyPB \n within 1000 bp of ESR1 ChIP-Seq peaks from sample A90") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


A90_Figure
########################################################################
