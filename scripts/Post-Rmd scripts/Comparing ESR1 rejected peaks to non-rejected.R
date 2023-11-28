############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#Load software
library(ggplot2)
library(ggthemes)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/Calling Cards analysis/Stenning_analysis_followup/Rejectedpeaks/results") 


#Set column names
col_names = c("Chromosome", "Position Start", "Position End", "Number of Overlaps with ChIP-Seq Peaks")
Genocol_names = c("Chromosome", "Position Start", "Position End", 'Name', 'Score', 'Strand',"Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_EK_HyPB_A90 = read.table("1kbp_A_ER_HyPB_EKDL230008858_rejected_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
Undir_EK_HyPB_A90 = read.table("1kbp_A_WT_HyPB_EKDL230008858_rejected_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

ESR1_accept_HyPB_A90 = read.table("accepted/1kbp_A_ER_HyPB_EKDL230008858_fdr_bh_0.05_peaks_r.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
Undir_Accept_HyPB_A90 = read.table("accepted/1kbp_A_WT_HyPB_EKDL230008858_p30_peaks.bed_B_GSM2970404_A90_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)

#Counting the total reads
ESR1_A90_total=nrow(ESR1_EK_HyPB_A90)
HyPB_A90_total=nrow(Undir_EK_HyPB_A90)


ESR1_accept_total=nrow(ESR1_accept_HyPB_A90)
HyPB_accept_total=nrow(Undir_Accept_HyPB_A90)


#Count the number of reads with no overlap
ESR1_A90_overlap = as.numeric(length(which(ESR1_EK_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_A90_overlap = as.numeric(length(which(Undir_EK_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

ESR1_accept_overlap = as.numeric(length(which(ESR1_accept_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_accept_overlap = as.numeric(length(which(Undir_Accept_HyPB_A90$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_A90 = (1-(ESR1_A90_overlap/ESR1_A90_total))*100
Percentage_HyPB_A90 = (1-(HyPB_A90_overlap/HyPB_A90_total))*100

Percentage_ESR1_accept = (1-(ESR1_accept_overlap/ESR1_accept_total))*100
Percentage_HyPB_accept = (1-(HyPB_accept_overlap/HyPB_accept_total))*100


#Create table to plot
RejectVsAccept_Finalpercentage_table = data.frame(Name=c('ER_Reject', 'ER_Accept', 'HyPB_Reject','HyPB_Accept'),
                                       Total=rep(c(ESR1_A90_total, ESR1_accept_total, HyPB_A90_total, HyPB_accept_total)),
                                       Number_Not_Overlapping_withChIP= c(ESR1_A90_overlap, ESR1_accept_overlap, HyPB_A90_overlap, HyPB_accept_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_A90, Percentage_ESR1_accept, Percentage_HyPB_A90, Percentage_HyPB_accept))

#Then we plot
graphorder= c('ER_Reject', 'ER_Accept', 'HyPB_Reject','HyPB_Accept')
A90_Figure <- ggplot(data = RejectVsAccept_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "Calling Cards and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlap from ESR1-HyPB and HyPB \n accepted and rejected peaks within 1000 bp of ESR1 ChIP-Seq peaks from sample A90") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


A90_Figure
write.csv(RejectVsAccept_Finalpercentage_table, "//userfs/jps558/w2k/Desktop/Calling Cards analysis/Stenning_analysis_followup/Rejectedpeaks/results/Comparing ESR1 rejected peaks to non-rejected.csv", row.names=FALSE)
########################################################################
