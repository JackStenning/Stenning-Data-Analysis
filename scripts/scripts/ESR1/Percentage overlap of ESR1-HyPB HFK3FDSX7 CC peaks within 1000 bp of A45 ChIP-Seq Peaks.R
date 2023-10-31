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
setwd("C:/Users/JackP/Documents/GitHub/Stenning_data_analysis/Data") 
#For markdown
#setwd("C:/<Your Directory>/Stenning_data_analysis/Data") 

#Set column names
col_names = c("Chromosome", "Position Start", "Position End", "Number of Overlaps with ChIP-Seq Peaks")
Genocol_names = c("Chromosome", "Position Start", "Position End", 'Name', 'Score', 'Strand',"Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_HyPB_A45 = read.table("20231017outputs/1kbp_A_ER_HyPB_HFK3FDSX7_L2_1_fdr_bh_0.05_no-r_peaks.bed_B_GSM2970403_A45_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
Undir_HyPB_A45 = read.table("20231017outputs/1kbp_A_WT_HyPB_HFK3FDSX7_L2_1_p30_peaks.bed_B_GSM2970403_A45_peaks.bed", header = TRUE, sep = "\t", col.names = col_names)
TTAA_A45 = read.table("20231017outputs/1kbp_A_hg38_TTAA.bed_B_GSM2970403_A45_peaks.bed", header = TRUE, sep = "\t", col.names = Genocol_names)

#Counting the total reads
ESR1_total=nrow(ESR1_HyPB_A45)
HyPB_total=nrow(Undir_HyPB_A45)
TTAA_total=nrow(TTAA_A45)

#Count the number of reads with no overlap
ESR1_overlap = as.numeric(length(which(ESR1_HyPB_A45$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_overlap = as.numeric(length(which(Undir_HyPB_A45$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
TTAA_overlap = as.numeric(length(which(TTAA_A45$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1 = (1-(ESR1_overlap/ESR1_total))*100
Percentage_HyPB = (1-(HyPB_overlap/HyPB_total))*100
Percentage_TTAA = (1-(TTAA_overlap/TTAA_total))*100

#Create table to plot
A45_Finalpercentage_table = data.frame(Name=c('ESR1-HyPB 1kbp', 'Undirected HyPB 1kbp', 'Genomic TTAA 1kbp'),
                                       Total=rep(c(ESR1_total, HyPB_total, TTAA_total)),
                                       Overlap_withChIP= c(ESR1_overlap, HyPB_overlap, TTAA_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1, Percentage_HyPB, Percentage_TTAA))

#Then we plot
graphorder= c('ESR1-HyPB 1kbp', 'Undirected HyPB 1kbp', 'Genomic TTAA 1kbp')
A45_Figure <- ggplot(data = A45_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "Calling Cards and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from ESR1-HyPB and HyPB \nsample HFK3FDSX7 within 1000 bp of ESR1 ChIP-Seq peaks from \nsample A45") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed", linewidth = 1) 


A45_Figure
########################################################################
