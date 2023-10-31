############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages

#install.packages("ggthemes")

#Load software
library(rgl)
library(ggplot2)
library(ggthemes)
##Process data

setwd ("C:/Users/JackP/Documents/GitHub/Stenning_data_analysis/Data/20230925Windowoutputs")

#For markdown
#setwd("C:/<Your Directory>/Stenning_data_analysis/Data") 


# In order to plot the percentage of reads within X base pairs of a CC peak, we needed to run an analysis
#that increased the window size of reads by X bp, and then search for overlaps between two files

#This has already been done on viking with 0.25 and 1 kbp, so we need to turn the outputs into dataframe
col_names = c("Chromosome", "Position Start", "Position End", "Number of insertions", "Number of Overlaps with ChIP-Seq Peaks")
Genocol_names = c("Chromosome", "Position Start", "Position End", 'notsure', 'notsure', 'notsure',"Number of Overlaps with ChIP-Seq Peaks")
SP1_HyPB_M_1kbp_narrow <- read.table("SP1_HyPB/ENCODE/Mitrapeaks/1kb_Windowoutput_MitraSP1-HyPB_CCpeaks_p9_SP1_ChIP_fq_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)
HyPB_1kbp_M_narrow <- read.table("SP1_HyPB/ENCODE/Mitrapeaks/undirected/1kb_Windowoutput_MitraHyPB_CCpeaks_p9_SP1_ChIP_fq_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)
GenoTTAA <- read.table("TTAA_win_SP1ChIP/ENCODE/1kb_Windowoutput_GenomicTTAA__SP1_ChIP_fq_ENCODE.bed", header = TRUE, sep = "\t", col.names = Genocol_names)

#Once in the dataframe, we can count how many reads are in each file by the number of lines in each dataframe
SP1_HyPB_M_TotalPeaks <- nrow(SP1_HyPB_M_1kbp_narrow)
HyPB_M_TotalPeaks <- nrow(HyPB_1kbp_M_narrow)
TTAAtotal <- nrow(GenoTTAA)
#Then we count the number of reads with no overlap - this is denoted by a '0' in the overlap column

SP1_HyPB_M_1kbp_a_peaks_b_overlap <- as.numeric(length(which(SP1_HyPB_M_1kbp_narrow$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
HyPB_M_1kbp_a_peaks_b_overlap <- as.numeric(length(which(HyPB_1kbp_M_narrow$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
TTAA_1kbp_a_peaks_b_overlap <- as.numeric(length(which(GenoTTAA$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Then we can work out percentage of overlapped reads from the number of non-overlapped reads
Percent_SP1_HyPB_1kbp_M <- (1-(SP1_HyPB_M_1kbp_a_peaks_b_overlap/SP1_HyPB_M_TotalPeaks))*100
Percent_HyPB_1kbp_M <- (1-(HyPB_M_1kbp_a_peaks_b_overlap/HyPB_M_TotalPeaks))*100
Percent_TTAA_1kbp <- (1-(TTAA_1kbp_a_peaks_b_overlap/TTAAtotal))*100

#and we can create a table to plot, Col1=Total reads,Col2= non-overlapped reads
Order=c('SP1-HyPB 1kbp', 'Undirected HyPB 1kbp', 'Genomic TTAA 1kbp')
Final_Percentage_table <- data.frame(Name=Order,
                                     Total=rep(c(SP1_HyPB_M_TotalPeaks, HyPB_M_TotalPeaks, TTAAtotal)),
                                     Overlap_withChIP= c(SP1_HyPB_M_1kbp_a_peaks_b_overlap, HyPB_M_1kbp_a_peaks_b_overlap, TTAA_1kbp_a_peaks_b_overlap),
                                     PercentageOverlapping= c(Percent_SP1_HyPB_1kbp_M, Percent_HyPB_1kbp_M, Percent_TTAA_1kbp))

#Finally, we can now plot
Figure <- ggplot(data = Final_Percentage_table, 
                 aes(x= factor(Name, Order), 
                     y = PercentageOverlapping,
                     fill = Name)) +  geom_col(color = "black",
                                               show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak", 
                                                                                           y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                           title = "Calling Cards and ChIP-Seq Overlap",
                                                                                           subtitle = "Determining percentage of SP1-HyPB and HyPB peaks from Mitra \npeak files within 1000 bp of SP1 ChIP-Seq peaks using ENCODE control") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed", linewidth = 1) 


Figure


########################################################################

