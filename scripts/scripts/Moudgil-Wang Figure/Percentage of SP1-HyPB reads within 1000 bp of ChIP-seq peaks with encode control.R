############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages

#install.packages("ggthemes")

#Load software
library(rgl)
library(ggplot2)
library(ggthemes)
##Process data
###work directory
setwd ("C:/Users/JackP/Documents/GitHub/Stenning_data_analysis/Data/20230912Windowoutputs")

#For markdown
#setwd("C:/<Your Directory>/Stenning_data_analysis/Data") 

# In order to plot the percentage of reads within X base pairs of a CC peak, we needed to run an analysis
#that increased the window size of reads by X bp, and then search for overlaps between two files

#This has already been done on viking with 0.25 and 1 kbp, so we need to turn the outputs into dataframe
col_names = c("Chromosome", "Position Start", "Position End", "Number of Overlaps with ChIP-Seq Peaks")
Genocol_names = c("Chromosome", "Position Start", "Position End", 'notsure', 'notsure', 'notsure',"Number of Overlaps with ChIP-Seq Peaks")
SP1_HyPB_1kbp_narrow <- read.table("SP1-HyPB_win_SP1ChIP/ENCODE/1kb_Windowoutput_SP1-HyPB_CCpeaks_p9_SP1_ChIP_fq_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)
HyPB_1kbp_narrow <- read.table("HyPB_win_SP1ChIP/ENCODE/1kb_Windowoutput_HyPB_CCpeaks_p9_SP1_ChIP_fq_ENCODE.bed", header = TRUE, sep = "\t", col.names = col_names)
GenoTTAA <- read.table("TTAA_win_SP1ChIP/ENCODE/1kb_Windowoutput_GenomicTTAA__SP1_ChIP_fq_ENCODE.bed", header = TRUE, sep = "\t", col.names = Genocol_names)

#Once in the dataframe, we can count how many reads are in each file by the number of lines in each dataframe
SP1_HyPBTotalPeaks <- nrow(SP1_HyPB_1kbp_narrow)
HyPBTotalPeaks <- nrow(HyPB_1kbp_narrow)
TTAAtotal <- nrow(GenoTTAA)
#Then we count the number of reads with no overlap - this is denoted by a '0' in the overlap column

SP1_HyPB_1kbp_a_peaks_b_overlap <- as.numeric(length(which(SP1_HyPB_1kbp_narrow$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
PBase_1kbp_a_peaks_b_overlap <- as.numeric(length(which(HyPB_1kbp_narrow$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
TTAA_1kbp_a_peaks_b_overlap <- as.numeric(length(which(GenoTTAA$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Then we can work out percentage of overlapped reads from the number of non-overlapped reads
Percent_HyPB_1kbp <- (1-(SP1_HyPB_1kbp_a_peaks_b_overlap/SP1_HyPBTotalPeaks))*100
Percent_PBase_1kbp <- (1-(PBase_1kbp_a_peaks_b_overlap/HyPBTotalPeaks))*100
Percent_TTAA_1kbp <- (1-(TTAA_1kbp_a_peaks_b_overlap/TTAAtotal))*100

#and we can create a table to plot, Col1=Total reads,Col2= non-overlapped reads
Order=c('SP1-HyPB 1kbp', 'Undirected HyPB 1kbp', 'Genomic TTAA 1kbp')
Final_Percentage_table <- data.frame(Name=Order,
                                     Total=rep(c(SP1_HyPBTotalPeaks, HyPBTotalPeaks, TTAAtotal)),
                                     Overlap_withChIP= c(SP1_HyPB_1kbp_a_peaks_b_overlap, PBase_1kbp_a_peaks_b_overlap, TTAA_1kbp_a_peaks_b_overlap),
                                     PercentageOverlapping= c(Percent_HyPB_1kbp, Percent_PBase_1kbp, Percent_TTAA_1kbp))

#Finally, we can now plot
Figure <- ggplot(data = Final_Percentage_table, 
                 aes(x= factor(Name, Order), 
                     y = PercentageOverlapping,
                     fill = Name)) +  geom_col(color = "black",
                                               show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak", 
                                                                                           y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                           title = "Calling Cards and ChIP-Seq Overlap",
                                                                                           subtitle = "Determining percentage of SP1-HyPB and HyPB peaks within 1000 bp \nof SP1 ChIP-Seq peaks using ENCODE control") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed", linewidth = 1) 


Figure


########################################################################

