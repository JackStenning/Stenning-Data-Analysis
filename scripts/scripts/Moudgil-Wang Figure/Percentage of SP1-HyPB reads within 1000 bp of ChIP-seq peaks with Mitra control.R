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
MSP1_MHyPB_1kbp_narrow <- read.table("SP1-HyPB_win_SP1ChIP/Mitra/1kb_Windowoutput_SP1-HyPB_CCpeaks_p9_SP1_ChIP_fq_Mitra.bed", header = TRUE, sep = "\t", col.names = col_names)
MHyPB_1kbp_narrow <- read.table("HyPB_win_SP1ChIP/Mitra/1kb_Windowoutput_HyPB_CCpeaks_p9_SP1_ChIP_fq_Mitra.bed", header = TRUE, sep = "\t", col.names = col_names)
MGenoTTAA <- read.table("TTAA_win_SP1ChIP/Mitra/1kb_Windowoutput_GenomicTTAA__SP1_ChIP_fq_Mitra.bed", header = TRUE, sep = "\t", col.names = Genocol_names)

#Once in the dataframe, we can count how many reads are in each file by the number of lines in each dataframe
MSP1_MHyPBTotalPeaks <- nrow(MSP1_MHyPB_1kbp_narrow )
MHyPBTotalPeaks <- nrow(MHyPB_1kbp_narrow)
MTTAAtotal <- nrow(MGenoTTAA)
#Then we count the number of reads with no overlap - this is denoted by a '0' in the overlap column

MSP1_HyPB_1kbp_a_peaks_b_overlap <- as.numeric(length(which(MSP1_MHyPB_1kbp_narrow$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
MPBase_1kbp_a_peaks_b_overlap <- as.numeric(length(which(MHyPB_1kbp_narrow$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
MTTAA_1kbp_a_peaks_b_overlap <- as.numeric(length(which(MGenoTTAA$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Then we can work out percentage of overlapped reads from the number of non-overlapped reads
MPercent_HyPB_1kbp <- (1-(MSP1_HyPB_1kbp_a_peaks_b_overlap/MSP1_MHyPBTotalPeaks))*100
MPercent_PBase_1kbp <- (1-(MPBase_1kbp_a_peaks_b_overlap/MHyPBTotalPeaks))*100
MPercent_TTAA_1kbp <- (1-(MTTAA_1kbp_a_peaks_b_overlap/MTTAAtotal))*100

#and we can create a table to plot, Col1=Total reads,Col2= non-overlapped reads
Order=c('SP1-HyPB 1kbp', 'Undirected HyPB 1kbp', 'Genomic TTAA 1kbp')
MFinal_Percentage_table <- data.frame(Name=Order,
                                     Total=rep(c(MSP1_MHyPBTotalPeaks, MHyPBTotalPeaks, MTTAAtotal)),
                                     Overlap_withChIP= c(MSP1_HyPB_1kbp_a_peaks_b_overlap, MPBase_1kbp_a_peaks_b_overlap, MTTAA_1kbp_a_peaks_b_overlap),
                                     PercentageOverlapping= c(MPercent_HyPB_1kbp, MPercent_PBase_1kbp, MPercent_TTAA_1kbp))

#Finally, we can now plot
MFigure <- ggplot(data = MFinal_Percentage_table, 
                 aes(x= factor(Name, Order), 
                     y = PercentageOverlapping,
                     fill = Name)) +  geom_col(color = "black",
                                               show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name and Distance from ChIP-Seq Peak", 
                                                                                           y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                           title = "Calling Cards and ChIP-Seq Overlap",
                                                                                           subtitle = "Determining percentage of SP1-HyPB and HyPB peaks within 1000 bp \nof SP1 ChIP-Seq peaks using Mitra control") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed", linewidth = 1) 


MFigure


########################################################################

