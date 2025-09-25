message("load libs")
library(ggplot2)
library(zoo)
library(dplyr)

message("read the input files")
file_paths <- list.files(path = "../", pattern = "*.pansel.txt", full.names = TRUE)
# Function to read and add filename as a column
read_and_label <- function(file) 
{
  data <- read.delim(file, header = FALSE)
  data$file_name <- basename(file)
  return(data)
}

# Read and combine all files into one dataframe
df <- do.call(rbind, lapply(file_paths, read_and_label))

colnames(df) <- c("bin_ID","bin_targetedStartPos","bin_targetedEndPos","Jaccard_index","numDiffPathsBetwAnchors","totalNumPathsBetwAnchors","startPosUsed","endPosUsed","leftMostAnchorNode_ID","leftMostAnchorNode_start","leftMostAnchorNode_end","rightMostAnchorNode_ID","rightMostAnchorNode_start","rightMostAnchorNode_end", "chrom")

#change the file names to chromosome names
df$chrom <- gsub(".pansel.txt", "", df$chrom)

#order the chromosomes
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")
#apply this as a factor
df$chrom <- factor(df$chrom,levels=chromosomeOrder)

#read the file with the centromere coordinates
centromereCoords <- read.table("MorexV3_centromereCoords.txt", header=T, sep="\t")
centromereCoords$chrom <- factor(centromereCoords$CHROM,levels=chromosomeOrder)

#drop rows with nan - the nans are all in the Jaccard column so we can't plot these values
df2 <- na.omit(df)

message("calculate a rolling mean")
df2$rollingMeans <- rollmean(df2$Jaccard_index, k=1000, fill=NA, align='left')

#high density area plot (in lieu of barplot where large numbers of values need to be plotted)
#axis label number is calculated automatically
message("make the plot")
divPlot <-
ggplot(df2,aes(x=startPosUsed)) + 
geom_line(aes(y=rollingMeans), color="blue") + 
geom_vline(data = centromereCoords, mapping = aes(xintercept = POS), colour = "red", linetype = "dashed", size = 2) +
facet_wrap(~ chrom,ncol=1) + 
ggtitle("Jaccard index") + 
xlab("Position") + 
ylab("Jaccard index") + 
theme(axis.text.x = element_text(angle = 90), text = element_text(size = 25)) + 
scale_x_continuous(breaks = pretty(df2$startPosUsed, n = 20))

message("export plot to a PNG file")
png("jaccard.png", width=1000,height=1500)
print(divPlot)
dev.off()

message("done")


