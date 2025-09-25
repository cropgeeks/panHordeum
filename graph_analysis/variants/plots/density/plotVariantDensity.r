args<-commandArgs(TRUE)

#read the BED file with the variant coordinates
variants <- read.table(args[1],sep="\t",header=F)
variantType <- args[2]

#set the column names from the BED file
colnames(variants)<-c("chr","start","end")

#order the chromosomes
chromosomeOrder <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")
#apply this as a factor
variants$chr <- factor(variants$chr,levels=chromosomeOrder)

#read the file with the centromere coordinates
coords <- read.table("MorexV3_centromereCoords.txt", header=T, sep="\t")
coords$chr <- factor(coords$CHROM,levels=chromosomeOrder)

#plot the density as an array of individual plots, one per chromosome, and stack them
library(ggplot2)
variantDensityPlot<-
ggplot(variants) + 
geom_histogram(aes(x=start),binwidth=1e6) + 
geom_vline(data = coords, mapping = aes(xintercept = POS), colour = "red", linetype = "dashed") +
facet_wrap(~ chr,ncol=1) + 
ggtitle(paste(variantType, "density, 1 Mbp bins")) + 
xlab("Position") + 
ylab("# variants") + 
theme_bw(base_size = 40)

#export to a PNG file
png(args[3],width=2000,height=2185)
print(variantDensityPlot)
dev.off()

#also export the counts in the bins to a file
binCounts <- ggplot_build(variantDensityPlot)$data[[1]]
write.table(binCounts,file=paste(variantType,"_binCounts.txt", sep=""),row.names = T,col.names = T,quote = F,sep='\t')

