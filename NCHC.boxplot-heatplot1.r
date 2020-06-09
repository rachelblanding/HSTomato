if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
Fulldata<-import("C:/Users/Rache/Downloads/counts.txt")
head(Fulldata)
ShortData<-Fulldata[,-(1:6)]
NCHCData<-ShortData[,c(1,3,5,13,15,17)]
NCData<-ShortData[,c(7,9,11)]
HCData<-ShortData[,c(1,3,5)]
myCPM<-cpm(NCData)
#Create barplot for DNA sizes
?DGEList
dgeSD<-DGEList(NCHCData)
dgeSD$samples
dgeSD$samples$lib.size
barplot(dgeSD$samples$lib.size, names=colnames(dgeSD),col="red")
title("Barplot of N-C and H-C DNA sizes")

logSD<-cpm(dgeSD, log = TRUE)

#Create logged boxplot of DNA sizes
boxplot(logSD, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logSD),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#Create MDS Plot
plotMDS(dgeSD, main="N-C and H-C MDS plot")

#Make graphs side by side
par(mfrow=c(1,2))

# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(NCData)

## Choose color for each cell
col.cell <- c("purple","orange", "darkgreen","grey","red","blue")[NCHCData$`H1-C.bam`]
?data.frame
data.frame(Fulldata$Geneid, col.cell)
infodata<-Fulldata[,1:6]

# Redo the MDS with cell type colouring
plotMDS(dgeSD,col=col.cell)
title("N-C and H-C DNA MDS Plot")

# Let's add a legend to the plot so we know which colours correspond to which cell type
?legend
legend(legend = Fulldata$Geneid,"topleft",fill=c("purple","orange", "darkgreen", "grey") )

# Add a title
title("M-C and H-C DNA MDS Plot")

# Similarly for status
levels(NCHCData$`M1-C.bam`)
col.status <- c("blue","red")[Fulldata$`N1-C.bam`]
col.status
plotMDS(dgeSD,col=col.status)
legend("bottomright",fill=c("blue","red"),legend=Fulldata$Strand,cex=0.8)
title("N-C and H-C Strand")

# Dimension 3 appears to separate pregnant samples from the rest. Dim4?
plotMDS(dgeSD,dim=c(3,4))
title("3 and 4D MDS plot")
Fulldata2<-Fulldata[1:24, ]
labelx <- paste(Fulldata2$Geneid, Fulldata2$Strand, Fulldata2$Length)
groupx <- paste(Fulldata2$Length,Fulldata2$Strand, sep = ".")
groupx <- factor(groupx)

#Create MDS file
glMDSPlot(dgeSD, labels=labelx, groups=groupx, folder="mds")

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logSD, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logSD[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[ShortData$`H1-C.bam`]
col.cell

# Plot the heatmap
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Heatmap of Tomato DNA", scale="row")

# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes\nacross samples",ColSideColors=col.cell,scale="row")
dev.off()
??ColSideColors

# Apply normalisation to DGEList object
dgeObj <- calcNormFactors(dgeSD)
dgeObj$samples
par(mfrow=c(1,2))

#Compare MD plots for each strand
plotMD(logSD,column = 1)
abline(h=0,col="grey")
plotMD(logSD,column = 4)
abline(h=0,col="lightblue")
plotMD(logSD,column = 2)
abline(h=0,col="grey")
plotMD(logSD,column = 5)
abline(h=0,col="lightblue")
plotMD(logSD,column = 3)
abline(h=0,col="grey")
plotMD(logSD,column = 6)
abline(h=0,col="lightblue")
par(mfrow=c(1,2))
plotMD(dgeObj,column = 1)
abline(h="grey")
plotMD(dgeObj,column = 4)
abline(h=0,col="grey")
plotMD(dgeObj,column = 2)
abline(h="grey")
plotMD(dgeObj,column = 5)
abline(h=0,col="grey")
plotMD(dgeObj,column = 3)
abline(h="grey")
plotMD(dgeObj,column = 6)
abline(h=0,col="grey")
