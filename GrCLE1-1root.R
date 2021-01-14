
#Using DESEq2 from the count matrix
setwd("/Users/m.achom/Documents/cornell/WangLab/CLE_project/fromCounts")
library(DESeq2)
library(limma)
countData <- data.matrix(read.delim("counts_allReps.txt", sep="\t", header=T, quote= "", fill=F, dec = ",", row.names = 1))
# for all the lines and control together to apply DESeq2 globally
countData <- countData[, -c(7:9,16:18)] #kept only GrCLE1-1_12roots, GrCLE4-3_10roots, GrCLE1-6_7roots and VDL
colnames(countData) <- c("GrCLE1_12roots1", "GrCLE1_12roots2", "GrCLE1_12roots3", "GrCLE1_22roots1", "GrCLE1_22roots2", "GrCLE1_22roots3", "GrCLE1_6roots1", "GrCLE1_6roots2", "GrCLE1_6roots3", "GrCLE4_3roots1", "GrCLE4_3roots2", "GrCLE4_3roots3", "VDL_4roots1", "VDL_4roots2", "VDL_4roots3")
#preparing metadata
coldata <- data.frame(row.names = colnames(countData), group = rep(c("GrCLE1_12roots", "GrCLE1_22roots", "GrCLE1_6roots", "GrCLE4_3roots", "VDL"), each=3), treatment = c(rep("overexpressed", 12), rep("control", 3)))
coldata$treatment <- factor(x = coldata$treatment,levels = c("overexpressed","control"))
coldata
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ treatment)
dds$treatment <- relevel(dds$treatment, ref = "control") #read note on factor levels on bioconductor vignettes.
# the function relevel also makes that the default log2 fold changes are calculated as treatment over control and not the other way around.
dds
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ] #keeping only the genes with counts > 1
nrow(dds)
vsd <- vst(dds, blind = FALSE) #applying normalisation with VST
head(assay(vsd), 3)
vsd_dat <- assay(vst(dds, blind = FALSE))
write.table(vsd_dat, file="vst_filtered_normalised_counts.txt", sep = "\t", quote = FALSE)

library("dplyr")
library("ggplot2")
plot <- plotPCA(vsd, intgroup = c("treatment", "group")) #PCA for all reps based on normalised VST count matrix
dev.off();
plot

#differential expression analysis and QC check
dds <- DESeq(dds) 
resultsNames(dds)
res <- results(dds, name = "treatment_overexpressed_vs_control")
write.csv(res, file= "DESeq_generated_exp.csv")
summary(res)
head(res)
mcols(res, use.names=TRUE)

#####
coldata <- data.frame(row.names = colnames(countData), group= rep(c("VDL", "GrCLE1_12roots", "GrCLE1_22roots", "GrCLE1_6_5roots", "GrCLE4_10roots"), each=3), condition = c(rep("control",3), rep("overexpressed", 12)))
coldata$condition <- factor(x = coldata$condition,levels = c("control","overexpressed"))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ condition)
dds
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "condition_overexpressed_vs_control")
write.csv(res, file= "DESeq_generated_exp.csv")
summary(res)
plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
nrow(dds)
try <- dds[ rowSums(counts(dds)) > 1, ] #filtering out low counts
nrow(try)


### analysing separately for each of the overexpressed lines as compared to control VDL
#GrCLE1-1
countData <- data.matrix(read.delim("counts_allReps.txt", sep="\t", header=T, quote= "", fill=F, dec = ",", row.names = 1))
countData <- countData[, -c(7:18)] #keeping only two line of GrCLE1-1 to plot PCA
colnames(countData) <- c("GrCLE1_12roots1", "GrCLE1_12roots2", "GrCLE1_12roots3", "GrCLE1_22roots1", "GrCLE1_22roots2", "GrCLE1_22roots3", "VDL_4roots1", "VDL_4roots2", "VDL_4roots3")
coldata <- data.frame(row.names = colnames(countData), group = rep(c("GrCLE1_12roots", "GrCLE1_22roots", "VDL"), each=3), treatment = c(rep("overexpressed", 6), rep("control", 3)))
coldata$treatment <- factor(x = coldata$treatment,levels = c("overexpressed","control"))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ treatment)
dds$treatment <- relevel(dds$treatment, ref = "control")
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ] #keeping only the genes with counts > 1
nrow(dds)
vsd <- vst(dds, blind = FALSE) #applying normalisation with VST
head(assay(vsd), 3)
vsd_dat <- assay(vst(dds, blind = FALSE))
write.table(vsd_dat, file="GrCLE1_12roots_vst_normalised_counts.txt", sep = "\t", quote = FALSE)

plot <- plotPCA(vsd, intgroup = c("treatment", "group")) #PCA for all reps based on normalised VST count matrix
plot

# now I will be further considering the two GrCLE1-1 lines separately
#GrCLE1-1-12_roots
countData <- countData[,-c(4:6)] # keeping only GrCLE1-1_12_roots
colnames(countData) <- c("GrCLE1_12roots1", "GrCLE1_12roots2", "GrCLE1_12roots3", "VDL_4roots1", "VDL_4roots2", "VDL_4roots3")
coldata <- data.frame(row.names = colnames(countData), group = rep(c("GrCLE1_12roots", "VDL"), each=3), treatment = c(rep("overexpressed", 3), rep("control", 3)))
coldata$treatment <- factor(x = coldata$treatment,levels = c("overexpressed","control"))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ treatment)
dds$treatment <- relevel(dds$treatment, ref = "control")
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ] #keeping only the genes with counts > 1
nrow(dds)
vsd <- vst(dds, blind = FALSE) #applying normalisation with VST
head(assay(vsd), 3)
vsd_dat <- assay(vst(dds, blind = FALSE))
write.table(vsd_dat, file="vst_filtered_GrCLE1_1_12roots.txt", sep = "\t", quote = FALSE)


plot <- plotPCA(vsd, intgroup = c("treatment", "group")) #PCA for all reps based on normalised VST count matrix
plot
# starting from here DEG analysis will be done separately on each of the vsd applied normalisation
dds <- DESeq(dds)
res_GrCLE1_12roots <- results(dds, name = "treatment_overexpressed_vs_control")
write.csv(res, file= "DESeq2_res_GrCLE1_12roots.csv")
summary(res) #importnat to note
head(res)
mcols(res, use.names=TRUE)
resSig <- subset(res_GrCLE1_12roots, padj < 0.05)
write.csv(resSig, file = "sig_GrCLE1_12roots.csv")

#To read the file of shared genes between GrCLE1-1_12roots and GrCLE1-1_22roots
setwd("./separate/GrCLE1-1/")

#volcano plot of GrCLE1-12_roots
# better to do this plot at last as a comparison between GrCLE1-1 and GrCLE4-3
aplha <- 0.05
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(), main="Volcano plot: GrCLE1-1_12roots", xlab="Log2 fold change", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=-1.5, col="black", lty=3, lwd=1.5)
abline(v= 1.5, col="black", lty=3, lwd=1.5)
abline(h=-log10(max(res$pvalue[res$padj<0.05], na.rm=TRUE)), col="black", lty=3, lwd=1.5)
with(subset(res , padj<.05), points(log2FoldChange, -log10(pvalue), pch=20,col="grey"))
with(subset(res, (log2FoldChange)> 1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, (log2FoldChange) < -1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))

# for future reference of GrCLE1-1_12 roots
res_GrCLE1_12roots <- res
GrCLE1_12roots_vsdDat <- vsd_dat
GrCLE1_12roots_vsd <- vsd


#############3
##GrCLE1-1_22_roots
countData <- data.matrix(read.delim("counts_allReps.txt", sep="\t", header=T, quote= "", fill=F, dec = ",", row.names = 1))
countData <- countData[, -c(1:3, 7:18)]
colnames(countData) <- c("GrCLE1_22roots1", "GrCLE1_22roots2", "GrCLE1_22roots3", "VDL_4roots1", "VDL_4roots2", "VDL_4roots3")
coldata <- data.frame(row.names = colnames(countData), group = rep(c("GrCLE1_22roots", "VDL"), each=3), treatment = c(rep("overexpressed", 3), rep("control", 3)))
coldata$treatment <- factor(x = coldata$treatment,levels = c("overexpressed","control"))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ treatment)
dds$treatment <- relevel(dds$treatment, ref = "control")
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ] #keeping only the genes with counts > 1
nrow(dds)
vsd <- vst(dds, blind = FALSE) #applying normalisation with VST
head(assay(vsd), 3)
vsd_dat <- assay(vst(dds, blind = FALSE))
write.table(vsd_dat, file="vst_filtered_GrCLE1_1_22roots.txt", sep = "\t", quote = FALSE)

dds <- DESeq(dds)
res <- results(dds, name = "treatment_overexpressed_vs_control")
write.csv(res, file= "DESeq2_res_GrCLE1_22roots.csv")
summary(res) #importnat to note
head(res)
mcols(res, use.names=TRUE)
resSig <- subset(res_GrCLE1_22roots, padj < 0.05)
nrow(resSig)
write.csv(resSig, file = "sig_GrCLE1_22roots.csv")
