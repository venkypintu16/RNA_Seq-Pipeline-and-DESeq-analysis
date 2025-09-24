###################################
########### PART 2 ################
###################################

#Load Data (RPKM File)
setwd("~/applied_comp_gen/project2/part2/")
dat_2 <- read.delim2(file = "rpkm_outfile_part2_ribo0_vs_diy.txt")
rpkm_2 <- dat_2[,9:14]

#Scatter Plot and 2D Density Plot (ggplot)
library(ggplot2)

rpkm_2 <- as.data.frame(apply(rpkm_2, 2, function(x) {as.numeric(as.character(x))})) #Convert values to numeric when N/A (These are taken as null values) are present

CTRL <- ggplot(log2(rpkm_2), aes(x = RiboZero_Ctrl_RPKM, y = DIY_Ctrl_RPKM )) + geom_point() + geom_smooth(method =lm, se = F) + geom_density_2d()
RIF <- ggplot(log2(rpkm_2), aes(x = RiboZero_Rif_RPKM, y = DIY_Rif_RPKM )) + geom_point() + geom_smooth(method =lm, se = F) + geom_density_2d()
CHLO <- ggplot(log2(rpkm_2), aes(x = RiboZero_Chlo_RPKM, y = DIY_Chlo_RPKM )) + geom_point() + geom_smooth(method =lm, se = F) + geom_density_2d()

plot(CTRL)
plot(RIF)
plot(CHLO)

#Differential Gene Expression Analysis
#Load Data
input_raw_count <- read.table(file = "feature_counts_part2_ribo0_vs_diy.txt", header = T, sep = '\t')
input_count <- input_raw_count[, 7:12]


#Perform DESeq2
library(DESeq2)
info_table <- data.frame(row.names = colnames(input_count),
                         Treatment = c("Control", "Rifampicin", "Chloramphenicol", "Control", "Rifampicin", "Chlorampheni"),
                         kit_type = c("RiboZero", "RiboZero", "RiboZero", "DIY", "DIY", "DIY"))

keep <- rowSums(input_count) >= 30 #Applying filter
filter_input_count <- input_count[keep, ] #Applying filter
dim(filter_input_count) #Applying filter

dds <- DESeqDataSetFromMatrix(filter_input_count, colData = info_table, design = ~ kit_type) #Construct a DESeq Data Set object

dds$kit_type <- relevel(dds$kit_type, ref = "RiboZero") #Set the factor level (set the reference)

dds <- DESeq(dds) #Run DESeq

norCounts <- counts(dds, normalized = TRUE) #Normalizing read counts
res <- results(dds) #DESeq Results
summary(res)

resSig <- res[ which(res$padj < 0.01), ] #Extract genes with adjusted P-value < 0.01
resSig.sorted.df <- as.data.frame(resSig[ order( resSig$log2FoldChange ), ])

#Heatmap
library(RColorBrewer)
library(gplots)

sigNorData <- norCounts[rownames(norCounts) %in% rownames(resSig),]
hmcol <-  colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatMapDF_nor <- t( apply(sigNorData, 1, function(x){(x-mean(x))/sd(x)}) )
colnames(heatMapDF_nor) <- c('RZ_Ctrl', 'RZ_Rif', 'RZ_Chlo', 'DIY_Ctrl', 'DIY_Rif', 'DIY_Chlo')
heatmap.2(heatMapDF_nor, col = hmcol, trace = 'none', margins = c(10, 10), labRow = F)

#Volcano Plot
res_plot      <- data.frame( res )
res_plot$col  <- 'gray40'

res_plot$col[res_plot$log2FoldChange > 1 & res_plot$padj < 0.01] <- 'red'
res_plot$col[res_plot$log2FoldChange < -1 & res_plot$padj < 0.01] <- 'cornflowerblue'

par(mar = c(5,5,1,1))
plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)', xlim = c(-8,8),)
     