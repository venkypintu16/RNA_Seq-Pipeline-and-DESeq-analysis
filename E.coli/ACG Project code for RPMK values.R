# Load required libraries
install.packages("data.table")
library(data.table)

# Read feature counts data
feature_counts <- fread("E:/feature_count.txt")

# Extract relevant columns (gene ID, read counts, and gene lengths)
gene_data <- feature_counts[, c("Geneid", "depleted", "total", "Length")]

# Rename columns for clarity
colnames(gene_data) <- c("Gene_ID", "Read_Counts_Depleted", "Read_Counts_Total", "Gene_Length")

# Calculate total mapped reads for depleted and total columns
total_mapped_reads_depleted <- sum(gene_data$Read_Counts_Depleted)
total_mapped_reads_total <- sum(gene_data$Read_Counts_Total)

# Calculate RPKM for depleted and total columns
rpkm_values_depleted <- (gene_data$Read_Counts_Depleted / (gene_data$Gene_Length / 1000)) / (total_mapped_reads_depleted / 1e6)
rpkm_values_total <- (gene_data$Read_Counts_Total / (gene_data$Gene_Length / 1000)) / (total_mapped_reads_total / 1e6)

# Add RPKM values as new columns to the gene data
gene_data$DepletedRPKM <- rpkm_values_depleted
gene_data$TotalRPKM <- rpkm_values_total

# Print or save the results
print(gene_data)

# If you want to save the results to a file:
# write.table(gene_data, "rpkm_values_with_columns.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(gene_data, file="E:/Ecoli_feature_count.RPMK_values.txt", sep="\t", quote=FALSE, row.names=FALSE)
#This script reads the feature counts data, calculates the RPKM values, adds them as new