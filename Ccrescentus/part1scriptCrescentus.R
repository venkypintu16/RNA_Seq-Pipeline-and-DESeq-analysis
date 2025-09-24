
# Alternatively, if you already have a data frame and want to assign column names
# You can use the colnames() function like this:
colnames(RPKM) <- c("Gene_ID", "Read_Counts_Depleted", "Read_Counts_Total", "Gene_Length", "DepletedRPKM", "TotalRPKM")

# Assuming 'data' is your data frame
RPKM <- RPKM[-1, ]

# Check the data type of 'DepletedRPKM' and 'TotalRPKM'
class(RPKM$DepletedRPKM)
class(RPKM$TotalRPKM)

# Convert the columns to numeric if they are not already
RPKM$DepletedRPKM <- as.numeric(as.character(RPKM$DepletedRPKM))
RPKM$TotalRPKM <- as.numeric(as.character(RPKM$TotalRPKM))


#STARTHERE

# Now, you can plot the log2-transformed values
plot(log2(RPKM$DepletedRPKM), log2(RPKM$TotalRPKM))

# GGPLOT.
library(ggplot2)
RPKM_c <- data.frame(depleted = as.numeric( log2(RPKM[,5] )),
                     total = as.numeric( log2(RPKM[,6] )) )
ggplot(RPKM_c, aes(x = depleted, y = total)) +
  geom_point() +
  geom_smooth(method = lm, se = F)

ggplot(RPKM_c, aes(x = depleted, y = total)) +
  geom_point() +
  geom_density_2d()

# Correlation bigger font
#RPKM_updated = RPKM[,5:6 ]
#corF <- cor(RPKM_updated)
#corF


#boxplot
# Read the data
norData <- read.table("~/Documents/Computational Genomics/RPKM.txt", row.names = 1, header = TRUE)

# Subset columns 5 and 6
subset_data <- norData[, c("DepletedRPKM", "TotalRPKM")]

# Create boxplot
boxplot(log2(subset_data), outline = FALSE)

