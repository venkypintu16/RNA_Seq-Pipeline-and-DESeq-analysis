# Normalization and Potting part.

setwd("C:/Users/DELL/OneDrive/Documents/Applied computational genomics/Midterm Project/Midterm/ECOLI PROJECT FILES/ECOLI PROJECT FILES")

df <- read.table("C:/Users/DELL/OneDrive/Documents/Applied computational genomics/Midterm Project/Midterm/ECOLI PROJECT FILES/ECOLI PROJECT FILES/Ecoli_feature_count.txt",
                      header = T, sep = '\t', row.names = 1)

RPKM = read.table("C:/Users/DELL/OneDrive/Documents/Applied computational genomics/Midterm Project/Midterm/ECOLI PROJECT FILES/ECOLI PROJECT FILES/Ecoli_feature_count.RPMK_values.txt")

colnames(RPKM) <- c("Gene_ID", "Read_Counts_Depleted", "Read_Counts_Total", "Gene_Length", "DepletedRPKM", "TotalRPKM")

# Assuming 'data' is your data frame
RPKM <- RPKM[-1, ]

# Check the data type of 'DepletedRPKM' and 'TotalRPKM'
class(RPKM$DepletedRPKM)
class(RPKM$TotalRPKM)

# Convert the columns to numeric if they are not already

RPKM$DepletedRPKM <- as.numeric(as.character(RPKM$DepletedRPKM))
RPKM$TotalRPKM <- as.numeric(as.character(RPKM$TotalRPKM))

#Normal Scatter Plot

plot(log2(RPKM$DepletedRPKM), log2(RPKM$TotalRPKM))

# Scatter plot using ggplot library.
library(ggplot2)
RPKM_c <- data.frame(c1 = as.numeric( log2(RPKM[,5] )),
                     c2 = as.numeric( log2(RPKM[,6] )) )
ggplot(RPKM_c, aes(x = c1, y = c2)) +
  geom_point() +
  xlab("Depleted") +
  ylab("Total") +
  geom_smooth(method = lm, se = F)

# Density plot.
ggplot(RPKM_c, aes(x = c1, y = c2)) +
  geom_point() +
  xlab("Depleted") +
  ylab("Total") +
  geom_density_2d()


# Scatter plot with linear regression line & Density plot
ggplot(RPKM_c, aes(x = c1, y = c2)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  geom_density_2d() +  # Adding 2D density plot
  xlab("Depleted") +
  ylab("Total")


#Boxplot

boxplot(log2(RPKM[,5]),outline = F, xlab = "Depleted", ylab = "mRNA expression")
boxplot(log2(RPKM[,6]), outline = F, xlab = "Total", ylab = "mRNA expression")


# Correlation analysis & Heatmaps.
RPKM_updated = RPKM[,5:6 ]
corF <- cor(RPKM_updated)
corF

heatmap(corF)

library(gplots)
heatmap.2(corF, Colv= F, dendrogram = 'row',
          trace = 'n',      ## hist plot line inside heatmap
          denscol = 'navy'  ## hist plot line in color key
)


# Histogram Plots,

df <- read.table("C:/Users/DELL/OneDrive/Documents/Applied computational genomics/Midterm Project/Midterm/ECOLI PROJECT FILES/ECOLI PROJECT FILES/Ecoli_feature_count.txt",
                 header = T, sep = '\t', row.names = 1)
df_DE = df[,6:7]

print( apply(df_DE, 2, function(x) quantile(as.numeric(x))) )

par(mfrow = c(1,2), mar = c(4,4,1,1))

# histplot before trimming
hist.plots <- apply(df_DE, 2, function(x) {hist(log2(x), breaks = 50)})

# Data filtering step.
expressed.ids <- apply(df_DE, 1, function(x) any(x > 20))
dfExp <- df_DE[expressed.ids, ]

par(mfrow = c(1,2), mar = c(4,4,1,1))
print( apply(dfExp, 2, function(x) quantile(as.numeric(x))) )

# hisplot after trimming.
hist.plots <- apply(dfExp, 2, function(x) {hist(log2(x), breaks = 50)})

