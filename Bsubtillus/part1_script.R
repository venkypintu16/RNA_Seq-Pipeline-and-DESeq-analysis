###################################
########### PART 1 ################
###################################

#Load Data (RPKM File)
setwd("~/applied_comp_gen/project2/part1/")
dat_1 <- read.delim2(file = "rpkm_outfile_B.subtillis.txt")
rpkm_1 <- dat_1[,5:6]

#Scatter Plot and 2D Density Plot (ggplot)
library(ggplot2)

rpkm_1$DepletedRPKM <- as.numeric(as.character(rpkm_1$DepletedRPKM)) #Convert values to numeric when N/A (These are taken as null values) are present
rpkm_1$TotalRPKM <- as.numeric(as.character(rpkm_1$TotalRPKM)) #Convert values to numeric when N/A (These are taken as null values) are present

ggplot(log2(rpkm_1), aes(x = DepletedRPKM, y = TotalRPKM)) + 
  geom_point() + geom_smooth(method =lm, se = F) + geom_density_2d()