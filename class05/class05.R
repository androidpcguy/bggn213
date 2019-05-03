# Author: Akshara Balachandra
# Date: 04/17/2019
# Class 5 R graphics and plots

# read data
weight.table <- read.table(file = 'bimm143_05_rstats/weight_chart.txt',
                           header = TRUE,
                           sep = '\t')

# plot scatter of age vs weight
plot(weight.table$Age, weight.table$Weight, xlab = 'Age', ylab = 'Weight', pch = 15, cex= 1.5, lwd = 2, typ = 'o')
