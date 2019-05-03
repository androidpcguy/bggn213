# colorsssss

mf_cnts <- read.table(file = 'bimm143_05_rstats/male_female_counts.txt',
                      header = TRUE,
                      sep = '\t')

barplot(mf_cnts$Count,
        names.arg = mf_cnts$Sample,
        las = 2,
        col = topo.colors(10))


up_down <- read.table(file = 'bimm143_05_rstats/up_down_expression.txt',
                      header = T,
                      sep = '\t')
table(up_down$State)

# change the color of each state
palette(value = c('red', 'gray', 'green'))

plot(up_down$Condition1, up_down$Condition2, xlab = "Expression C1",
     ylab = "Expression C2", col = up_down$State, pch = 15)

