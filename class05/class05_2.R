# bar plotssss


# read data
features <- read.table(file = 'bimm143_05_rstats/feature_counts.txt',
                       header = TRUE,
                       sep = '\t')

par(mar = c(5.1, 12, 4.1, 2.1))
# plotssss
barplot(features$Count,
        horiz = TRUE,
        xlab = 'Feature Count',
        names.arg = features$Feature,
        las = 1,
        main = 'Number of features in mouse genome')

