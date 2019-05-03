---
author: Akshara Balachandra
title: Class 08 Unsupervised Learning
date: 04/26/19
output: github_document
---

## KNN clustering examples

```{r}
tmp <- c(rnorm(30, -3), rnorm(30, 3))
x <- cbind(x = tmp, y = rev(tmp))
plot(x)
```

```{r}
# running kmeans
ans <- kmeans(x, 2, 5)
ans
print(paste("There are", ans$size[1], "points in cluster 1 and",
	    ans$size[2], "points in cluster 2"))

```



### Questions

1. How many points are in each cluster?
	- 30 each
2. What component of result object details
	- cluster size: `size`
	- cluster assigment/membership: `cluster`
	- cluster center: `centers`
3. Plot colored kmeans clustered assigment

```{r}
plot(x, typ = 'p', col = ans$cluster)
points(ans$centers, col = 'blue', pch = 18, cex = 3)
```

## Heirarchical clustering
No need to specify number of clusters beforehand....

```{r}
d <- dist(x) # get distance matrix
hc <- hclust(d)
hc # not very useful to print object
```

So, let's plot instead...

```{r}
plot(hc)
abline(h = 10, col = 'red')
gp1 <- cutree(hc, k = 2) # cut to give two clusters
gp2 <- cutree(hc, h = 10) # cut at height 10
gp3 <- cutree(hc, k = 3) # cut to give three clusters

table(gp2, gp3)
```

### Example HClustering

```{r}
# Step 1. Generate some example data for clustering
x <- rbind(
	    matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2), # c1
	    matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
	    matrix(c(rnorm(50, mean = 1, sd = 0.3), # c3
		     rnorm(50, mean = 0, sd = 0.3)), ncol = 2))

colnames(x) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x)
# Step 3. Generate colors for known clusters
# (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x, col=col)
```

#### Questions
Use dist(), hclust(), plot, cutree functions to return 2 and 3 clusters

```{r}
distance <- dist(x)
hc <- hclust(distance)
plot(hc)

clust2 <- cutree(hc, k = 2)
clust3 <- cutree(hc, k = 3)

table(clust2, clust3)
```

Plotting the three clusters found...

```{r}
plot(x, col = clust3, main = 'HClustering clusters')
plot(x, col = col, main = 'Original clusters')
```

## PCA
PCA with `prcomp` from `{base}` R package.

```{r}
mydata <- read.csv('https://tinyurl.com/expression-CSV', row.names = 1)
head(mydata)
dim(mydata)
```

There are `r nrow(mydata)` genes in the data, with `r ncol(mydata)` experimental
conditions.

```{r}
pca <- prcomp(t(mydata), scale = T)
summary(pca)

attributes(pca)
```

Let's make our first PCA plot.

```{r}

## Precent variance is often more informative to look at
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pca.var.per

# make some colors yo
mycols <- c(rep('red', 5), rep('blue', 5))

xlab = paste('PC1 (', pca.var.per[1], '%)', sep = '')
ylab = paste('PC2 (', pca.var.per[2], '%)', sep = '')
plot(pca$x[,1], pca$x[,2], xlab = xlab, ylab = ylab, col = mycols)
text(pca$x[,1], pca$x[,2], colnames(mydata))
```

## PCA on UK foods dataset

Read the data...
```{r}
food.data <- read.csv('data/UK_foods.csv', row.names = 1)
head(food.data)
```
There are `r nrow(food.data)` rows and `r ncol(food.data)` columns in the
dataset. You should read the data in properly instead of setting manually
the first column as row names.

Exploratory data analysis with barplots!

```{r}
barplot(as.matrix(food.data), beside = T, col = rainbow(nrow(food.data)))

# weirdo stacked plot...
barplot(as.matrix(food.data), col = rainbow(nrow(food.data)))
```

You need `beside = T` to get horizontally "stacked" bars. We can
also generate pairwise plots between every pair of countries and
visualize.

```{r}
pairs(food.data, col = rainbow(nrow(food.data)), pch = 16)
```

Northern Ireland has blue and orange points that are way off the diagonal
unlike the other countries.

PCA:

```{r}
pca <- prcomp(t(food.data))
summary(pca)

plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500))
#            ENGLAND	WALES	SCOT	NI
colors <- c('orange',	'red',	'blue',	'green')
text(pca$x[,1], pca$x[,2], colnames(food.data), col = colors)
```

Calculating/plotting percent variation captured by each PC...

```{r}
sdev <- pca$sdev
v <- round(sdev ^ 2 / sum(sdev^2) * 100)

barplot(v, xlab = 'Principal component', ylab = 'Percent variation')
```

### PC Weights
Calculating loading scores....

```{r}
## Lets focus on PC1 as it accounts for > 90% of variance
par(mar=c(10, 3, 5, 0))
barplot( pca$rotation[,1], las=2, main = 'PC1 loading scores')
```

```{r}
par(mar=c(10, 3, 5, 0))
barplot( pca$rotation[,2], las=2, main = 'PC2 loading scores')
```

PC2 tells mainly about potatoes and soft drinks (largest weights by absolute
value)

Biplotting...

```{r}
biplot(pca)
```




