---
author: Akshara Balachandra
title: Class 09 Unsupervised Learning
date: 05/01/19
output: rmarkdown::github_document
---

## Unsupervised learning analysis of breast cancer cells

```{r}
wisc.df <- read.csv('data/WisconsinCancer.csv', header = T)
nrow(wisc.df)
ncol(wisc.df)
table(wisc.df$diagnosis)
```

Number of features containing "_mean"

```{r}
sum(grepl("_mean", colnames(wisc.df)))
length(grep("_mean", colnames(wisc.df)))
```

Converting dataframe to matrix taking only columns 3 to 32...

```{r}
final.col <- ncol(wisc.df) - 1
wisc.data <- as.matrix(wisc.df[,3:final.col])
row.names(wisc.data) <- wisc.df$id

diagnosis <- wisc.df$diagnosis

#head(wisc.data)
```

## PCA

Check column means and standard deviation of features.

```{r}
colMeans(wisc.data)
apply(wisc.data, 2, sd)
```

Column means are very different between different features, so I need
to scale the data before pca.

```{r}
wisc.pr <- prcomp(wisc.data, scale = T)

# summary of results
summary(wisc.pr)
```

Q4: 44% of original variance from PC-1
Q5: 3 PCs (72.6% variance)
Q6: 7 PCs (91.0% variance)

### Interpreting PCA results

```{r}
biplot(wisc.pr)
```

This plot is pretty much useless, we can't see anything. So, let's create
a scatter plot of the first two PCs.

```{r}
plot(x = wisc.pr$x[,1], y = wisc.pr$x[,2], col = diagnosis,
     xlab = "PC1", ylab = "PC2")
```

```{r}
plot(x = wisc.pr$x[,1], y = wisc.pr$x[,3], col = diagnosis,
     xlab = "PC1", ylab = "PC3")
```

PC1 vs PC2 gives a more clear-cut boundary between the two groups
than PC1 vs PC3 because PC2 explains more variance in the data
than PC3.


### Variance explained

```{r}
pc.var <- wisc.pr$sdev ^ 2

# variance explained by each pc: pve
total_var <- sum(pc.var)

pve <- pc.var / total_var

plot(pve, xlab = "Principal Component",
     ylab = 'Proportion of Variance Explained',
     ylim = c(0, 1), type = 'o')

# barplot version
barplot(pve, ylab = "Percent of Variance Explained",
	names.arg=paste0("PC", 1:length(pve)), las=2, axes = FALSE)

axis(2, at=pve, labels=round(pve,2 )*100)
```

### Communicating PCA results

Loading values for PCs

```{r}
sort(abs(wisc.pr$rotation[, 1]))
```

## Heirarchical Clustering

Scale and find distances of data...

```{r}
data.scaled <- scale(wisc.data)
data.dist <- dist(data.scaled)
```

Create hclust model...

```{r}
wisc.hclust <- hclust(data.dist, method = 'complete')
plot(wisc.hclust)
abline(h = 19, col = 'red', lty = 2)

```

### Selecting number of clusters

```{r}
wisc.hclust.clusters <- cutree(wisc.hclust, k = 4)
table(wisc.hclust.clusters, diagnosis)
```

4 clusters gives best separation.

## Combining methods (pca + hclust)

Do hclust on pca data


```{r}
wisc.pr.hclust <- hclust( dist( wisc.pr$x[, 1:7] ), method = 'ward.D2')
plot(wisc.pr.hclust)

grps <- cutree(wisc.pr.hclust, k = 2)
table(grps, diagnosis)

plot(wisc.pr$x[, 1:2], col = grps)
plot(wisc.pr$x[, 1:2], col = diagnosis)
```

Recolor the hclust groupings PC plot...

```{r}
g <- as.factor(grps)
levels(g)

g <- relevel(g, 2)
levels(g)

# plot with reordered levels
plot(wisc.pr$x[, 1:2], col = g)
```

## Precition

Project new data onto PCA basis vectors

```{r}
url <- 'https://tinyurl.com/new-samples-CSV'
new <- read.csv(url)
npc <- predict(wisc.pr, newdata = new)
head(npc)


# plot these patients with existing data...
plot(wisc.pr$x[, 1:2], col = g)
points(npc[, 1], npc[, 2], col = 'blue', pch = 16, cex=3)
text(npc[, 1], npc[, 2], col = 'white')
```

Patient 2 should be followed up on because they have features that align them
closely wiht other cancer patients.

## Optional (PCA on protein data)







