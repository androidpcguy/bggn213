---
title: 'Class 6: R Functions'
author: "Akshara Balachandra"
date: "April 19, 2019"
output: rmarkdown::github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = T)
```

## Overview
Today we'll be going over **R functions**, but will start off with some
**file reading first**.

```{r}
plot(1:10, type = 'l', col = 'blue')
```


# test1 file
```{r}
test1 <- read.table('test1.txt', header = T, sep = ',')
test1
```

# test2 file
```{r}
test2 <- read.table('./test2.txt', sep = '$', header = T)
test2
```

# test3 file
```{r}
test3 <- read.table('./test3.txt', header = F)
test3
```





## First function
```{r}
add <- function(x, y = 1) {
  # the body
  x + y
}
```


```{r}
add(c(1, 3, 5))
```

```{r}
# add(1, 3, 5)
# add(x = 1, y = "barry")
```


```{r}
rescale <- function(x) {
  rng <- range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}

```

```{r}
rescale2 <- function(x, na.rm = T) {
  rng <- range(x, na.rm = na.rm)
  (x - rng[1]) / (rng[2] - rng[1])
}
```

```{r}
x <- c(1:10, NA)
rng <- range(x)
(x - rng[1]) / (rng[2] - rng[1])
rng
```


```{r}
print(rescale2(c(1:10, NA)))
print(rescale2(c(1:10, NA), na.rm = F))
```



Another example extension...

```{r}
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
  rng <-range(x, na.rm=na.rm)

  print("Hello")
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  print("is it me you are looking for?")
  if(plot) {
    plot(answer, typ="b", lwd=4)
    print("Don't sing please!!!!")
  }
  print("I can see it in ...")
  return(answer)
}
```

```{r}
rescale3(1:10)
```

```{r}
rescale3(1:10, plot = T)
```



## Practical protein interaction example

```{r}
source('./bio3d_example.R')
prot_interact('4AKE')
prot_interact('1AKE')
prot_interact('1E4Y')
```






















