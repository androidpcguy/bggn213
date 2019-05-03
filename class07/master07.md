Class 07: R functions and packages
================
Akshara Balachandra
April 24, 2019

## More on function writing

First we will revisit our function from last day

``` r
source('http://tinyurl.com/rescale-R')
```

``` r
x <- c(1:10, 'foobar')
is.numeric(x)
```

    ## [1] FALSE

``` r
#rescale2(x)
```

## Function practice

Write a function to identify NA elements in two vectors

Simple example where I know what the answer should be

``` r
x <- c(1,2,NA, 3, NA)
y <- c(NA,3,NA,3,4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
both.na <- is.na(x) & is.na(y) # this gives true when both are NA
which(both.na) # index of true
```

    ## [1] 3

``` r
print(paste('There are', sum(both.na), 'instances when both vectors are NA'))
```

    ## [1] "There are 1 instances when both vectors are NA"

This is my working snippet of code\!

``` r
both.na <- function(vec1, vec2) {
  sum(is.na(vec1) & is.na(vec2))
}
```

``` r
both.na(x,y)
```

    ## [1] 1

Stress testing
    now….

``` r
both.na(rep(NA, 5), c(rep(NA, 4), 3))
```

    ## [1] 4

``` r
both.na(rep(NA, 5), rep(NA,2))
```

    ## Warning in is.na(vec1) & is.na(vec2): longer object length is not a
    ## multiple of shorter object length

    ## [1] 5

Check that lengths of inputs are equal

``` r
x <- rep(NA, 5)
y <- rep(NA, 3)

length(x) != length(y) # check if lengths NOT equal
```

    ## [1] TRUE

Try both.na3() function with extra features

``` r
x <- c(1,2,NA, 3, NA)
y <- c(NA,3,NA,3,4)

both_na3(x, y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

## Calculate grades

``` r
source('grade.R')

student1 <- c(rep(100, 7), 90)
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)

print('Student 1...')
```

    ## [1] "Student 1..."

``` r
summary(grade(student1))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     100     100     100     100     100     100

``` r
print('Student 2...')
```

    ## [1] "Student 2..."

``` r
summary(grade(student2)) # hooray, both work!
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    80.0    90.0    90.0    91.0    93.5   100.0

``` r
rows = function(x) lapply(seq_len(nrow(x)), function(i) lapply(x,"[",i))

url <- "https://tinyurl.com/gradeinput"
all.grades <- read.csv(url, row.names = 1)

final.grades <- apply(all.grades, MARGIN = 1, FUN = grade)
final.mean.grades <- apply(final.grades,2, mean)
final.mean.grades
```

    ##  student-1  student-2  student-3  student-4  student-5  student-6 
    ##      91.75      82.50      84.25      84.25      88.25      89.00 
    ##  student-7  student-8  student-9 student-10 student-11 student-12 
    ##      94.00      93.75      87.75      79.00      86.00      91.75 
    ## student-13 student-14 student-15 student-16 student-17 student-18 
    ##      92.25      87.75      78.75      89.50      88.00      94.50 
    ## student-19 student-20 
    ##      82.75      82.75

Find the best-performing students
    now…

``` r
sort(final.mean.grades, decreasing = T)
```

    ## student-18  student-7  student-8 student-13  student-1 student-12 
    ##      94.50      94.00      93.75      92.25      91.75      91.75 
    ## student-16  student-6  student-5 student-17  student-9 student-14 
    ##      89.50      89.00      88.25      88.00      87.75      87.75 
    ## student-11  student-3  student-4 student-19 student-20  student-2 
    ##      86.00      84.25      84.25      82.75      82.75      82.50 
    ## student-10 student-15 
    ##      79.00      78.75

Some gene stuff

``` r
x <- df1$IDs
y <- df2$IDs

common.ids <- intersect(x, y)
```

``` r
y[y %in% x]
```

    ## [1] "gene2" "gene3"

``` r
cbind(x[x %in% y],
y[y %in% x])
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

``` r
gene_intersect2(df1, df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

``` r
merge(df1, df2, by='IDs')
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

``` r
library(shiny)
library(rsconnect)
runApp('shiny.R')
```
