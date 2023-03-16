Class07: Machine Learning I
================
Patric Young

## 

In this class, we will explore clustering and dimensionalitly reduction
methods.

# K-means

Make up some data input where we know what the answer should be.

``` r
rnorm(10)
```

     [1] -0.5446561  1.8299574  1.4549246  0.6292615  0.3707596 -0.6934319
     [7]  0.2248745  1.3262848  0.1766236  0.1882518

``` r
hist(rnorm(1000000))
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-2-1.png)

``` r
hist(rnorm(1000000, -3))
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-3-1.png)

``` r
tmp <- c(rnorm(30, -3), rnorm(30, +3))
rev(tmp)
```

     [1]  3.5478508  2.0417097  2.6515231  0.8597969  0.6507474  2.4811517
     [7]  3.6141852  1.9525781  4.1334206  2.0553489  1.9983666  1.1962784
    [13]  2.0588018  2.4756248  1.9943805  2.4447155  4.0210305  2.7786669
    [19]  3.9515607  2.7304075  3.6161669  2.2722666  4.0948985  3.2076295
    [25]  2.9985298  4.4460256  3.6930090  4.7534582  2.3887685  2.6645098
    [31] -4.1314507 -3.2690350 -2.0192256 -2.3600379 -3.0247802 -2.8351145
    [37] -3.0906398 -2.8133165 -2.1852169 -4.2976923 -4.0378027 -3.9279342
    [43] -4.3511399 -2.7614300 -3.2941946 -3.6278492 -2.5235388 -3.2341771
    [49] -3.5903037 -2.1338342 -3.6737545 -2.2195428 -1.1499767 -1.7245212
    [55] -3.6470725 -2.6888128 -3.2256293 -2.6794424 -2.0574409 -1.7491168

``` r
#x <- cbind(tmp, rev(tmp))
```

``` r
tmp <- c(rnorm(30, -3), rnorm(30, +3))
x <- cbind(x=tmp, y=rev(tmp))
head(x)
```

                 x        y
    [1,] -2.160853 3.261030
    [2,] -3.517031 3.115322
    [3,] -3.966554 2.866127
    [4,] -2.866308 2.524019
    [5,] -3.407723 3.966252
    [6,] -0.513691 4.662615

Quick plot of x to see the two groups at -3, +3 and +3, -3

``` r
plot(x)
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-6-1.png)

Use the `kmeans()` function setting k to 2 and nstart=20. How many
things in here are not default. In this case its x and centers. How many
clusters we want goes inside the parenthesis.

``` r
km <- kmeans(x, centers = 2, nstart=20)
km
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1 -3.014188  3.089700
    2  3.089700 -3.014188

    Clustering vector:
     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

    Within cluster sum of squares by cluster:
    [1] 59.68348 59.68348
     (between_SS / total_SS =  90.4 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

Q1. How many points are in each cluster?

``` r
km$size
```

    [1] 30 30

Q2. What component if your result details

- cluster assignment/membership?

  ``` r
  km$cluster
  ```

       [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
      [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

- cluster center?

``` r
km$center
```

              x         y
    1 -3.014188  3.089700
    2  3.089700 -3.014188

Q3. Plot x colored by the kmeans cluster assignment and add cluster
centers as blue points

``` r
plot(x, col=c("red", "blue"))
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
plot(x, col=km$cluster)
points(km$centers, col="blue", pch=15, cex=2)
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-12-1.png)

Play with kmeans and ask for different number of clusters. kmeans
calculates distance between all the points, then updates where the
center is.

``` r
km <- kmeans(x, centers = 6, nstart=20)
plot(x, col=km$cluster)
points(km$centers, col="red", pch=6, cex=2)
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-13-1.png)

# Hierarchical Clustering

This is another very useful and widely employed clustering method which
has the advantage over `kmeans()` in that it can help reveal the
\_\_\_\_\_\_ of the true grouping in your data.

The `hclust()` function wants a distance matrix as input. Use `dist()`

kmeans(x, centers=2)

hclust(d)

cutree(hc, k=2)

``` r
d <- dist(x)
hc <- hclust(d)
hc
```


    Call:
    hclust(d = d)

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

There is a plot method for `hclust()` results.

``` r
plot(hc)
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-15-1.png)

``` r
plot(hc)
abline(h=10, col="red")
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-16-1.png)

To get my cluster membership vector, I had to “cut” my tree to yield
sub-trees or branches with all the members of a given cluster residing
on the same cut branch. The function to do this is called `cutree()`.

``` r
cutree(hc, h=10)
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
grps <- cutree(hc, h=10)
grps
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

It is often helpful to use the `k=` argument to cutree rather than the
`h=` height of cutting with `cutree()`. This will cut the tree to yield
the number of clusters you want.

``` r
grps <- cutree(hc, k=4)
grps
```

     [1] 1 2 2 2 1 1 1 1 1 2 1 2 2 2 1 2 2 1 1 1 2 2 1 1 2 1 2 2 1 2 3 4 3 3 4 3 4 4
    [39] 3 3 4 4 4 3 3 4 3 3 3 4 3 4 4 4 4 4 3 3 3 4

KM \<- kmeans(x, centers=2)

HC \<- hclust(dist(x))

GRPS \<- cutree(hc, k=2)

PCA \<-prcomp(t(x))

# Principal Component Analysis (PCA)

The base R function for PCA is called `prcomp()`. We use it to reduce
dimensionality, which is everything you can measure about a data set. We
want to visualize the most important things without losing information.
Principal components are new low dimensional axis or surfaces closest to
the observations. A line of best fit.

# PCA of UK Foods

``` r
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)
head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

**\##Q1**. How many rows and columns are in your new data frame named
`x`? What R functions could you use to answer this questions?

I could use these three below functions to find out the number of rows
and columns: `dim(x)` or `ncol(x)` or `nrow(x)`

``` r
rownames(x) <- x[,1]
x <- x[,-1]
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
dim(x)
```

    [1] 17  4

``` r
x <- read.csv(url, row.names=1)
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

**Q2.** Which approach to solving the ‘row-names problem’ mentioned
above do you prefer and why? Is one approach more robust than another
under certain circumstances?

I believe the second approach to the “row-names problem” is the superior
method because it is more efficient. Compared to the first method, which
began deleting columns upon running because of the `x, -1`, the second
method does not delete anything. Furthermore, there are less lines of
text involved.

``` r
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-24-1.png)

**Q3**: Changing what optional argument in the above **barplot()**
function results in the following plot?

If we change `besides=T` to `besides=F`, the following plot is created.
The bars for each country are wider and furthermore stacked on top of
one another.

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-25-1.png)

**Q5**: Generating all pairwise plots may help somewhat. Can you make
sense of the following code and resulting figure? What does it mean if a
given point lies on the diagonal for a given plot?

The following code further transforms the data and maps the data for all
four countries onto pairwise plots. This allows us to see the
distribution of single variables and the relationship between multiple
variables. If a given point lies on the diagonal for a given plot, it
means that that particular point has association with other points and
clusters. It is part of a trend.

``` r
pairs(x, col=rainbow(10), pch=16)
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-26-1.png)

**Q6**. What is the main differences between N. Ireland and the other
countries of the UK in terms of this data-set?

N. Ireland has a dramatically lower standard deviation compared to the
other countries of the UK in terms of this data-set. Its STDEV is less
than 1 while the other countries’ STDEV is over 50. N. Ireland also has
a Proportion of Variance of 0 while the other countries have values
greater than 0.

``` r
pca <- prcomp( t(x))
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

A “PCA plot” (aka “Score Plot, PC1vPC2 plot, etc) is a plot comparing
two plots

``` r
pca$x
```

                     PC1         PC2         PC3           PC4
    England   -144.99315    2.532999 -105.768945  2.842865e-14
    Wales     -240.52915  224.646925   56.475555  7.804382e-13
    Scotland   -91.86934 -286.081786   44.415495 -9.614462e-13
    N.Ireland  477.39164   58.901862    4.877895  1.448078e-13

``` r
plot(pca$x[,1], pca$x[,2], col=c("orange", "red", "blue", "darkgreen"), pch=18)
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-29-1.png)

**Q7**. Complete the code below to generate a plot of PC1 vs PC2. The
second line adds text labels over the data points

``` r
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500))
text(pca$x[,1], pca$x[,2], colnames(x))
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-30-1.png)

**Q8.** Customize your plot so that the colors of the country names
match the colors in our UK and Ireland map and table at start of this
document.

``` r
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500))
text(pca$x[,1], pca$x[,2], colnames(x),col=c("orange", "magenta", "blue", "darkgreen"), pch=18 )
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-31-1.png)

``` r
v <- round( pca$sdev^2/sum(pca$sdev^2) * 100 )
v
```

    [1] 67 29  4  0

``` r
z <- summary(pca)
z$importance
```

                                 PC1       PC2      PC3          PC4
    Standard deviation     324.15019 212.74780 73.87622 4.188568e-14
    Proportion of Variance   0.67444   0.29052  0.03503 0.000000e+00
    Cumulative Proportion    0.67444   0.96497  1.00000 1.000000e+00

``` r
barplot(v, xlab="Principal Component", ylab="Percent Variation")
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-34-1.png)

``` r
par(mar=c(10, 3, 0.35, 0))
barplot( pca$rotation[,1], las=2 )
```

![](BIMM143-Lab-07_files/figure-commonmark/unnamed-chunk-35-1.png)
