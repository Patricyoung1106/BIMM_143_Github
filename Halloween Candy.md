Halloween Candy Mini Project
================
Patric Young

# 1. Importing candy data

``` r
candy_file <- "candy-data.csv"
candy = read.csv("candy-data.csv", row.names=1)
head(candy)
```

                 chocolate fruity caramel peanutyalmondy nougat crispedricewafer
    100 Grand            1      0       1              0      0                1
    3 Musketeers         1      0       0              0      1                0
    One dime             0      0       0              0      0                0
    One quarter          0      0       0              0      0                0
    Air Heads            0      1       0              0      0                0
    Almond Joy           1      0       0              1      0                0
                 hard bar pluribus sugarpercent pricepercent winpercent
    100 Grand       0   1        0        0.732        0.860   66.97173
    3 Musketeers    0   1        0        0.604        0.511   67.60294
    One dime        0   0        0        0.011        0.116   32.26109
    One quarter     0   0        0        0.011        0.511   46.11650
    Air Heads       0   0        0        0.906        0.511   52.34146
    Almond Joy      0   1        0        0.465        0.767   50.34755

- **Q1**. How many different candy types are in this dataset?

  ``` r
  nrow(candy)
  ```

      [1] 85

  There are 85 different candy types in this data set.

- **Q2**. How many fruity candy types are in the data set?

  ``` r
  table(candy$fruit)
  ```


       0  1 
      47 38 

There are 38 fruity candy types are in this data set.

# 2. What is your favorate candy?

- **Q3**. What is your favorite candy in the data set and what is it’s
  `winpercent` value?

``` r
candy["Twix", ]$winpercent
```

    [1] 81.64291

Favorite candy: Twix

Winpercent: 81.64%

- **Q4**. What is the `winpercent` value for “Kit Kat”?

  ``` r
  candy["Kit Kat", ]$winpercent
  ```

      [1] 76.7686

  76.77%

- **Q5**. What is the `winpercent` value for “Tootsie Roll Snack Bars”?

``` r
candy["Tootsie Roll Snack Bars", ]$winpercent
```

    [1] 49.6535

49.65%

``` r
library("skimr")
skim(candy)
```

|                                                  |       |
|:-------------------------------------------------|:------|
| Name                                             | candy |
| Number of rows                                   | 85    |
| Number of columns                                | 12    |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |       |
| Column type frequency:                           |       |
| numeric                                          | 12    |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |       |
| Group variables                                  | None  |

Data summary

**Variable type: numeric**

| skim_variable    | n_missing | complete_rate |  mean |    sd |    p0 |   p25 |   p50 |   p75 |  p100 | hist  |
|:-----------------|----------:|--------------:|------:|------:|------:|------:|------:|------:|------:|:------|
| chocolate        |         0 |             1 |  0.44 |  0.50 |  0.00 |  0.00 |  0.00 |  1.00 |  1.00 | ▇▁▁▁▆ |
| fruity           |         0 |             1 |  0.45 |  0.50 |  0.00 |  0.00 |  0.00 |  1.00 |  1.00 | ▇▁▁▁▆ |
| caramel          |         0 |             1 |  0.16 |  0.37 |  0.00 |  0.00 |  0.00 |  0.00 |  1.00 | ▇▁▁▁▂ |
| peanutyalmondy   |         0 |             1 |  0.16 |  0.37 |  0.00 |  0.00 |  0.00 |  0.00 |  1.00 | ▇▁▁▁▂ |
| nougat           |         0 |             1 |  0.08 |  0.28 |  0.00 |  0.00 |  0.00 |  0.00 |  1.00 | ▇▁▁▁▁ |
| crispedricewafer |         0 |             1 |  0.08 |  0.28 |  0.00 |  0.00 |  0.00 |  0.00 |  1.00 | ▇▁▁▁▁ |
| hard             |         0 |             1 |  0.18 |  0.38 |  0.00 |  0.00 |  0.00 |  0.00 |  1.00 | ▇▁▁▁▂ |
| bar              |         0 |             1 |  0.25 |  0.43 |  0.00 |  0.00 |  0.00 |  0.00 |  1.00 | ▇▁▁▁▂ |
| pluribus         |         0 |             1 |  0.52 |  0.50 |  0.00 |  0.00 |  1.00 |  1.00 |  1.00 | ▇▁▁▁▇ |
| sugarpercent     |         0 |             1 |  0.48 |  0.28 |  0.01 |  0.22 |  0.47 |  0.73 |  0.99 | ▇▇▇▇▆ |
| pricepercent     |         0 |             1 |  0.47 |  0.29 |  0.01 |  0.26 |  0.47 |  0.65 |  0.98 | ▇▇▇▇▆ |
| winpercent       |         0 |             1 | 50.32 | 14.71 | 22.45 | 39.14 | 47.83 | 59.86 | 84.18 | ▃▇▆▅▂ |

- **Q6**. Is there any variable/column that looks to be on a different
  scale to the majority of the other columns in the dataset?

  The winpercent variable is on a different scale to the majority of the
  other columns in the data set. While the other variables are on a
  scale of 0 to 1, winpercent is on a 0 to 100 percent scale.

- **Q7**. What do you think a zero and one represent for the
  `candy$chocolate` column?

  A zero represents no “vote” was casted for Chocolate whereas a one
  represents that a vote was casted.

- **Q8**. Plot a histogram of `winpercent` values

``` r
library("ggplot2")
```

``` r
ggplot(candy, aes(x = winpercent)) +
geom_histogram()
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-9-1.png)

- **Q9**. Is the distribution of `winpercent` values symmetrical?

  The distribution of `winpercent` values are not symmetrical.

- **Q10**. Is the center of the distribution above or below 50%?

The center of distribution is slightly below 50%

- **Q11**. On average is chocolate candy higher or lower ranked than
  fruit candy?

  On average, Chocolate candy is ranked higher than fruit candy

- **Q12**. Is this difference statistically significant?

  Since the p-value falls below the significance level, the result is
  statistically significant.

# 3. Overall Candy Rankings

- **Q13**. What are the five least liked candy types in this set?

``` r
head(candy[order(candy$winpercent),], n=5)
```

                       chocolate fruity caramel peanutyalmondy nougat
    Nik L Nip                  0      1       0              0      0
    Boston Baked Beans         0      0       0              1      0
    Chiclets                   0      1       0              0      0
    Super Bubble               0      1       0              0      0
    Jawbusters                 0      1       0              0      0
                       crispedricewafer hard bar pluribus sugarpercent pricepercent
    Nik L Nip                         0    0   0        1        0.197        0.976
    Boston Baked Beans                0    0   0        1        0.313        0.511
    Chiclets                          0    0   0        1        0.046        0.325
    Super Bubble                      0    0   0        0        0.162        0.116
    Jawbusters                        0    1   0        1        0.093        0.511
                       winpercent
    Nik L Nip            22.44534
    Boston Baked Beans   23.41782
    Chiclets             24.52499
    Super Bubble         27.30386
    Jawbusters           28.12744

The five least liked candy types in this set are Nik L Nip, Boston Baked
Beans, Chiclets, Super Bubble, and Jawbusters.

- **Q14**. What are the top 5 all time favorite candy types out of this
  set?

The top 5 all time favorite candy types out of this set are Reeses
Peanut Butter Cup, Reeses Minatures, Twix, Kit Kat, and Reeses Pieces.

**Q15**. Make a first barplot of candy ranking based on `winpercent`
values.

``` r
library(ggplot2)

ggplot(candy) + 
  aes(winpercent, rownames(candy)) +
  geom_col()
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-11-1.png)

- **Q16**. This is quite ugly, use the `reorder()` function to get the
  bars sorted by `winpercent`?

``` r
library(ggplot2)

ggplot(candy) + 
  aes(winpercent, reorder(rownames(candy),winpercent)) +
  geom_col()
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
my_cols=rep("black", nrow(candy))
my_cols[as.logical(candy$chocolate)] = "chocolate"
my_cols[as.logical(candy$bar)] = "brown"
my_cols[as.logical(candy$fruity)] = "pink"
```

**Q17**. What is the worst ranked chocolate candy?

The worst ranked chocolate candy are Sixlets.

**Q18**. What is the best ranked fruity candy?

The best ranked fruity candy are Starbursts.

# 4. Taking a look at pricepercent

``` r
ord <- order(candy$pricepercent, decreasing = TRUE)
head( candy[ord,c(11,12)], n=5 )
```

                             pricepercent winpercent
    Nik L Nip                       0.976   22.44534
    Nestle Smarties                 0.976   37.88719
    Ring pop                        0.965   35.29076
    Hershey's Krackel               0.918   62.28448
    Hershey's Milk Chocolate        0.918   56.49050

**Q19**. Which candy type is the highest ranked in terms of `winpercent`
for the least money - i.e. offers the most bang for your buck?

Hershey’s Krackel is the highest `winpercent` candy (62.28%) with the
lowest `pricepercent`(0.918%).

**Q20**. What are the top 5 most expensive candy types in the dataset
and of these which is the least popular?

The top 5 most expensive candy types in the data set are Nik L Nip,
Nestle Smarties, Ring Pop, Hershey’s Krackel, and Hershey’s Milk
Chocolate. The least popular of these are the Nik L Nips.

- **Q21**. Make a bar plot again with `geom_col()` this time using
  `pricepercent` and then improve this step by step, first ordering the
  x-axis by value and finally making a so called “dot chat” or
  “lollipop” chart by swapping `geom_col()` for `geom_point()` +
  `geom_segment()`.

``` r
ggplot(candy) +
  aes(pricepercent, reorder(rownames(candy), pricepercent)) +
  geom_col()
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-15-1.png)

``` r
# Make a lollipop chart of pricepercent
ggplot(candy) +
  aes(pricepercent, reorder(rownames(candy), pricepercent)) +
  geom_segment(aes(yend = reorder(rownames(candy), pricepercent), 
                   xend = 0), col="gray40") +
    geom_point()
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-16-1.png)

**Q22**. Examining this plot what two variables are anti-correlated
(i.e. have minus values)?

Chocolate and Fruity candies are the two most anti-correlated variables.

**Q23**. Similarly, what two variables are most positively correlated?

Chocolate and Chocolate are the two most positively correlated
variables, followed by Chocolate and Bar.

# 6. Principal Component Analysis

``` r
pca <- prcomp(candy, scale = TRUE )
summary(pca)
```

    Importance of components:
                              PC1    PC2    PC3     PC4    PC5     PC6     PC7
    Standard deviation     2.0788 1.1378 1.1092 1.07533 0.9518 0.81923 0.81530
    Proportion of Variance 0.3601 0.1079 0.1025 0.09636 0.0755 0.05593 0.05539
    Cumulative Proportion  0.3601 0.4680 0.5705 0.66688 0.7424 0.79830 0.85369
                               PC8     PC9    PC10    PC11    PC12
    Standard deviation     0.74530 0.67824 0.62349 0.43974 0.39760
    Proportion of Variance 0.04629 0.03833 0.03239 0.01611 0.01317
    Cumulative Proportion  0.89998 0.93832 0.97071 0.98683 1.00000

``` r
plot(pca$x[,1:2])
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
plot(pca$x[,1:2], col=my_cols, pch=16)
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-19-1.png)

``` r
# Make a new data-frame with our PCA results and candy data
my_data <- cbind(candy, pca$x[,1:3])
```

``` r
p <- ggplot(my_data) + 
        aes(x=PC1, y=PC2, 
            size=winpercent/100,  
            text=rownames(my_data),
            label=rownames(my_data)) +
        geom_point(col=my_cols)

p
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-21-1.png)

``` r
library(ggrepel)

p + geom_text_repel(size=3.3, col=my_cols, max.overlaps = 7)  + 
  theme(legend.position = "none") +
  labs(title="Halloween Candy PCA Space",
       subtitle="Colored by type: chocolate bar (dark brown), chocolate other (light brown), fruity (red), other (black)",
       caption="Data from 538")
```

    Warning: ggrepel: 39 unlabeled data points (too many overlaps). Consider
    increasing max.overlaps

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-22-1.png)

``` r
par(mar=c(8,4,2,2))
barplot(pca$rotation[,1], las=2, ylab="PC1 Contribution")
```

![](Halloween-Candy_files/figure-commonmark/unnamed-chunk-23-1.png)

- **Q24**. What original variables are picked up strongly by PC1 in the
  positive direction? Do these make sense to you?

  Fruity, Hard, and Pluribus are the 3 variables picked up strongly by
  PC1 in the positive direction. They make sense to me because it means
  these 3 variables are correlated just like Chocolate and Bar.
