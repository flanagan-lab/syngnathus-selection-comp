Selection pressures in *Syngnathus scovelli*
================



- <a href="#cleaning-the-datasets" id="toc-cleaning-the-datasets">Cleaning
  the datasets</a>
- <a href="#calculating-the-degree-of-sexual-dimorphism"
  id="toc-calculating-the-degree-of-sexual-dimorphism">Calculating the
  degree of sexual dimorphism</a>
  - <a href="#checking-the-assumptions-for-a-pairwise-comparison"
    id="toc-checking-the-assumptions-for-a-pairwise-comparison">Checking the
    assumptions for a pairwise comparison</a>
  - <a href="#investigate-distributions-and-run-the-tests"
    id="toc-investigate-distributions-and-run-the-tests">Investigate
    distributions and run the tests</a>
- <a href="#summary-statistics-for-successfully-mated-individuals"
  id="toc-summary-statistics-for-successfully-mated-individuals">Summary
  statistics for successfully mated individuals</a>
  - <a href="#males" id="toc-males">Males</a>
  - <a href="#females" id="toc-females">Females</a>
- <a href="#differences-between-mated-individuals-and-unmated-individuals"
  id="toc-differences-between-mated-individuals-and-unmated-individuals">Differences
  between mated individuals and unmated individuals</a>
  - <a href="#males-1" id="toc-males-1">Males</a>
    - <a href="#visual-comparison" id="toc-visual-comparison">Visual
      Comparison</a>
    - <a href="#testing-the-difference"
      id="toc-testing-the-difference">Testing the difference</a>
  - <a href="#females-1" id="toc-females-1">Females</a>
    - <a href="#visual-comparison-1" id="toc-visual-comparison-1">Visual
      Comparison</a>
    - <a href="#testing-the-difference-1"
      id="toc-testing-the-difference-1">Testing the difference</a>
- <a
  href="#looking-into-the-opportunity-for-selection-in-males-and-females"
  id="toc-looking-into-the-opportunity-for-selection-in-males-and-females">Looking
  into the Opportunity for Selection in Males and Females</a>
  - <a
    href="#generating-the-total-opportunity-for-selection-i-and-the-opportunity-for-sexual-selection-i_s"
    id="toc-generating-the-total-opportunity-for-selection-i-and-the-opportunity-for-sexual-selection-i_s">Generating
    the total opportunity for selection (<span
    class="math inline"><em>I</em></span>) and the opportunity for sexual
    selection (<span
    class="math inline"><em>I</em><sub><em>S</em></sub></span>)</a>
  - <a href="#partitioning-the-total-opportunity-for-selection-i"
    id="toc-partitioning-the-total-opportunity-for-selection-i">Partitioning
    the Total Opportunity for Selection (<span
    class="math inline"><em>I</em></span>)</a>
- <a
  href="#mate-success-versus-reproductive-success-bateman-gradient-beta_ss"
  id="toc-mate-success-versus-reproductive-success-bateman-gradient-beta_ss">Mate
  success versus Reproductive success (Bateman Gradient, <span
  class="math inline"><em>β</em><sub><em>S</em><em>S</em></sub></span>)</a>
  - <a href="#investigating-the-impact-of-zeros-on-the-bateman-gradient"
    id="toc-investigating-the-impact-of-zeros-on-the-bateman-gradient">Investigating
    the impact of “zeros” on the Bateman Gradient</a>
    - <a href="#removing-the-zeros-from-the-plot"
      id="toc-removing-the-zeros-from-the-plot">Removing the zeros from the
      plot</a>
    - <a href="#removing-the-zeros-from-the-calculation-of-relative-fitness"
      id="toc-removing-the-zeros-from-the-calculation-of-relative-fitness">Removing
      the zeros from the calculation of relative fitness</a>
- <a
  href="#investing-selection-differentials-on-snout-vent-length-s-and-s"
  id="toc-investing-selection-differentials-on-snout-vent-length-s-and-s">Investing
  selection differentials on snout-vent-length (<span
  class="math inline"><em>s</em></span> and <span
  class="math inline"><em>s</em>′</span>)</a>
  - <a href="#looking-into-the-maximum-sexual-selection-differential"
    id="toc-looking-into-the-maximum-sexual-selection-differential">Looking
    into the Maximum Sexual Selection Differential</a>
- <a href="#visualizing-post-copulatory-selection"
  id="toc-visualizing-post-copulatory-selection">Visualizing
  post-copulatory selection</a>

``` r
#This is a cohesive list of all the libraries used in this document
library(ggplot2)
library(cowplot)
library(fBasics)
library(pwr)
library(dplyr)
library(tidyr)
library(knitr)
```

``` r
#Metadata for males and females from the mesocosm experiments
fem_mesoSS <- read.csv("data/all_fem_meso_scovelli.csv")
mal_mesoSS <- read.csv("data/all_mal_meso_scovelli.csv")
```

This document will follow the same analysis as was outlined in
`selection_analysis_floridae.Rmd`. For more thorough details refer back
to that document.

The datasets used in this document were pulled from a publically
accessible manuscript that ran similar experimental breeding populations
(Rose et al., 2013).

# Cleaning the datasets

For the male dataset I want to do a few things before using it: 1.
Remove the two males that died (C1M5 and C6M2). 2. Replace all of the
“NAs” present in the males who didn’t mate with 0’s when appropriate. 3.
Add a column of mating success to the dataset. This will be either a 0
or 1 as *S. scovelli* males only mate once.

In the female dataset I am going to remove C1F2 since she doesn’t have
any data related to her reproductive success.

For both the male and the female datasets I also want to subset them out
to only include the control fish and not fish that were exposed to
estrogen.

``` r
#Subset the datasets to remove the fish exposed to estrogen
fem_succSS <- fem_mesoSS[grep("C", fem_mesoSS$trial_num),]
mal_succSS <- mal_mesoSS[grep("C", fem_mesoSS$trial_num),]

#Adding full fishIDs to make removing individuals easier
fem_succSS$femID <- paste0(fem_succSS$trial_num, "F",
                         fem_succSS$fishID)
mal_succSS$maleID <- paste0(mal_succSS$trial_num, "M",
                          mal_succSS$fishID)
#Removing the one female
fem_succSS <- subset(fem_succSS, !(femID %in% "C1F2"))

#Removing the two males who died
mal_succSS <- subset(mal_succSS, !(maleID %in% c("C1M5", "C6M2")))

#Replace NAs with 0s in the columns related to fitness
mal_succSS[, c("NumDeveloped", 
               "NumUndeveloped",
               "totalEggs")] <- sapply(mal_succSS[, c("NumDeveloped",
                                                      "NumUndeveloped",
                                                      "totalEggs")],
                                       function(x)
                                         as.numeric(ifelse(is.na(x), 0, x)))

#Add a column for females and males to denote mating success
mal_succSS$MatingSuccess <- ifelse(mal_succSS$totalEggs > 0, 
                                   1, 
                                   0)

fem_succSS$mated <- ifelse(fem_succSS$MatingSuccess > 0,
                           1,
                           0)
```

# Calculating the degree of sexual dimorphism

The original datasets contained only information about standard length.
To allow for a more complete comparison between all three species within
the genus *Syngnathus* I obtained the original photographs of those
pipefish and re-measured them to obtain measurements of both torso depth
and snout-vent length for males and females. For males I also generated
additional measurements of pouch area and pouch length.

Now, we can see if males and females differ in standard length, torso
depth, and snout-vent length. First, I need to see if assumptions are
met, i.e. variances are equal and data is normally distributed.

## Checking the assumptions for a pairwise comparison

The main two things that I will be looking into include:

1.  Equal variances between groups (using `var.test()`).
2.  Normal distribution of the data (using `normalTest()`).

To account for the fact that fish who are longer may just inherently be
deeper as well, I am going to adjust the depth by the snout-vent length
of the pipefish prior to running any analyses. This svl-adjusted depth
will then be used in all consequent analyses.

``` r
#Adjust the torso depth
fem_succSS$depth_adj <- fem_succSS$depth/fem_succSS$svl
mal_succSS$depth_adj <- mal_succSS$depth/mal_succSS$svl

#Testing to see if the variances are equal
var.test(fem_succSS$length, mal_succSS$length) #not equal
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_succSS$length and mal_succSS$length
    ## F = 0.47745, num df = 54, denom df = 53, p-value = 0.007693
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.2775688 0.8201320
    ## sample estimates:
    ## ratio of variances 
    ##          0.4774502

``` r
var.test(fem_succSS$depth_adj, mal_succSS$depth_adj) #not equal
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_succSS$depth_adj and mal_succSS$depth_adj
    ## F = 2.2275, num df = 54, denom df = 53, p-value = 0.00406
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.294973 3.826254
    ## sample estimates:
    ## ratio of variances 
    ##           2.227502

``` r
var.test(fem_succSS$svl, mal_succSS$svl) #equal
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_succSS$svl and mal_succSS$svl
    ## F = 0.60547, num df = 54, denom df = 53, p-value = 0.06888
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.3519922 1.0400308
    ## sample estimates:
    ## ratio of variances 
    ##          0.6054671

``` r
#Testing for normal distribution
normalTest(fem_succSS$length, method = "da") #not normal
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 8.6573
    ##     Z3  | Skewness: 2.4792
    ##     Z4  | Kurtosis: 1.5846
    ##   P VALUE:
    ##     Omnibus  Test: 0.01319 
    ##     Skewness Test: 0.01317 
    ##     Kurtosis Test: 0.1131

``` r
normalTest(mal_succSS$length, method = "da") #normal
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 1.1276
    ##     Z3  | Skewness: 0.4421
    ##     Z4  | Kurtosis: -0.9655
    ##   P VALUE:
    ##     Omnibus  Test: 0.569 
    ##     Skewness Test: 0.6584 
    ##     Kurtosis Test: 0.3343

``` r
normalTest(fem_succSS$depth_adj, method = "da") #normal
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 4.1601
    ##     Z3  | Skewness: -0.9031
    ##     Z4  | Kurtosis: -1.8288
    ##   P VALUE:
    ##     Omnibus  Test: 0.1249 
    ##     Skewness Test: 0.3665 
    ##     Kurtosis Test: 0.06743

``` r
normalTest(mal_succSS$depth_adj, method = "da") #normal
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 1.4518
    ##     Z3  | Skewness: -1.0588
    ##     Z4  | Kurtosis: 0.5752
    ##   P VALUE:
    ##     Omnibus  Test: 0.4839 
    ##     Skewness Test: 0.2897 
    ##     Kurtosis Test: 0.5652

``` r
normalTest(fem_succSS$svl, method = "da") #not normal
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 11.019
    ##     Z3  | Skewness: 2.8516
    ##     Z4  | Kurtosis: 1.6993
    ##   P VALUE:
    ##     Omnibus  Test: 0.004048 
    ##     Skewness Test: 0.00435 
    ##     Kurtosis Test: 0.08927

``` r
normalTest(mal_succSS$svl, method = "da") #normal
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 0.7681
    ##     Z3  | Skewness: 0.2322
    ##     Z4  | Kurtosis: -0.8451
    ##   P VALUE:
    ##     Omnibus  Test: 0.6811 
    ##     Skewness Test: 0.8164 
    ##     Kurtosis Test: 0.3981

## Investigate distributions and run the tests

I will run a Wilcoxon test for standard length and snout-vent length and
a Welch’s two-sample t-test for depth.

![*Histograms of male and female pipefish body
sizes.*](selection_analysis_scovelli_files/figure-gfm/histogram_sizes-1.png)

``` r
wilcox.test(fem_succSS$length, mal_succSS$length) #Sig. difference
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  fem_succSS$length and mal_succSS$length
    ## W = 2531, p-value = 2.353e-10
    ## alternative hypothesis: true location shift is not equal to 0

``` r
wilcox.test(fem_succSS$svl, mal_succSS$svl) #Sig. difference
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  fem_succSS$svl and mal_succSS$svl
    ## W = 2801, p-value = 1.552e-15
    ## alternative hypothesis: true location shift is not equal to 0

``` r
t.test(fem_succSS$depth_adj, mal_succSS$depth_adj, 
       var.equal = FALSE) #Sig. difference
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  fem_succSS$depth_adj and mal_succSS$depth_adj
    ## t = 14.983, df = 94.535, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.02186095 0.02853951
    ## sample estimates:
    ##  mean of x  mean of y 
    ## 0.11729070 0.09209047

For the Gulf pipefish, there are significant differences between males
and females in terms of standard length, snout-vent length, and torso
depth.

``` r
#Checking the power - length
d_mean_len <- abs(mean(fem_succSS$length, na.rm = TRUE) - 
                    mean(mal_succSS$length, na.rm = TRUE))
pool_sd_len <- sqrt((var(fem_succSS$length, na.rm = TRUE) + 
                       var(mal_succSS$length, na.rm = TRUE))/ 2)
d_len <- d_mean_len/pool_sd_len

pwr.t.test(n = nrow(fem_succSS), 
           d = d_len,
           sig.level = 0.05,
           type = 'two.sample',
           alternative = 'two.sided')
```

    ## 
    ##      Two-sample t test power calculation 
    ## 
    ##               n = 55
    ##               d = 1.492238
    ##       sig.level = 0.05
    ##           power = 1
    ##     alternative = two.sided
    ## 
    ## NOTE: n is number in *each* group

``` r
#Checking the power - SVL
d_mean_svl <- abs(mean(fem_succSS$svl, na.rm = TRUE) - 
                    mean(mal_succSS$svl, na.rm = TRUE))
pool_sd_svl <- sqrt((var(fem_succSS$svl, na.rm = TRUE) + 
                       var(mal_succSS$svl, na.rm = TRUE))/ 2)
d_svl <- d_mean_svl/pool_sd_svl

pwr.t.test(n = nrow(fem_succSS), 
           d = d_svl,
           sig.level = 0.05,
           type = 'two.sample',
           alternative = 'two.sided')
```

    ## 
    ##      Two-sample t test power calculation 
    ## 
    ##               n = 55
    ##               d = 2.16815
    ##       sig.level = 0.05
    ##           power = 1
    ##     alternative = two.sided
    ## 
    ## NOTE: n is number in *each* group

``` r
#Checking the power - Depth
d_mean_depth <- abs(mean(fem_succSS$depth_adj, na.rm = TRUE) - 
                      mean(mal_succSS$depth_adj, na.rm = TRUE))
pool_sd_depth <- sqrt((var(fem_succSS$depth_adj, na.rm = TRUE) + 
                         var(mal_succSS$depth_adj, na.rm = TRUE))/ 2)
d_depth <- d_mean_depth/pool_sd_depth
pwr.t.test(n = nrow(fem_succSS), 
           d = d_depth,
           sig.level = 0.05,
           type = 'two.sample',
           alternative = 'two.sided')
```

    ## 
    ##      Two-sample t test power calculation 
    ## 
    ##               n = 55
    ##               d = 2.865307
    ##       sig.level = 0.05
    ##           power = 1
    ##     alternative = two.sided
    ## 
    ## NOTE: n is number in *each* group

For all variables we have a power of over 0.9 or over 90% so we can be
confident in our interpretation.

# Summary statistics for successfully mated individuals

## Males

Across all 7 trials and 54 total males, there were 49 males that mated
at least one time and 0 of those males had two mates.

Looking across all males, including the ones that did not mate, this is
what we find as the mean, sd, and se for the number of embryos
transferred and how many of those developed versus didn’t:

|                     |       mean |         SD |        SE | max | min |
|:--------------------|-----------:|-----------:|----------:|----:|----:|
| Number of Embryos   | 25.8333333 | 12.6203639 |  1.717414 |  53 |   0 |
| Developed Embryos   | 23.1111111 | 13.3213783 |   1.81281 |  52 |   0 |
| Undeveloped Embryos |  2.7222222 |  4.7201562 | 0.6423319 |  22 |   0 |

These values will be influenced by the number of 0s coming from males
who did not mate. So let’s look at the same thing, but this time for
only males who had at least one successful mating:

|                     |       mean |         SD |        SE | max | min |
|:--------------------|-----------:|-----------:|----------:|----:|----:|
| Number of Embryos   | 28.4693878 |  9.9626428 | 1.4232347 |  53 |   7 |
| Developed Embryos   | 25.4693878 | 11.6029846 | 1.6575692 |  52 |   0 |
| Undeveloped Embryos |          3 |  4.8733972 | 0.6961996 |  22 |   0 |

I know want to see if there are correlations between morphometrics such
as standard length and the size of the brood pouch and brood size (i.e.,
are larger males securing larger broods).

![*Scatterplot of the relationship between brood pouch size metrics and
the number of embryos a male
had.*](selection_analysis_scovelli_files/figure-gfm/em-v-bp-1.png)

There may be some correlation happening here, but it doesn’t look
particularly strong. Let’s run some correlations tests to see what they
say.

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(mated_malSS$bp_area) and mated_malSS$totalEggs
    ## t = 1.228, df = 47, p-value = 0.2256
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1103484  0.4359036
    ## sample estimates:
    ##       cor 
    ## 0.1763193

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(mated_malSS$bp_length) and mated_malSS$totalEggs
    ## t = 1.7502, df = 47, p-value = 0.08661
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.03636651  0.49418250
    ## sample estimates:
    ##       cor 
    ## 0.2473596

There is no significant correlation with either of the brood pouch size
metrics. Let’s see if there is some relationship with just the overall
size of the male.

![*Scatterplot of the relationship between standard length (mm) and the
number of embryos a male
had.*](selection_analysis_scovelli_files/figure-gfm/em-v-sl-1.png)

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_malSS$length and mated_malSS$totalEggs
    ## t = 2.8197, df = 47, p-value = 0.007016
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.1110634 0.5976513
    ## sample estimates:
    ##       cor 
    ## 0.3803805

There is a significant correlation! Larger males have bigger broods in
terms of standard length.

## Females

Across all 7 trials and 55 total females, there were 29 females that
mated at least one time, 10 females that mated twice, and 4 that mated 3
times.

Looking across all females, including the ones that did not mate, this
is what we find as the mean, sd, and se for the total number of embryos
transferred from each female (across all of her mates if applicable) and
how many of those developed versus didn’t:

|                     |       mean |         SD |        SE | max | min |
|:--------------------|-----------:|-----------:|----------:|----:|----:|
| Number of Embryos   | 24.3818182 | 27.4364524 | 3.6995305 | 104 |   0 |
| Developed Embryos   | 22.5818182 | 25.5384176 | 3.4435995 |  95 |   0 |
| Undeveloped Embryos |        1.8 |   3.713838 | 0.5007738 |  14 |   0 |

These values will be influenced by the number of 0s coming from females
who did not mate. So let’s look at the same thing, but this time for
only females who had at least one successful mating:

|                     |       mean |         SD |        SE | max | min |
|:--------------------|-----------:|-----------:|----------:|----:|----:|
| Number of Embryos   | 46.2413793 | 20.1204501 | 3.7362738 | 104 |  18 |
| Developed Embryos   | 42.8275862 |  18.968148 | 3.5222967 |  95 |  15 |
| Undeveloped Embryos |  3.4137931 |  4.5710052 | 0.8488144 |  14 |   0 |

I want to see what relationship there may be between female body size
(in terms of standard length, depth, and SVL) and the number of eggs she
transferred. I also want to see on average how many eggs were
transferred per mating. I’m going to calculate this by taking the total
number of eggs and dividing it by the number of mates.

    ## [1] 30.11494

    ## [1] 1.293215

![Scatterplot of the relationship between female size metrics and the
number of eggs
transferred.](selection_analysis_scovelli_files/figure-gfm/em-v-fem-size-1.png)

There also appears to be a correlation between female body size and the
number of eggs transferred, especially in terms of depth and snout-vent
length. Let’s run some correlations tests to see what they say.

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_femSS$length and as.numeric(mated_femSS$totalEggs)
    ## t = 0.51174, df = 27, p-value = 0.613
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2784988  0.4484091
    ## sample estimates:
    ##        cor 
    ## 0.09801098

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_femSS$depth_adj and as.numeric(mated_femSS$totalEggs)
    ## t = 2.3662, df = 27, p-value = 0.02541
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.05649854 0.67795500
    ## sample estimates:
    ##       cor 
    ## 0.4144227

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_femSS$svl and as.numeric(mated_femSS$totalEggs)
    ## t = 1.4645, df = 27, p-value = 0.1546
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1057441  0.5801045
    ## sample estimates:
    ##      cor 
    ## 0.271275

There is no sig. correlation between length or svl and the number of
eggs transferred but we do see a significantly positive relationship
between depth and number of eggs transferred! It seems to be depth is a
good indicator of fecundity for females.

# Differences between mated individuals and unmated individuals

I want to now see if there are any significant differences in the sizes
of individuals who mated vs individuals that didn’t mate in males and
females. I am going to be focusing on the same morphometrics outlined
above.

## Males

### Visual Comparison

Before conducting any analyses, let’s see if we can visually detect any
differences between males who mated and unmated individuals.

![*Six different morphometrics compared between males who sucessfully
mated versus those that didn’t. Orange represents unmated and blue
represents mated
males.*](selection_analysis_scovelli_files/figure-gfm/mat-status-morph-mal-1.png)

I don’t notice many differences, however, it appears that males who
mated may be slightly larger.

### Testing the difference

Let’s now put some statistical power behind the difference in various
morphometrics between mated and unmated individuals. I am first going to
test the assumptions and then run the appropriate version of a t-test.

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  mal_succSS$length by as.factor(mal_succSS$MatingSuccess)
    ## F = 1.108, num df = 4, denom df = 48, p-value = 0.7274
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.3613617 9.2916889
    ## sample estimates:
    ## ratio of variances 
    ##           1.108019

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  mal_succSS$depth_adj by as.factor(mal_succSS$MatingSuccess)
    ## F = 1.0676, num df = 4, denom df = 48, p-value = 0.7658
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.3481692 8.9524694
    ## sample estimates:
    ## ratio of variances 
    ##           1.067568

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  mal_succSS$svl by as.factor(mal_succSS$MatingSuccess)
    ## F = 1.4997, num df = 4, denom df = 48, p-value = 0.4343
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##   0.4891137 12.5765727
    ## sample estimates:
    ## ratio of variances 
    ##           1.499737

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  mal_succSS$bp_area by as.factor(mal_succSS$MatingSuccess)
    ## F = 1.2012, num df = 4, denom df = 48, p-value = 0.6451
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##   0.3917452 10.0729378
    ## sample estimates:
    ## ratio of variances 
    ##           1.201182

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  mal_succSS$bp_length by as.factor(mal_succSS$MatingSuccess)
    ## F = 0.32826, num df = 4, denom df = 48, p-value = 0.2848
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.1070578 2.7527761
    ## sample estimates:
    ## ratio of variances 
    ##          0.3282642

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 1.1276
    ##     Z3  | Skewness: 0.4421
    ##     Z4  | Kurtosis: -0.9655
    ##   P VALUE:
    ##     Omnibus  Test: 0.569 
    ##     Skewness Test: 0.6584 
    ##     Kurtosis Test: 0.3343

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 1.4518
    ##     Z3  | Skewness: -1.0588
    ##     Z4  | Kurtosis: 0.5752
    ##   P VALUE:
    ##     Omnibus  Test: 0.4839 
    ##     Skewness Test: 0.2897 
    ##     Kurtosis Test: 0.5652

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 0.7681
    ##     Z3  | Skewness: 0.2322
    ##     Z4  | Kurtosis: -0.8451
    ##   P VALUE:
    ##     Omnibus  Test: 0.6811 
    ##     Skewness Test: 0.8164 
    ##     Kurtosis Test: 0.3981

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 2.3329
    ##     Z3  | Skewness: 1.3565
    ##     Z4  | Kurtosis: 0.702
    ##   P VALUE:
    ##     Omnibus  Test: 0.3115 
    ##     Skewness Test: 0.1749 
    ##     Kurtosis Test: 0.4827

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 6.2768
    ##     Z3  | Skewness: 2.1637
    ##     Z4  | Kurtosis: 1.263
    ##   P VALUE:
    ##     Omnibus  Test: 0.04335 
    ##     Skewness Test: 0.03048 
    ##     Kurtosis Test: 0.2066

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  mal_succSS$length by as.factor(mal_succSS$MatingSuccess)
    ## t = -1.3823, df = 52, p-value = 0.1728
    ## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
    ## 95 percent confidence interval:
    ##  -13.609115   2.507156
    ## sample estimates:
    ## mean in group 0 mean in group 1 
    ##        81.47800        87.02898

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  mal_succSS$depth_adj by as.factor(mal_succSS$MatingSuccess)
    ## t = -2.3573, df = 52, p-value = 0.02221
    ## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.013611257 -0.001093714
    ## sample estimates:
    ## mean in group 0 mean in group 1 
    ##      0.08541877      0.09277125

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  mal_succSS$svl by as.factor(mal_succSS$MatingSuccess)
    ## t = -0.41393, df = 52, p-value = 0.6806
    ## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
    ## 95 percent confidence interval:
    ##  -4.421717  2.909439
    ## sample estimates:
    ## mean in group 0 mean in group 1 
    ##        32.91580        33.67194

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  mal_succSS$bp_area by as.factor(mal_succSS$MatingSuccess)
    ## t = -0.61472, df = 52, p-value = 0.5414
    ## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
    ## 95 percent confidence interval:
    ##  -11.397150   6.051836
    ## sample estimates:
    ## mean in group 0 mean in group 1 
    ##        39.32220        41.99486

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  mal_succSS$bp_length by as.factor(mal_succSS$MatingSuccess)
    ## W = 100, p-value = 0.5115
    ## alternative hypothesis: true location shift is not equal to 0

There were no significant differences in terms of the sizes of males who
mated versus those who did not. This may also be a result of the small
sample size for unmated individuals.

Let’s explore this a bit more and overlay the distribution of all males
(mated and unmated) with the males who did mate and see how it varies
for snout-vent length and depth.

![*Overlay of the size range of males who mated on top of the size range
of all males for snout-vent length (left) and torso depth
(right).*](selection_analysis_scovelli_files/figure-gfm/mated-unmated-hist-1.png)

## Females

Similarly, now let’s see if we can identify any significant differences
in the morphometrics of females who were able to obtain mates versus
those who were unsuccessful.

### Visual Comparison

![*Four different morphometrics compared between females who sucessfully
mated versus those that didn’t. Orange represents unmated and blue
represents mated
females.*](selection_analysis_scovelli_files/figure-gfm/mat-status-morph-fem-1.png)

For all categories it appears that females are larger in length, depth,
and snout-vent length.

### Testing the difference

Let’s now put some statistical power behind the difference in various
morphometrics between mated and unmated individuals.

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_succSS$length by fem_succSS$mated
    ## F = 0.59253, num df = 25, denom df = 28, p-value = 0.1898
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.2741342 1.3030855
    ## sample estimates:
    ## ratio of variances 
    ##          0.5925279

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_succSS$depth_adj by fem_succSS$mated
    ## F = 0.55912, num df = 25, denom df = 28, p-value = 0.1458
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.2586767 1.2296090
    ## sample estimates:
    ## ratio of variances 
    ##          0.5591173

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_succSS$svl by fem_succSS$mated
    ## F = 0.72165, num df = 25, denom df = 28, p-value = 0.413
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.3338732 1.5870523
    ## sample estimates:
    ## ratio of variances 
    ##          0.7216508

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 8.6573
    ##     Z3  | Skewness: 2.4792
    ##     Z4  | Kurtosis: 1.5846
    ##   P VALUE:
    ##     Omnibus  Test: 0.01319 
    ##     Skewness Test: 0.01317 
    ##     Kurtosis Test: 0.1131

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 4.1601
    ##     Z3  | Skewness: -0.9031
    ##     Z4  | Kurtosis: -1.8288
    ##   P VALUE:
    ##     Omnibus  Test: 0.1249 
    ##     Skewness Test: 0.3665 
    ##     Kurtosis Test: 0.06743

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 11.019
    ##     Z3  | Skewness: 2.8516
    ##     Z4  | Kurtosis: 1.6993
    ##   P VALUE:
    ##     Omnibus  Test: 0.004048 
    ##     Skewness Test: 0.00435 
    ##     Kurtosis Test: 0.08927

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  fem_succSS$length by fem_succSS$mated
    ## W = 251, p-value = 0.03361
    ## alternative hypothesis: true location shift is not equal to 0

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  fem_succSS$depth_adj by fem_succSS$mated
    ## t = -2.8757, df = 53, p-value = 0.005794
    ## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.012789147 -0.002279209
    ## sample estimates:
    ## mean in group 0 mean in group 1 
    ##       0.1133181       0.1208523

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  fem_succSS$svl by fem_succSS$mated
    ## W = 237, p-value = 0.01783
    ## alternative hypothesis: true location shift is not equal to 0

Females who mated are significantly larger in terms of standard length,
torso depth, and snout-vent length.

Let’s explore this a bit more and overlay the distribution of sizes in
all females (mated and unmated) with the sizes of females who did mate
and see how it varies.

![*Overlay of the standard length, torso depth, and snout-vent length of
females who mated on top of the size range of all
females.*](selection_analysis_scovelli_files/figure-gfm/mated-unmated-hist-depth-1.png)

There is still quite a spread in the sizes of females who mated, but in
every case the female with the largest value secured a mate.

# Looking into the Opportunity for Selection in Males and Females

One of the benefits of using genetic parentage analysis is that we can
now calculate the opportunity for selection and the opportunity for
sexual selection in male and female pipefish.

## Generating the total opportunity for selection ($I$) and the opportunity for sexual selection ($I_S$)

Because each trial provides an independent “population” (i.e., pipefish
from one trial **cannot mate** with pipefish from another trial), I am
going to calculate these metrics for each trial individually and then I
will average it. With these I can then also generate 95% confidence
intervals which I will investigate for indications of significance in
two ways:

- If the confidence intervals **DO NOT** cross 0 -\> significant
  selection.
- If the confidence intervals between the sexes **DO NOT** cross -\>
  significantly different selection between the two sexes.

``` r
##FEMALES
#Create a dataframe to store the calculations of I and I_S in
fem_opp_selection <- data.frame(matrix(ncol = 3,
                                       nrow = 0))

colnames(fem_opp_selection) <- c("trial_num", "I", "I_s")

#Loop through the different trials and calculate I and I_S
for (trial in unique(fem_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- fem_succSS[fem_succSS$trial_num == trial, ]
  
  #Calculate opportunity selection
  I <- var(tmp$NumDeveloped)/(mean(tmp$NumDeveloped)^2)
  
  I_s <- var(tmp$MatingSuccess)/(mean(tmp$MatingSuccess)^2)
  
  #Combining all of the selection values (Is) and save the output
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  selection <- cbind(trial_num, I, I_s)
  
  fem_opp_selection <- rbind(fem_opp_selection, selection)
  
}


##MALES
#Create a dataframe to store the calculations of I and I_S in
mal_opp_selection <- data.frame(matrix(ncol = 3,
                                       nrow = 0))

colnames(mal_opp_selection) <- c("trial_num", "I", "I_s")

#Loop through the different trials and calculate I and I_S
for (trial in unique(mal_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- mal_succSS[mal_succSS$trial_num == trial, ]
  
  #Calculate opportunity selection
  I <- var(tmp$NumDeveloped)/(mean(tmp$NumDeveloped)^2)
  
  I_s <- var(tmp$MatingSuccess)/(mean(tmp$MatingSuccess)^2)
  
  #Combining all of the selection values (Is) and save the output
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  selection <- cbind(trial_num, I, I_s)
  
  mal_opp_selection <- rbind(mal_opp_selection, selection)
  
}

#Merge the selection coefficients from males and females into one dataset to 
#make life easier
fem_opp_selection$Sex <- "F"
mal_opp_selection$Sex <- "M"

opp_selection_all <- rbind(fem_opp_selection, mal_opp_selection)
```

Now that I have calculated the opportunity for selection and sexual
selection, I want to generate my averages and 95% CI for both.

``` r
#List the columns of interest
columns <- c("I", "I_s")

#Create a dataframe to store the final values in
opp_average <- data.frame(matrix(ncol = 4,
                                 nrow = 0))
colnames(opp_average) <- c("Average", "Interval", "Episode_sel", "Sex")

#Calculate the critical value
crit <- qt(p = 0.975, df = (nrow(fem_opp_selection) - 1))

for (j in 1:length(columns)) {
    
    col_name <- columns[[j]]
    
    #Calculate the means
    mean <- t(t(tapply(opp_selection_all[, colnames(opp_selection_all) == col_name], 
                       opp_selection_all$Sex, 
                       mean)))
    
    #Calculate standard error
    se <- t(t(tapply(opp_selection_all[, colnames(opp_selection_all) == col_name], 
                     opp_selection_all$Sex, 
                 function(x){
                   sqrt(var(x))/sqrt(length(x))
                 })))
    
    #Calculate the value that is added and subtracted from the mean
    int <- se*crit
    
    #Combine the data together
    episode <- as.data.frame(cbind(mean, int))
    colnames(episode) <- c("Average", "Interval")
    
    episode$Episode_sel <- col_name
    episode$Sex <- rownames(episode)
    
    rownames(episode) <- NULL
    
    opp_average <- rbind(opp_average, episode)
    
  }
```

Let’s now explore some results:

| Episode_sel | F                 | M                  |
|:------------|:------------------|:-------------------|
| I           | 1.64 (0.53, 2.76) | 0.35 (0.11, 0.58)  |
| I_s         | 1.58 (0.7, 2.46)  | 0.15 (-0.08, 0.38) |

Average Opportunity of Selection (95% CI) for Males and Females

![*Average opportunity for selection and opportunity for sexual
selection for male (purple) and female (green) S. scovelli. Errorbars
represent the 95% confidence intervals around the
mean*](selection_analysis_scovelli_files/figure-gfm/opp-selection-figure-1.png)

We can see that for female *S. scovelli* there is a significant
opportunity for selection and opportunity for sexual selection, however,
males only experience a significant opportunity for selection, but not
sexual selection. Additionally, we can see a significant difference in
the opportunity for sexual selection between males and females.

## Partitioning the Total Opportunity for Selection ($I$)

Once again for partitioning the opportunity for selection, I am going to
calculate selection for the trials individually in males and females and
then average across all trials to get the final values and the 95% CIs.
For pre-mating processes I am focusing on mating success and then for
post-mating processes I am looking at the total number of eggs
transferred/received and the proportion of those eggs which developed
(showing fertilization success).

``` r
#Create a dataframe to store all of the intermediate values of fitness in
fem_succ_fitness <- data.frame(matrix(ncol = ncol(fem_succSS) + 9,
                                      nrow = 0))
colnames(fem_succ_fitness) <- c(colnames(fem_succSS),
                                "w1", "w1_squared",
                                "W2", "W2_bar", "w2",
                                "W3", "W3_bar", "w3", "i3")

#Create a dataframe to store the final calculations of I in
opp_selection_episodes_fem <- data.frame(matrix(ncol = 12,
                                            nrow = 0))
colnames(opp_selection_episodes_fem) <- c("trial_num", "I_1", "I_1per", "I_2", "I_2per", 
                                          "I_3", "I_3per", "I_12", "I_12per",
                                          "I", "Iper")

for (trial in unique(fem_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- fem_succSS[fem_succSS$trial_num == trial, ]
  
  #Calculate the absolute pre-copulatory fitness (Eq. 14 Arnold & Wade 1984)
  #This is the same as the calculation of I_s
  tmp$w1 <- tmp$MatingSuccess/mean(tmp$MatingSuccess) #Relative mating success
  tmp$w1_squared <- (tmp$w1)^2
  
  I_1 <- var(tmp$w1) #Variance in relative mating success
  
  #Post-copulatory selection event 1 (Number of eggs transferred) (Eq. 15 Arnold & Wade 1984)
  tmp$W2 <- ifelse(tmp$MatingSuccess > 0,
                   tmp$totalEggs/tmp$MatingSuccess,
                   0) #Number of eggs per mate
  tmp$W2_bar <- tmp$W2 * (tmp$w1/nrow(tmp)) #Number of eggs per mate adjusted by the # of individuals with fitness W
  tmp$w2 <- tmp$W2/sum(tmp$W2_bar)
  
  I_2 <- (sum((tmp$w1 * (tmp$w2)^2))/nrow(tmp) - 1) * nrow(tmp)/(nrow(tmp) - 1)
  
  #Post-copulatory selection event 2 (Number of eggs developed) (Eq. 16 Arnold & Wade 1984)
  tmp$W3 <- ifelse(tmp$totalEggs > 0,
                   tmp$NumDeveloped/tmp$totalEggs,
                   0) #Proportion of transferred eggs that developed
  tmp$W3_bar <- tmp$W3 * ((tmp$totalEggs/mean(tmp$totalEggs))/nrow(tmp)) #Prop. of eggs developed adjusted by the # of individuals with fitness W
  tmp$w3 <- tmp$W3/sum(tmp$W3_bar)
  tmp$i3 <- ((tmp$totalEggs/mean(tmp$totalEggs))/nrow(tmp)) * ((tmp$w3 - 1)^2)
  
  I_3 <- sum(tmp$i3) * nrow(tmp)/(nrow(tmp) - 1)

  I_12 <- var(tmp$totalEggs)/(mean(tmp$totalEggs)^2)
  
  #Total opportunity for selection
  I <- var(tmp$NumDeveloped)/(mean(tmp$NumDeveloped)^2)
  
  #Calculating percentages for each selection event
  I_1per <- (I_1/I)*100
  I_2per <- (I_2/I)*100
  I_3per <- (I_3/I)*100
  I_12per <- (I_12/I)*100
  Iper <- (I/I)*100
  
  #Combining all of the selection values (Is) and saving the output
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  selection <- cbind(trial_num, I_1, I_1per, I_2, I_2per, I_3, I_3per,
                     I_12, I_12per, I, Iper)
  
  opp_selection_episodes_fem <- rbind(opp_selection_episodes_fem, selection)
  
  #Save the intermediate values
  fem_succ_fitness <- rbind(fem_succ_fitness, tmp)
}

#Exporting the data
#write.csv(fem_succ_fitness, "data/scovelli_int_I_fem.csv", row.names = FALSE)
```

``` r
#Create a dataframe to store all of the intermediate values of fitness in
mal_succ_fitness <- data.frame(matrix(ncol = ncol(mal_succSS) + 9,
                                      nrow = 0))
colnames(mal_succ_fitness) <- c(colnames(mal_succSS),
                                "w1", "w1_squared",
                                "W2", "W2_bar", "w2",
                                "W3", "W3_bar", "w3", "i3")

#Create a dataframe to store the final calculations of I in
opp_selection_episodes_mal <- data.frame(matrix(ncol = 12,
                                            nrow = 0))
colnames(opp_selection_episodes_mal) <- c("trial_num", "I_1", "I_1per", "I_2", "I_2per", 
                                          "I_3", "I_3per", "I_12", "I_12per",
                                          "I", "Iper", "I_s")

for (trial in unique(mal_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- mal_succSS[mal_succSS$trial_num == trial, ]
  
  #Calculate the absolute pre-copultory fitness (Eq. 14 Arnold & Wade 1984)
  tmp$w1 <- tmp$MatingSuccess/mean(tmp$MatingSuccess) #Relative mating success
  tmp$w1_squared <- (tmp$w1)^2
  
  I_1 <- var(tmp$w1) #Variance in relative mating success
  
  #Post-copulatory selection event 1 (Number of eggs transferred) (Eq. 15 Arnold & Wade 1984)
  tmp$W2 <- ifelse(tmp$MatingSuccess > 0,
                   tmp$totalEggs/tmp$MatingSuccess,
                   0) #Number of eggs per mate
  tmp$W2_bar <- tmp$W2 * (tmp$w1/nrow(tmp)) #Number of eggs per mate adjusted by the # of individuals with fitness W
  tmp$w2 <- tmp$W2/sum(tmp$W2_bar)
  
  I_2 <- (sum((tmp$w1 * (tmp$w2)^2))/nrow(tmp) - 1) * nrow(tmp)/(nrow(tmp) - 1)
  
  #Post-copulatory selection event 2 (Number of eggs developed) (Eq. 16 Arnold & Wade 1984)
  tmp$W3 <- ifelse(tmp$totalEggs > 0,
                   tmp$NumDeveloped/tmp$totalEggs,
                   0) #Proportion of transferred eggs that developed
  tmp$W3_bar <- tmp$W3 * ((tmp$totalEggs/mean(tmp$totalEggs))/nrow(tmp)) #Prop. of eggs developed adjusted by the # of individuals with fitness W
  tmp$w3 <- tmp$W3/sum(tmp$W3_bar)
  tmp$i3 <- ((tmp$totalEggs/mean(tmp$totalEggs))/nrow(tmp)) * ((tmp$w3 - 1)^2)
  
  I_3 <- sum(tmp$i3) * nrow(tmp)/(nrow(tmp) - 1)

  I_12 <- var(tmp$totalEggs)/(mean(tmp$totalEggs)^2)
  
  #Total opportunity for selection
  I <- var(tmp$NumDeveloped)/(mean(tmp$NumDeveloped)^2)

  #Calculating percentages for each selection event
  I_1per <- (I_1/I)*100
  I_2per <- (I_2/I)*100
  I_3per <- (I_3/I)*100
  I_12per <- (I_12/I)*100
  Iper <- (I/I)*100
  
  #Combining all of the selection values (Is) and saving the output
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  selection <- cbind(trial_num, I_1, I_1per, I_2, I_2per, I_3, I_3per,
                     I_12, I_12per, I, Iper)
  
  opp_selection_episodes_mal <- rbind(opp_selection_episodes_mal, selection)
  
  #Save the intermediate values
  mal_succ_fitness <- rbind(mal_succ_fitness, tmp)
}

#Exporting the data
#write.csv(mal_succ_fitness, "data/scovelli_int_I_mal.csv", row.names = FALSE)
```

``` r
#Merge the selection coefficients from males and females into one dataset to 
#make life easier
opp_selection_episodes_fem$Sex <- "F"
opp_selection_episodes_mal$Sex <- "M"

opp_selection_episodes_all <- rbind(opp_selection_episodes_fem, opp_selection_episodes_mal)

#Exporting the data
#write.csv(opp_selection_episodes_all, "data/scovelli_opp_selection.csv", row.names = FALSE)

#List the columns of interest
columns <- c("I_1", "I_2", "I_12", "I_3","I")

#Create a dataframe to store the final values in
opp_episodes_average <- data.frame(matrix(ncol = 4,
                                    nrow = 0))
colnames(opp_episodes_average) <- c("Average", "Interval", 
                                    "Episode_sel", "Sex")

#Calculate the critical value
crit <- qt(p = 0.975, df = (nrow(opp_selection_episodes_fem) - 1))

for (j in 1:length(columns)) {
    
    col_name <- columns[[j]]
    
    #Calculate the means
    mean <- t(t(tapply(opp_selection_episodes_all[, colnames(opp_selection_episodes_all) 
                                                  == col_name], 
                       opp_selection_episodes_all$Sex, 
                       mean)))
    
    #Calculate standard error
    se <- t(t(tapply(opp_selection_episodes_all[, colnames(opp_selection_episodes_all) 
                                                == col_name], 
                     opp_selection_episodes_all$Sex, 
                 function(x){
                   sqrt(var(x))/sqrt(length(x))
                 })))
    
    #Calculate the value that is added and subtracted from the mean
    int <- se*crit
    
    #Combine the data together
    episode <- as.data.frame(cbind(mean, int))
    colnames(episode) <- c("Average", "Interval")
    
    episode$Episode_sel <- col_name
    episode$Sex <- rownames(episode)
    
    rownames(episode) <- NULL
    
    opp_episodes_average <- rbind(opp_episodes_average, episode)
    
  }
```

Let’s now explore some results:

| Episode_sel | F                  | M                   |
|:------------|:-------------------|:--------------------|
| I_1         | 1.582 (0.7, 2.46)  | 0.149 (-0.08, 0.38) |
| I_2         | 0.065 (0.01, 0.12) | 0.098 (0.04, 0.16)  |
| I_12        | 1.574 (0.62, 2.53) | 0.256 (0.03, 0.49)  |
| I_3         | 0.009 (0, 0.02)    | 0.06 (-0.01, 0.13)  |
| I           | 1.641 (0.53, 2.76) | 0.345 (0.11, 0.58)  |

Average Episode of Selection (95% CI) for Males and Females

![*Average opportunity for selection for the different episodes for male
(purple) and female (green) S. scovelli. Errorbars represent the 95%
confidence intervals around the
mean*](selection_analysis_scovelli_files/figure-gfm/opp-select-episodes-figure-1.png)

From the table and the plot we can see that once again there are some
significant differences in the selection between males and females.

Let’s now look more into the percentage of the overall opportunity for
selection made up for by each individual episode of selection:

``` r
sexes <- c("M", "F")

#Create a dataframe to store the final values in
opp_percents <- data.frame(matrix(ncol = 3,
                                 nrow = 0))
colnames(opp_percents) <- c("Percent", "Episode_sel", "Sex")

#Calculate the percentage for each episode in the two sexes
for (sex in sexes) {
    
  #subset dataset based on sex
  tmp_sex <- opp_episodes_average[opp_episodes_average$Sex == sex, ]
  
  #Pull out the overall opp. for selection value
  I_average <- tmp_sex$Average[tmp_sex$Episode_sel == "I"]
  
  #Calculate what percentage of the I is represented by each episode
  percents <- as.data.frame(t(t(apply(tmp_sex, 1, function(x){
      
      (as.numeric(x[1])/I_average)*100
      
    }))))
    
    colnames(percents) <- "Percent"
    percents$Episode_sel <- tmp_sex$Episode_sel
    percents$Sex <- sex
    
    rownames(percents) <- NULL
    
    opp_percents <- rbind(opp_percents, percents)
    
  }
```

![*The proportion of the total opportunity for selection that is
represented by each episode of selection for males and
females.*](selection_analysis_scovelli_files/figure-gfm/generate-figa-opp-selection-1.png)

Matching the previous plots, most of the opportunity for selection in
*Syngnathus scovelli* females can be attributed to variance in mating
success ($I_1$) rather than variance in eggs transferred/received
($I_2$) or variance in the proportion of eggs developed ($I_3$).
However, for males it appears that post-copulatory selection may be more
important. This may be a result of males only having the capacity to
mate once, meaning there is less variation in the mating success (i.e,
they either mate or they don’t).

# Mate success versus Reproductive success (Bateman Gradient, $\beta_{SS}$)

To calculate $\beta_{SS}$ we use *relative* measures of fitness:
($\frac{indvidual's fitness}{mean(fitness)}$)

I am going to generate the measurements of relative fitness once again
for the trials individually rather than across all of the trials.

``` r
#Calculating relative fitness as a metric for reproductive success
#Create a dataframe to store all of the calculations of relative fitness in
fem_bateman <- data.frame(matrix(ncol = 3,
                                 nrow = 0))
colnames(fem_bateman) <- c("trial", "MatingSuccess","rel_repo_fitness")

#Loop through each trial to calculate relative fitness
for (trial in unique(fem_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- fem_succSS[fem_succSS$trial_num == trial, ]
  
  #Calculate relative fitness
  rel_repo_fitness <- tmp$totalEggs/mean(tmp$totalEggs)
  
  #Calculte mating fitness
  rel_mate_succuess <- tmp$MatingSuccess/mean(tmp$MatingSuccess)
  
  #Column-bind the trial #, Mating success, and calculated rel. fitness
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  fitness <- cbind("trial" = rep(trial_num, nrow(tmp)), 
                   "MatingSuccess" = rel_mate_succuess, 
                   rel_repo_fitness)
  
  #Add this chunk of data to the dataframe we created
  fem_bateman <- rbind(fem_bateman, fitness)
}


#Repeat process for the Male mating data
mal_bateman <- data.frame(matrix(ncol = 3,
                                 nrow = 0))
colnames(mal_bateman) <- c("trial", "MatingSuccess","rel_repo_fitness")

for (trial in unique(mal_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- mal_succSS[mal_succSS$trial_num == trial, ]
  
  #Calculate relative fitness
  rel_repo_fitness <- tmp$totalEggs/mean(tmp$totalEggs)
  
  #Calculate mating fitness
  rel_mate_succuess <- tmp$MatingSuccess/mean(tmp$MatingSuccess)
  
  #Column-bind the trial #, Mating success, and calculated rel. fitness
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  fitness <- cbind("trial" = rep(trial_num, nrow(tmp)), 
                   "MatingSuccess" = rel_mate_succuess, 
                   rel_repo_fitness)
  
  #Add this chunk of data to the dataframe we created
  mal_bateman <- rbind(mal_bateman, fitness)
}
```

Once we have the measures of relative fitness we can use them to run the
weighted least-squares regression for males and females separately.

``` r
#Generating Bateman's gradient
#Define the model
fem_model <- lm(fem_bateman$rel_repo_fitness ~ fem_bateman$MatingSuccess)
mal_model <- lm(mal_bateman$rel_repo_fitness ~ mal_bateman$MatingSuccess)

#define weights to use
wt_fem <- 1 / lm(abs(fem_model$residuals) ~ fem_model$fitted.values)$fitted.values^2
wt_mal <- 1 / lm(abs(mal_model$residuals) ~ mal_model$fitted.values)$fitted.values^2

#perform weighted least squares regression
wls_model_fem <- lm(fem_bateman$rel_repo_fitness ~ fem_bateman$MatingSuccess,
                    weights=wt_fem)
wls_model_mal <- lm(mal_bateman$rel_repo_fitness ~ mal_bateman$MatingSuccess,
                    weights=wt_mal)

#Investigate the results
summary(wls_model_fem) #significant
```

    ## 
    ## Call:
    ## lm(formula = fem_bateman$rel_repo_fitness ~ fem_bateman$MatingSuccess, 
    ##     weights = wt_fem)
    ## 
    ## Weighted Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2071 -0.2194 -0.0661  0.0434  3.2156 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               0.005085   0.016835   0.302    0.764    
    ## fem_bateman$MatingSuccess 1.012927   0.036178  27.999   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.126 on 53 degrees of freedom
    ## Multiple R-squared:  0.9367, Adjusted R-squared:  0.9355 
    ## F-statistic: 783.9 on 1 and 53 DF,  p-value: < 2.2e-16

``` r
summary(wls_model_mal) #significant
```

    ## 
    ## Call:
    ## lm(formula = mal_bateman$rel_repo_fitness ~ mal_bateman$MatingSuccess, 
    ##     weights = wt_mal)
    ## 
    ## Weighted Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1059 -0.9700  0.0000  0.6263  2.9816 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -3.230e-16  3.560e-02    0.00        1    
    ## mal_bateman$MatingSuccess  1.000e+00  5.242e-02   19.08   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.224 on 52 degrees of freedom
    ## Multiple R-squared:  0.875,  Adjusted R-squared:  0.8726 
    ## F-statistic:   364 on 1 and 52 DF,  p-value: < 2.2e-16

For both males and females there is a significant slope, meaning both
sexes see increases in fitness with each additional mating. I am
interested to run the model with the two datasets combined to see if
there is an interaction of sex.

``` r
#Combine the two datasets
fem_bateman$Sex <- "F"
mal_bateman$Sex <- "M"

all_bateman <- rbind(fem_bateman, mal_bateman)

#Running a weighted least squares regression between MS and Sex
MS_sex_model <- lm(all_bateman$rel_repo_fitness ~
                     all_bateman$MatingSuccess*all_bateman$Sex)

wt_all <- 1 / lm(abs(MS_sex_model$residuals) ~
                   MS_sex_model$fitted.values)$fitted.values^2

wls_MS_sex_model <- lm(all_bateman$rel_repo_fitness ~ 
                         all_bateman$MatingSuccess*all_bateman$Sex, weights = wt_all)

summary(wls_MS_sex_model)
```

    ## 
    ## Call:
    ## lm(formula = all_bateman$rel_repo_fitness ~ all_bateman$MatingSuccess * 
    ##     all_bateman$Sex, weights = wt_all)
    ## 
    ## Weighted Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4732 -0.5445 -0.0676  0.5699  3.3343 
    ## 
    ## Coefficients:
    ##                                             Estimate Std. Error t value
    ## (Intercept)                                 0.006342   0.021924   0.289
    ## all_bateman$MatingSuccess                   1.009830   0.042332  23.855
    ## all_bateman$SexM                           -0.006342   0.051927  -0.122
    ## all_bateman$MatingSuccess:all_bateman$SexM -0.009830   0.070392  -0.140
    ##                                            Pr(>|t|)    
    ## (Intercept)                                   0.773    
    ## all_bateman$MatingSuccess                    <2e-16 ***
    ## all_bateman$SexM                              0.903    
    ## all_bateman$MatingSuccess:all_bateman$SexM    0.889    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.205 on 105 degrees of freedom
    ## Multiple R-squared:  0.9093, Adjusted R-squared:  0.9067 
    ## F-statistic: 350.8 on 3 and 105 DF,  p-value: < 2.2e-16

When we combine all of the data together we can see two additional
things:

1.  There is no significant interactions effect of Mating success and
    sex
2.  The slopes for males and females are not significantly different
    from each other

Let’s visually look at this pattern now:

![*Relationship between reproductive success and mating success for male
(purple) and female (green) *Syngnathus scovelli*. Reproductive success
is shown as relative fitness (i.e. number of offspring produced divided
by the mean number of offspring produced). Bateman’s gradient is shown
as the weighted least-squares regression line (dashed) for males and
females.*](selection_analysis_scovelli_files/figure-gfm/plot-bateman3-1.png)

The plot confirms the results from the two models. We can see that there
is a steep slope for both males and females, however, the two lines are
not different.

## Investigating the impact of “zeros” on the Bateman Gradient

It has been shown previously that the inclusion or exclusion of
individuals who were unable to achieve a mate. When the non-mated
individuals are included, the gradient is influenced by the relative
fitness gain that is experienced from gaining a mate (i.e, moving from 0
to 1). If we do not include the non-mated individuals, we can more
clearly look at the fitness increase associated with obtaining more than
one mate.

For *S. scovelli*, because males only mate once it does not make sense
to re-do the analysis without the non-mated individuals. Therefore, I
will only be focusing on the females.

### Removing the zeros from the plot

The first way I am addressing this is by using the same datasets as
before and just excluding the 0’s from the plot and the model:

``` r
#Generating Bateman's gradient
#Define the model
fem_model2 <- lm(fem_bateman$rel_repo_fitness[fem_bateman$MatingSuccess != 0] ~
                   fem_bateman$MatingSuccess[fem_bateman$MatingSuccess != 0])

#define weights to use
wt_fem2 <- 1 / lm(abs(fem_model2$residuals) ~
                    fem_model2$fitted.values)$fitted.values^2

#perform weighted least squares regression
wls_model_fem2 <- lm(fem_bateman$rel_repo_fitness[fem_bateman$MatingSuccess != 0] ~
                       fem_bateman$MatingSuccess[fem_bateman$MatingSuccess != 0],
                    weights=wt_fem2)

#Investigate the results
summary(wls_model_fem2) #significant
```

    ## 
    ## Call:
    ## lm(formula = fem_bateman$rel_repo_fitness[fem_bateman$MatingSuccess != 
    ##     0] ~ fem_bateman$MatingSuccess[fem_bateman$MatingSuccess != 
    ##     0], weights = wt_fem2)
    ## 
    ## Weighted Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8981 -0.4245 -0.1745  0.8043  3.4571 
    ## 
    ## Coefficients:
    ##                                                           Estimate Std. Error
    ## (Intercept)                                                 0.2879     0.1715
    ## fem_bateman$MatingSuccess[fem_bateman$MatingSuccess != 0]   0.8465     0.1005
    ##                                                           t value Pr(>|t|)    
    ## (Intercept)                                                 1.679    0.105    
    ## fem_bateman$MatingSuccess[fem_bateman$MatingSuccess != 0]   8.426 4.88e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.361 on 27 degrees of freedom
    ## Multiple R-squared:  0.7245, Adjusted R-squared:  0.7143 
    ## F-statistic:    71 on 1 and 27 DF,  p-value: 4.884e-09

![*Relationship between reproductive success and mating success for
female *Syngnathus scovelli* who achieved at least one mate.
Reproductive success is shown as relative fitness (i.e. number of
offspring produced divided by the mean number of offspring produced).
Bateman’s gradient is shown as the weighted least-squares regression
line
(dashed).*](selection_analysis_scovelli_files/figure-gfm/plot-bateman2-1.png)

When we exclude the non-mated individuals from the plot and the model,
we can see the results do not change for females, the slope is still
significant.

### Removing the zeros from the calculation of relative fitness

The other way to approach this is rather than just eliminating the zeros
from the scatter plot, I can remove the zeros in the calculation of
relative fitness.

``` r
#Calculating relative fitness as a metric for reproductive success
#Create a dataframe to store all of the calculations of relative fitness in
fem_bateman_nozero <- data.frame(matrix(ncol = 3,
                                 nrow = 0))
colnames(fem_bateman_nozero) <- c("trial", "MatingSuccess","rel_repo_fitness")

#Loop through each trial to calculate relative fitness
for (trial in unique(fem_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- fem_succSS[fem_succSS$trial_num == trial, ]
  
  #Calculate relative fitness
  rel_repo_fitness <- tmp$totalEggs[tmp$MatingSuccess != 0]/
    mean(tmp$totalEggs[tmp$MatingSuccess != 0])
  
  #Calculte mating fitness
  rel_mate_succuess <- tmp$MatingSuccess[tmp$MatingSuccess != 0]/
    mean(tmp$MatingSuccess[tmp$MatingSuccess != 0])
  
  #Column-bind the trial #, Mating success, and calculated rel. fitness
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  fitness <- cbind("trial" = rep(trial_num, nrow(tmp[tmp$MatingSuccess != 0,])), 
                   "MatingSuccess" = rel_mate_succuess, 
                   rel_repo_fitness)
  
  #Add this chunk of data to the dataframe we created
  fem_bateman_nozero <- rbind(fem_bateman_nozero, fitness)
}
```

Once we have the measures of relative fitness we can use them to run the
weighted least-squares regression for males and females separately.

``` r
#Generating Bateman's gradient
#Define the model
fem_model3 <- lm(fem_bateman_nozero$rel_repo_fitness ~ 
                  fem_bateman_nozero$MatingSuccess)

#define weights to use
wt_fem3 <- 1 / lm(abs(fem_model3$residuals) ~
                    fem_model3$fitted.values)$fitted.values^2

#perform weighted least squares regression
wls_model_fem3 <- lm(fem_bateman_nozero$rel_repo_fitness ~
                      fem_bateman_nozero$MatingSuccess,
                    weights=wt_fem3)


#Investigate the results
summary(wls_model_fem3) #significant
```

    ## 
    ## Call:
    ## lm(formula = fem_bateman_nozero$rel_repo_fitness ~ fem_bateman_nozero$MatingSuccess, 
    ##     weights = wt_fem3)
    ## 
    ## Weighted Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9390 -0.7931 -0.1342  0.7556  3.5060 
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        0.2202     0.1135   1.941   0.0628 .  
    ## fem_bateman_nozero$MatingSuccess   0.7797     0.1092   7.140 1.12e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.359 on 27 degrees of freedom
    ## Multiple R-squared:  0.6537, Adjusted R-squared:  0.6409 
    ## F-statistic: 50.97 on 1 and 27 DF,  p-value: 1.12e-07

![*Relationship between reproductive success and mating success for
female *Syngnathus scovelli* who achieved at least one mate.
Reproductive success is shown as relative fitness (i.e. number of
offspring produced divided by the mean number of offspring produced).
Relative fitness is calculated without the individuals who did not mate.
Bateman’s gradient is shown as the weighted least-squares regression
line
(dashed).*](selection_analysis_scovelli_files/figure-gfm/plot-bateman-1.png)

With this way of excluding the individuals who did not mate, there is
still a significant increase in relative fitness with each additional
mating.

# Investing selection differentials on snout-vent-length ($s$ and $s'$)

A selection differential is the covariance between a trait and relative
fitness and is often thought of as the difference in a trait mean before
and after selection.

Just as I did for the decomposition of the opportunity for selection, I
am going to decompose the total selection differentials into the
different pre- and post-mating episodes.

To keep things consistent with the analysis I conducted on *S.
floridae*, I will once again be investigating snout-vent length as the
trait of interest, even though standard length and torso depth was also
significantly different between males and females.

``` r
#Create a dataframe to store all of the intermediate values of fitness in
fem_succ_select_diff <- data.frame(matrix(ncol = ncol(fem_succSS) + 6,
                                          nrow = 0))
colnames(fem_succ_select_diff) <- c(colnames(fem_succSS),
                                    "fit1", "eggs_per_mate","fit2", 
                                    "prop_dev", "fit3", "StdLength")

#Create a dataframe to store the final calculations of I in
select_diff_fem <- data.frame(matrix(ncol = 11,
                                     nrow = 0))
colnames(select_diff_fem) <- c("trial", "s1", "s2", "s3", "s12", "s123",
                               "s1_prime", "s2_prime", "s3_prime", 
                               "s12_prime", "s123_prime")

for (trial in unique(fem_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- fem_succSS[fem_succSS$trial_num == trial, ]
  
  #Calculate fitness relating to pre-cop. selection (#matings)
  tmp$fit1 <- tmp$MatingSuccess/mean(tmp$MatingSuccess) #Relative mating success

  #Calculate fitness relating to post-mating selection (#eggs transferred)
  tmp$eggs_per_mate <- tmp$totalEggs/tmp$MatingSuccess
  ##If mating success = 0, eggs_per_mate = NA and it not included in the calculation
  ##of the relative fitness moving forward
  tmp$fit2 <- ifelse(tmp$MatingSuccess > 0,
                     tmp$eggs_per_mate/mean(tmp$eggs_per_mate, na.rm = TRUE),
                     0) #Relative eggs transferred

  #Calculate fitness relating to post-mating selection (eggs that developed)
  tmp$prop_dev <- (tmp$NumDeveloped/tmp$MatingSuccess)/tmp$eggs_per_mate
  tmp$fit3 <- ifelse(tmp$MatingSuccess > 0,
                     tmp$prop_dev/mean(tmp$prop_dev, na.rm = TRUE),
                     0)
  
  #Standardizing the trait value to have a mean of 0 and sd of unity
  tmp$StdLength <- (tmp$svl - mean(tmp$svl))/sd(tmp$svl)
  
  #Calculating the absolute selection differentials (s)
  s1 <- cov(tmp$svl, tmp$fit1)
  s12 <- cov(tmp$svl, tmp$fit2)
  s123 <- cov(tmp$svl, tmp$fit3)
  s2 <- s12 - s1
  s3 <- s123 - s12
  
  #Calculating the standardized selection differentials (s')
  s1_prime <- cov(tmp$StdLength, tmp$fit1)
  s12_prime <- cov(tmp$StdLength, tmp$fit2)
  s123_prime <- cov(tmp$StdLength, tmp$fit3)
  s2_prime <- s12_prime - s1_prime
  s3_prime <- s123_prime - s12_prime
  
  #Combining all of the selection differentials (s, s') and saving the output
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  selection <- cbind(trial_num, s1, s2, s3, s12, s123, 
                     s1_prime, s2_prime, s3_prime, s12_prime, s123_prime)
  
  select_diff_fem <- rbind(select_diff_fem, selection)
  
  #Save the intermediate values
  fem_succ_select_diff <- rbind(fem_succ_select_diff, tmp)
}

#Exporting the data
write.csv(fem_succ_select_diff, "data/scovelli_int_diff_fem.csv", row.names = FALSE)
```

``` r
#Create a dataframe to store all of the intermediate values of fitness in
mal_succ_select_diff <- data.frame(matrix(ncol = ncol(mal_succSS) + 6,
                                          nrow = 0))
colnames(mal_succ_select_diff) <- c(colnames(mal_succSS),
                                    "fit1", "eggs_per_mate","fit2", "prop_dev", 
                                    "fit3", "StdLength")

#Create a dataframe to store the final calculations of I in
select_diff_mal <- data.frame(matrix(ncol = 11,
                                     nrow = 0))
colnames(select_diff_mal) <- c("trial", "s1", "s2", "s3", "s12", "s123",
                               "s1_prime", "s2_prime", "s3_prime", 
                               "s12_prime", "s123_prime")

for (trial in unique(mal_succSS$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- mal_succSS[mal_succSS$trial_num == trial, ]
  
  #Calculate fitness relating to pre-cop. selection (#matings)
  tmp$fit1 <- tmp$MatingSuccess/mean(tmp$MatingSuccess) #Relative mating success

  #Calculate fitness relating to post-mating selection (#eggs transferred)
  tmp$eggs_per_mate <- tmp$totalEggs/tmp$MatingSuccess
  tmp$fit2 <- ifelse(tmp$MatingSuccess > 0,
                     tmp$eggs_per_mate/mean(tmp$eggs_per_mate, na.rm = TRUE),
                     0) #Relative eggs transferred

  #Calculate fitness relating to post-mating selection (eggs that developed)
  tmp$prop_dev <- (tmp$NumDeveloped/tmp$MatingSuccess)/tmp$eggs_per_mate
  tmp$fit3 <- ifelse(tmp$MatingSuccess > 0,
                     tmp$prop_dev/mean(tmp$prop_dev, na.rm = TRUE),
                     0)
  
  #Standardizing the trait value to have a mean of 0 and sd of unity
  tmp$StdLength <- (tmp$svl - mean(tmp$svl))/sd(tmp$svl)
  
  #Calculating the absolute selection differentials (s)
  s1 <- cov(tmp$svl, tmp$fit1)
  s12 <- cov(tmp$svl, tmp$fit2)
  s123 <- cov(tmp$svl, tmp$fit3)
  s2 <- s12 - s1
  s3 <- s123 - s12
  
  #Calculating the standardized selection differentials (s')
  s1_prime <- cov(tmp$StdLength, tmp$fit1)
  s12_prime <- cov(tmp$StdLength, tmp$fit2)
  s123_prime <- cov(tmp$StdLength, tmp$fit3)
  s2_prime <- s12_prime - s1_prime
  s3_prime <- s123_prime - s12_prime
  
  #Combining all of the selection differentials (s, s') and saving the output
  trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
                               trial))
  selection <- cbind(trial_num, s1, s2, s3, s12, s123, 
                     s1_prime, s2_prime, s3_prime, s12_prime, s123_prime)
  
  select_diff_mal <- rbind(select_diff_mal, selection)
  
  #Save the intermediate values
  mal_succ_select_diff <- rbind(mal_succ_select_diff, tmp)
}

#Exporting the data
write.csv(mal_succ_select_diff, "data/fuscus_int_diff_mal.csv", row.names = FALSE)
```

``` r
#Merge the male and female datasets together
select_diff_fem$Sex <- "F"
select_diff_mal$Sex <- "M"

select_diff_all <- rbind(select_diff_fem, select_diff_mal)

write.csv(select_diff_all, "data/scovelli_select_diff.csv", row.names = FALSE)

#List the columns of interest
columns <- c("s1", "s2", "s3", "s123",
             "s1_prime", "s2_prime", "s3_prime", "s123_prime")

#Create a dataframe to store the final values in
sd_average <- data.frame(matrix(ncol = 4,
                                nrow = 0))
colnames(sd_average) <- c("Average", "Interval", "Select_diff", "Sex")

#Calculate the critical value
crit <- qt(p = 0.975, df = (nrow(select_diff_fem) - 1))

#Calculating the averages and confidence intervals for each species and 
#selection differential
for (j in 1:length(columns)) {
    
    col_name <- columns[[j]]
    
    #Calculate the means
    mean <- t(t(tapply(select_diff_all[, colnames(select_diff_all) == col_name], 
                       select_diff_all$Sex, mean)))
    
    #Calculate standard error
    se <- t(t(tapply(select_diff_all[, colnames(select_diff_all) == col_name], 
                     select_diff_all$Sex, 
                 function(x){
                   
                   sqrt(var(x))/sqrt(length(x))
                   
                 })))
    
    #Calculate the value that is added and subtracted from the mean
    int <- se*crit
    
    #Combine the data together
    episode <- as.data.frame(cbind(mean, int))
    colnames(episode) <- c("Average", "Interval")
    
    episode$Select_diff <- col_name
    episode$Sex <- rownames(episode)
    
    rownames(episode) <- NULL
    
    sd_average <- rbind(sd_average, episode)
    
}
```

Now that we have the average select diff for males and females alongside
the 95% confidence intervals we can visualize some results:

| Select_diff | F                   | M                  |
|:------------|:--------------------|:-------------------|
| s1          | 1.02 (0, 2.05)      | 0.01 (-0.44, 0.46) |
| s2          | -0.76 (-1.62, 0.11) | 0.14 (-0.11, 0.39) |
| s3          | 0.02 (-0.15, 0.18)  | 0.04 (-0.28, 0.36) |
| s123        | 0.28 (-0.06, 0.62)  | 0.19 (-0.26, 0.63) |
| s1_prime    | 0.34 (0.01, 0.68)   | 0.05 (-0.11, 0.21) |
| s2_prime    | -0.24 (-0.53, 0.04) | 0.03 (-0.07, 0.13) |
| s3_prime    | 0 (-0.07, 0.07)     | 0.03 (-0.08, 0.14) |
| s123_prime  | 0.1 (-0.02, 0.23)   | 0.11 (-0.06, 0.27) |

Average Selection Differentials (95% CI) for Males and Females

![*Absolute (left) and standardized (right) selection differentials for
male (purple) and female (green) S. fuscus. Error bars represent the 95%
confidence intervals around the
mean.*](selection_analysis_scovelli_files/figure-gfm/generate-fig-select-diff-1.png)

We can see from these results that males and females are experiencing
largely similar selection on snout-vent length except for in the first
episode of post-mating selection (number of eggs transferred).

Overall, it appears larger size in svl is more beneficial. However,
while we can see these patterns, there is only significant selection on
standardized snout-vent length (all the 95% CI cross 0) for females in
pre-mating selection ($s_1$).

## Looking into the Maximum Sexual Selection Differential

We can calculate $s'_{max}$:

$$
s'_{max}\ =\ \beta_{SS}\sqrt{I_S}
$$

One major benefit of this is that both the variation in mating success
and the Bateman gradient are taken into effect. Additionally, $s'_{max}$
places an upper bound on the strength of sexual selection during a
sexual selection episode.

Let’s now generate $s'_{max}$ for male and female *S. scovelli*.

``` r
#Pull out the Bateman gradient (slopes of the models)
bateman_fem <- coefficients(wls_model_fem)[2]
bateman_mal <- coefficients(wls_model_mal)[2]

#Pull out the Average opp sexual selection for males and females
opp_ss_fem <- opp_average$Average[opp_average$Sex == "F" & 
                                    opp_average$Episode_sel == "I_s"]
opp_ss_mal <- opp_average$Average[opp_average$Sex == "M" & 
                                    opp_average$Episode_sel == "I_s"]

#Calculate select diff max
select_diff_max_fem <- bateman_fem * sqrt(opp_ss_fem)
select_diff_max_mal <- bateman_mal * sqrt(opp_ss_mal)
```

For females, the $s'_{max}$ is 1.2740814, while $s'$ is 0.1016966, which
is about 1/12 of the max selection that could be experienced. For males,
$s'_{max}$ is 0.3860836 and $s'$ is 0.1063331 which is only about 1/3 of
the max selection that could be experienced.

# Visualizing post-copulatory selection

As a way to visualize selection acting AFTER the mating event
(post-copulatory selection) I am plotting the proportion of eggs that
survived against mating success. Hopefully this will tell us if
acquiring more mates is having any affect on the ability for the eggs to
develop.

I am going to plot this relationship only looking into the individuals
who mated.

![*Plotting the relationship between the proportion of eggs that
developed and the number of mates aquired for both males (purple) and
females
(green).*](selection_analysis_scovelli_files/figure-gfm/surv-v-matings-1.png)

There may be a correlation for females, let’s investigate further:

``` r
cor.test(fem_succSS$MatingSuccess[fem_succSS$MatingSuccess != 0],
         fem_succSS$prop_surviving[fem_succSS$totalEggs != 0])
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  fem_succSS$MatingSuccess[fem_succSS$MatingSuccess != 0] and fem_succSS$prop_surviving[fem_succSS$totalEggs != 0]
    ## t = -2.4919, df = 27, p-value = 0.01914
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.68962962 -0.07832723
    ## sample estimates:
    ##        cor 
    ## -0.4324193

There is a significantly negative correlation between mating success and
reproductive success for female *S. scovelli*.
