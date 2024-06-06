Selection pressures in *Syngnathus floridae*
================



- <a href="#calculating-the-degree-of-sexual-dimorphism"
  id="toc-calculating-the-degree-of-sexual-dimorphism">Calculating the
  degree of sexual dimorphism</a>
  - <a href="#checking-the-assumptions-for-a-pairwise-comparison"
    id="toc-checking-the-assumptions-for-a-pairwise-comparison">Checking the
    assumptions for a pairwise comparison</a>
  - <a href="#investigate-distributions-and-run-the-tests"
    id="toc-investigate-distributions-and-run-the-tests">Investigate
    distributions and run the tests</a>
- <a
  href="#calculating-mating-and-reproductive-success-for-individuals-who-mated"
  id="toc-calculating-mating-and-reproductive-success-for-individuals-who-mated">Calculating
  mating and reproductive success for individuals who mated</a>
- <a href="#summary-statistics-for-successfully-mated-individuals"
  id="toc-summary-statistics-for-successfully-mated-individuals">Summary
  statistics for successfully mated individuals</a>
  - <a href="#males" id="toc-males">Males</a>
  - <a href="#females" id="toc-females">Females</a>

``` r
#This is a cohesive list of all the libraries used in this document
library(ggplot2)
library(cowplot)
library(fBasics)
library(pwr)
```

``` r
#MomIDs and embryo counts for each section of the male's brood pouch
em_dat <- read.csv("data/EmbryoParentage_floridae.csv")

#Metadata for males and females from the mesocosm experiments
fem_mesoFL <- read.csv("data/all_fem_meso_floridae.csv")
mal_mesoFL <- read.csv("data/all_mal_meso_floridae.csv")
```

# Calculating the degree of sexual dimorphism

Other papers have reported varying levels of significant or
non-significant size differences between males and females in this
species, particularly in terms of standard length (measured from the tip
of the snout to the base of the caudal fin). I want to investigate the
data that I have to see what sexual size dimorphism is like for this
population. I am doing this across all fish that were used, including
those trials that had no successful matings.

The aspects of size that I am interested in looking at include:

- **Standard length** (A; yellow line): measured from the tip of the
  snout to the base of the caudal fin.
- **Torso depth** (B; blue line): measured from the point in front of
  the dorsal fin down (perpendicular to standard length).
- **Snout-vent length** (C; red line): measured from the tip of the
  snout to the urogenital opening.

<p float="center">

<img src="../figs/morph_fig.png" style="width:950px;"/>

</p>

## Checking the assumptions for a pairwise comparison

Before comparing male size to female size I want to ensure the data
meets the assumptions. The main two things that I will be looking into
include:

1.  Equal variances between my groups (using `var.test()`).
2.  Normal distribution of the data (using `normalTest()`).

The null hypotheses for these two tests are that the variances are equal
and the data is normally distributed.

To account for the fact that fish who are longer may just inherently be
deeper as well, I am going to adjust the depth by the standard length of
the pipefish prior to funning any analyses

``` r
#Adjust the torso depth
fem_mesoFL$depth_adj <- fem_mesoFL$depth/fem_mesoFL$length
mal_mesoFL$depth_adj <- mal_mesoFL$depth/mal_mesoFL$length

#Testing to see if the variances are equal
var.test(fem_mesoFL$length, mal_mesoFL$length) #EQUAL
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_mesoFL$length and mal_mesoFL$length
    ## F = 1.4609, num df = 85, denom df = 85, p-value = 0.08239
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.9521741 2.2413401
    ## sample estimates:
    ## ratio of variances 
    ##           1.460872

``` r
var.test(fem_mesoFL$depth_adj, mal_mesoFL$depth_adj) #EQUAL
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_mesoFL$depth_adj and mal_mesoFL$depth_adj
    ## F = 1.3686, num df = 84, denom df = 85, p-value = 0.151
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.8911382 2.1030826
    ## sample estimates:
    ## ratio of variances 
    ##           1.368614

``` r
var.test(fem_mesoFL$svl, mal_mesoFL$svl) #NOT EQUAL
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  fem_mesoFL$svl and mal_mesoFL$svl
    ## F = 2.1705, num df = 85, denom df = 85, p-value = 0.0004335
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  1.414675 3.330029
    ## sample estimates:
    ## ratio of variances 
    ##           2.170463

``` r
#Testing for normal distribution - Females
normalTest(fem_mesoFL$length, method = "da") #NOT NORMAL
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 6.0342
    ##     Z3  | Skewness: 0.5032
    ##     Z4  | Kurtosis: -2.4044
    ##   P VALUE:
    ##     Omnibus  Test: 0.04894 
    ##     Skewness Test: 0.6148 
    ##     Kurtosis Test: 0.0162

``` r
normalTest(fem_mesoFL$depth_adj, method = "da") #NORMAL
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 0.6513
    ##     Z3  | Skewness: 0.7106
    ##     Z4  | Kurtosis: -0.3825
    ##   P VALUE:
    ##     Omnibus  Test: 0.7221 
    ##     Skewness Test: 0.4773 
    ##     Kurtosis Test: 0.7021

``` r
normalTest(fem_mesoFL$svl, method = "da") #NORMAL
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 3.0383
    ##     Z3  | Skewness: 0.7065
    ##     Z4  | Kurtosis: -1.5935
    ##   P VALUE:
    ##     Omnibus  Test: 0.2189 
    ##     Skewness Test: 0.4799 
    ##     Kurtosis Test: 0.1111

``` r
#Testing for normal distribution - Males
normalTest(mal_mesoFL$length, method = "da") #NORMAL
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 0.0501
    ##     Z3  | Skewness: -0.1257
    ##     Z4  | Kurtosis: 0.1851
    ##   P VALUE:
    ##     Omnibus  Test: 0.9753 
    ##     Skewness Test: 0.9 
    ##     Kurtosis Test: 0.8531

``` r
normalTest(mal_mesoFL$depth_adj, method = "da") #NORMAL
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 1.9819
    ##     Z3  | Skewness: 0.5381
    ##     Z4  | Kurtosis: 1.3009
    ##   P VALUE:
    ##     Omnibus  Test: 0.3712 
    ##     Skewness Test: 0.5905 
    ##     Kurtosis Test: 0.1933

``` r
normalTest(mal_mesoFL$svl, method = "da") #NORMAL
```

    ## 
    ## Title:
    ##  D'Agostino Normality Test
    ## 
    ## Test Results:
    ##   STATISTIC:
    ##     Chi2 | Omnibus: 3.4451
    ##     Z3  | Skewness: 0.5003
    ##     Z4  | Kurtosis: -1.7874
    ##   P VALUE:
    ##     Omnibus  Test: 0.1786 
    ##     Skewness Test: 0.6168 
    ##     Kurtosis Test: 0.07387

## Investigate distributions and run the tests

Looking at the distributions below I do not expect to see differences in
standard length, but maybe in depth and SVL. I will run one of three
tests depending on the results of the above assumption testing:

1.  If **both** assumptions are met: *Two Sample t-test*
2.  If only the **variances** are not equal, but the data is normal:
    *Welch Two Sample t-test*
3.  If the **data is not normal**, regardless of whether variances are
    equal: *Wilcoxon rank sum test with continuity correction*

Based on this, I will run a Wilcoxon test for standard length, a Welch
two sample t-test for snout-vent length, and a two sample t-test for
torso depth (adjusted).

![Histograms of male and female pipefish body
sizes.](selection_analysis_floridae_files/figure-gfm/histogram_sizes-1.png)

``` r
#Running the appropriate test
wilcox.test(fem_mesoFL$length, mal_mesoFL$length)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  fem_mesoFL$length and mal_mesoFL$length
    ## W = 3930, p-value = 0.4784
    ## alternative hypothesis: true location shift is not equal to 0

``` r
t.test(fem_mesoFL$depth_adj, mal_mesoFL$depth_adj, var.equal = TRUE)
```

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  fem_mesoFL$depth_adj and mal_mesoFL$depth_adj
    ## t = 4.1064, df = 169, p-value = 6.248e-05
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.001494861 0.004262830
    ## sample estimates:
    ##  mean of x  mean of y 
    ## 0.03537844 0.03249959

``` r
t.test(fem_mesoFL$svl, mal_mesoFL$svl, var.equal = FALSE)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  fem_mesoFL$svl and mal_mesoFL$svl
    ## t = 3.86, df = 149.61, p-value = 0.0001682
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  2.949949 9.137632
    ## sample estimates:
    ## mean of x mean of y 
    ##  80.01955  73.97576

From these results we can see that there are no significant differences
between male and female pipefish in terms of standard length. However,
female dusky pipefish are significantly deeper that males and possess a
significantly longer snout-vent length. For those two morphometrics I am
going to run a power test to ensure that we are confident in our
results.

``` r
#Checking the power - SVL
d_mean_svl <- abs(mean(fem_mesoFL$svl, na.rm = TRUE) - 
                    mean(mal_mesoFL$svl, na.rm = TRUE))
pool_sd_svl <- sqrt((var(fem_mesoFL$svl, na.rm = TRUE) + 
                       var(mal_mesoFL$svl, na.rm = TRUE))/ 2)
d_svl <- d_mean_svl/pool_sd_svl

pwr.t.test(n = nrow(fem_mesoFL), 
           d = d_svl,
           sig.level = 0.05,
           type = 'two.sample',
           alternative = 'two.sided')
```

    ## 
    ##      Two-sample t test power calculation 
    ## 
    ##               n = 86
    ##               d = 0.5886437
    ##       sig.level = 0.05
    ##           power = 0.9698186
    ##     alternative = two.sided
    ## 
    ## NOTE: n is number in *each* group

``` r
#Checking the power - Depth
d_mean_depth <- abs(mean(fem_mesoFL$depth_adj, na.rm = TRUE) - 
                      mean(mal_mesoFL$depth_adj, na.rm = TRUE))
pool_sd_depth <- sqrt((var(fem_mesoFL$depth_adj, na.rm = TRUE) + 
                         var(mal_mesoFL$depth_adj, na.rm = TRUE))/ 2)
d_depth <- d_mean_depth/pool_sd_depth
pwr.t.test(n = nrow(fem_mesoFL), 
           d = d_depth,
           sig.level = 0.05,
           type = 'two.sample',
           alternative = 'two.sided')
```

    ## 
    ##      Two-sample t test power calculation 
    ## 
    ##               n = 86
    ##               d = 0.627763
    ##       sig.level = 0.05
    ##           power = 0.9835455
    ##     alternative = two.sided
    ## 
    ## NOTE: n is number in *each* group

For both variables we have a power of over 0.9 or over 90% so we can be
confident in our interpretation.

# Calculating mating and reproductive success for individuals who mated

*Syngnathus floridae* (dusky pipefish) were sampled from three distinct
seagrass beds around Tampa Bay in Tampa, Florida. Sexually mature
females (standard length $\ge$ 120mm) and pregnant males were collected
and brought back to the University of Tampa and put into experimental
breeding populations. In these trials, 8 males and 8 females were housed
together in a 140L tank for a period of 14-days and allowed to mate
freely. Parentage analysis was done with all of the pregnant males from
the trials to figure out how many times each male and female mated, and
the number of eggs that were transferred. The results of that are here.

First I had to calculate the mating and reproductive success for each
male and female who mated based on the assigned mom for each genotyped
embryo.

``` r
#Row-by-Row analysis of parentage data by male brood pouch section

#Read in the data
#em_dat <- read.csv("~/EmbryoParentage.csv")

#For each row in the dataset(each section of the pouch) apply this function
mom_counts <- do.call(rbind,apply(em_dat, 1, function(one_section){
  
  #Save all of the momIDs into an object
  mom_ids<-c(one_section[grep("momID",names(one_section))])  
  
  #Calculate the number of eggs that belongs to each potential mom based on
  #the proportions and total number of developed and undeveloped embryos
  mom_props<-c(as.numeric(one_section[grep("prop",names(one_section))]))
  mom_counts_dev<-mom_props*as.numeric(one_section["num_embryos_dev"])
  mom_counts_und<-mom_props*as.numeric(one_section["num_embryos_non_dev"])
  
  #Create a dataframe that contains the maleID, pouch section number and the
  #number of eggs that belongs to each momID
  this_section<-data.frame(
    maleID=one_section["maleID"],
    section_num=one_section["section_num"],
    mom_ids[which((mom_counts_dev + mom_counts_und) > 0)],
    mom_counts_dev[which((mom_counts_dev + mom_counts_und)>0)],
    mom_counts_und[which((mom_counts_dev + mom_counts_und)>0)]
  )
  
  #Rename the columns
  colnames(this_section)[3:5]<-c("momID","num_dev","num_und")
  
  return(this_section)
  
}))

#Calculate female fitness
fem_fitness<-do.call(rbind,by(mom_counts, mom_counts$momID,function(dat){
  
  mom_fitness<-data.frame(
    momID=unique(dat$momID),
    MatingSuccess=length(unique(dat$maleID)),
    NumDeveloped=round(sum(dat$num_dev)),
    NumUndeveloped=round(sum(dat$num_und))
  )
  return(mom_fitness)
}))

fem_fitness$totalEggs <- fem_fitness$NumDeveloped + fem_fitness$NumUndeveloped

#Calculate Male Fitness 
mal_fitness<-do.call(rbind,by(mom_counts, mom_counts$maleID,function(dat){
 
  dad_fitness<-data.frame(
    maleID=unique(dat$maleID),
    MatingSuccess=length(unique(dat$momID)),
    NumDeveloped_Calc=round(sum(dat$num_dev)),
    NumUndeveloped_Calc=round(sum(dat$num_und))
  )
  return(dad_fitness)
}))

mal_fitness$totalEggs <- mal_fitness$NumDeveloped_Calc + mal_fitness$NumUndeveloped_Calc
```

After running the above R script we have two datasets, `mal_fitness` and
`fem_fitness`. These datasets include information about the mating
success (number of mates) and reproductive success (Number of embryos
transferred). We can split reproductive success up further later if we
want to from the total number of embryos transferred to the number of
embryos developed and the number that were undeveloped.

I want to include all of the other metadata that I have for these
individuals (traits, collection location, latency to pregnancy, etc.) as
well as tack on all of the information for the individuals who did not
mate. To do that I am going to need to merge the fitness datasets with
`fem_meso` and `mal_meso`.

``` r
#Make a column in *_meso that contains the full fishID (i.e. FL1M3) to match the 
#formatting in the fitness datasets (make sure they have the same name for merging purposes)
fem_mesoFL$momID <- paste0("FL", fem_mesoFL$trial_num, "F", fem_mesoFL$fishID)
mal_mesoFL$maleID <- paste0("FL", mal_mesoFL$trial_num, "M", mal_mesoFL$fishID)

#Merge the datasets based on the columns created above
fem_allFL <- merge(fem_mesoFL, fem_fitness, by = "momID", all.x = TRUE, all.y = TRUE)
mal_allFL <- merge(mal_mesoFL, mal_fitness, by = "maleID", all.x = TRUE, all.y = TRUE)
```

There are a few trials that I want to remove from the analysis:

1.  All trials where there were no successful matings (7, 9, 10, 11).

2.  Trial 1, a male gave birth and the babies were immediately eaten by
    the adults so the trial was ended early and therefore I was unable
    to get any parentage information for that trial.

I also want to replace the NAs that were automatically added to the
columns from the fitness dataset (MatingSuccess, NumDeveloped,
NumUndeveloped, totalEggs) with 0s and add a column to the female
dataset that tells me whether or not the female mated (with 1 or 0).

``` r
#Subset the merged datasets to remove trials without successful matings and Trial 1
fem_succFL <- subset(fem_allFL, !(trial_num %in% c(7, 9, 10, 11, 1)))
mal_succFL <- subset(mal_allFL, !(trial_num %in% c(7, 9, 10, 11, 1)))

#Replace NAs with 0s in the columns related to fitness
mal_succFL[, c("MatingSuccess", "NumDeveloped_Calc", 
               "NumUndeveloped_Calc", "totalEggs")] <- sapply(mal_succFL[, c("MatingSuccess", 
                                                                             "NumDeveloped_Calc", 
                                                                             "NumUndeveloped_Calc", 
                                                                             "totalEggs")],
                                                              function(x)
                                                                ifelse(is.na(x), 0, x))

fem_succFL[, c("MatingSuccess", "NumDeveloped", 
               "NumUndeveloped", "totalEggs")] <- sapply(fem_succFL[, c("MatingSuccess", 
                                                                        "NumDeveloped", 
                                                                        "NumUndeveloped", 
                                                                        "totalEggs")],
                                                         function(x)
                                                           ifelse(is.na(x), 0, x))

#Add a column for females to denote mated or unmated
fem_succFL$mated <- ifelse(fem_succFL$MatingSuccess > 0, 1, 0)
```

# Summary statistics for successfully mated individuals

## Males

Across all 7 trials and 56 total males, there were 24 males that mated
at least one time and 6 of those males had two mates.

Looking across all males, including the ones that did not mate, this is
what we find as the mean, sd, and se for the number of embryos
transferred and how many of those developed versus didn’t:

|                     |       mean |          SD |         SE | max | min |
|:--------------------|-----------:|------------:|-----------:|----:|----:|
| Number of Embryos   |      108.5 | 155.1902059 | 20.7381636 | 500 |   0 |
| Developed Embryos   | 94.1964286 | 139.7347383 | 18.6728398 | 463 |   0 |
| Undeveloped Embryos | 14.3035714 |  27.7907242 |  3.7136917 | 131 |   0 |

These values will be influenced by the number of 0s coming from males
who did not mate. So let’s look at the same thing, but this time for
only males who had at least one successful mating:

|                     |        mean |          SD |         SE | max | min |
|:--------------------|------------:|------------:|-----------:|----:|----:|
| Number of Embryos   | 253.1666667 | 139.1941611 | 28.4128892 | 500 |   9 |
| Developed Embryos   | 219.7916667 | 133.7427704 | 27.3001287 | 463 |   7 |
| Undeveloped Embryos |      33.375 |  34.3901054 |  7.0198509 | 131 |   0 |

We can see from the bottom table that even when we only include males
who mated there is still a wide range in the brood size. I want to see
what relationship there is between brood pouch size (in terms of both
total area and length) and brood size (total number of embryos).

![Scatterplot of the relationship between brood pouch size metrics and
the number of embryos a male
had.](selection_analysis_floridae_files/figure-gfm/em-v-bp-1.png)

There may be some correlation happening here, but it doesn’t look
particularly strong. Let’s run some correlations tests to see what they
say.

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(mated_malFL$bp_area) and mated_malFL$totalEggs
    ## t = 0.90485, df = 22, p-value = 0.3753
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2316769  0.5507364
    ## sample estimates:
    ##       cor 
    ## 0.1894228

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(mated_malFL$bp_length) and mated_malFL$totalEggs
    ## t = 0.51929, df = 22, p-value = 0.6087
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3069838  0.4916139
    ## sample estimates:
    ##       cor 
    ## 0.1100398

There is not a significant correlation between the number of eggs and
size of the brood pouch when we look at brood pouch area OR brood pouch
length.

Multiple of the wild study papers looked at correlations between body
size in terms of standard length and the number of embryos and found
significant positive correlations.

![Scatterplot of the relationship between standard length (mm) and the
number of embryos a male
had.](selection_analysis_floridae_files/figure-gfm/em-v-sl-1.png)

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(mated_malFL$length) and mated_malFL$totalEggs
    ## t = -0.35647, df = 22, p-value = 0.7249
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.4649642  0.3379463
    ## sample estimates:
    ##         cor 
    ## -0.07578142

There appears to be no significant correlation between standard length
and the number of eggs in males, unlike what has been found in the
species previously. This is unsurprising as we didn’t find any hint of a
relationship in the brood pouch size metrics and there was a
considerable number of males that did not have a full brood pouch.

## Females

Across all 7 trials and 56 total females, there were 21 females that
mated at least one time, 5 females that mated twice, and 2 that mated 3
times.

Looking across all females, including the ones that did not mate, this
is what we find as the mean, sd, and se for the total number of embryos
transferred from each female (across all of her mates if applicable) and
how many of those developed versus didn’t:

|                     |        mean |          SD |         SE | max | min |
|:--------------------|------------:|------------:|-----------:|----:|----:|
| Number of Embryos   | 108.4821429 | 195.2328206 |  26.089083 | 916 |   0 |
| Developed Embryos   |  94.1785714 | 177.0691808 | 23.6618646 | 785 |   0 |
| Undeveloped Embryos |  14.3035714 |  29.7033329 |  3.9692748 | 131 |   0 |

These values will be influenced by the number of 0s coming from females
who did not mate. So let’s look at the same thing, but this time for
only females who had at least one successful mating:

|                     |        mean |          SD |         SE | max | min |
|:--------------------|------------:|------------:|-----------:|----:|----:|
| Number of Embryos   | 289.2857143 | 223.3819919 |  48.745947 | 916 |  35 |
| Developed Embryos   | 251.1428571 | 211.7324457 | 46.2038076 | 785 |   7 |
| Undeveloped Embryos |  38.1428571 |   38.360508 |  8.3709491 | 131 |   0 |

We can see from the bottom table that even when we only include females
who mated there is still a wide range in the number of eggs transferred.
I want to see what relationship there may be between female body size
(in terms of standard length, depth, and SVL) and the number of eggs she
transferred. I also want to see on average how many eggs were
transferred per mating. I’m going to calculate this by taking the total
number of eggs and dividing it by the number of mates.

    ## [1] 202.3254

    ## [1] 14.68782

![Scatterplot of the relationship between female size metrics and the
number of eggs
transferred.](selection_analysis_floridae_files/figure-gfm/em-v-fem-size-1.png)

There may be some correlation happening here, but it doesn’t look
particularly strong. Let’s run some correlations tests to see what they
say.

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_femFL$length and as.numeric(mated_femFL$totalEggs)
    ## t = -0.92304, df = 19, p-value = 0.3676
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5864091  0.2465727
    ## sample estimates:
    ##        cor 
    ## -0.2071652

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_femFL$depth_adj and as.numeric(mated_femFL$totalEggs)
    ## t = 2.0729, df = 19, p-value = 0.05202
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.002716801  0.726473616
    ## sample estimates:
    ##       cor 
    ## 0.4294737

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mated_femFL$svl and as.numeric(mated_femFL$totalEggs)
    ## t = -0.5112, df = 19, p-value = 0.6151
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5219227  0.3318958
    ## sample estimates:
    ##        cor 
    ## -0.1164796

There is not a significant correlation between the number of eggs and
size of the female in terms of standard length, depth, or snout-vent
length. Interestingly, however, there is a negative correlation for
length and SVL and a positive correlation for depth (but they are all
overall weak).
