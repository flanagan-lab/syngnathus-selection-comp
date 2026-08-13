Bootstrapping selection estimates
================



- [Read in the datasets](#read-in-the-datasets)
- [Opportunity for selection](#opportunity-for-selection)
  - [Partitioning *I* for all
    individuals](#partitioning-i-for-all-individuals)
- [Selection differentials](#selection-differentials)

``` r
#This is a cohesive list of all the libraries used in this document
library(knitr)
# To install spfTools:
# devtools::install_github("https://github.com/spflanagan/spfTools/")
library(spfTools)
```

``` r
set.seed(2025)
num_bootstraps<-1000 #change this as needed
source("R/partition_I.R")
source("R/bootstrap_partition_I.R")
source("R/calc_selection_diffs.R")
source("R/bootstrap_s.R")
```

## Read in the datasets

``` r
fem_succFL <- read.csv("data/floridae_fem_succ.csv")
mal_succFL <- read.csv("data/floridae_mal_succ.csv")
```

``` r
fem_succFU <- read.csv("data/fuscus_fem_succ.csv")
mal_succFU <- read.csv("data/fuscus_mal_succ.csv")
```

``` r
fem_succSC <- read.csv("data/scovelli_fem_succ.csv")
mal_succSC <- read.csv("data/scovelli_mal_succ.csv")
```

## Opportunity for selection

### Partitioning *I* for all individuals

``` r
boot_floridae <- bootstrap_partition_I(num_bootstraps,
                                       fem_succFL,
                                       mal_succFL)
write.csv(boot_floridae, "data/floridae_bootstrapped_I_partitions.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_floridae<-read.csv("data/floridae_bootstrapped_I_partitions.csv")
```

``` r
boot_fuscus <- bootstrap_partition_I(num_bootstraps,
                                       fem_succFU,
                                       mal_succFU)
write.csv(boot_fuscus, "data/fuscus_bootstrapped_I_partitions.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_fuscus<-read.csv("data/fuscus_bootstrapped_I_partitions.csv")
```

``` r
boot_scovelli <- bootstrap_partition_I(num_bootstraps,
                                       fem_succSC,
                                       mal_succSC)
write.csv(boot_scovelli, "data/scovelli_bootstrapped_I_partitions.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_scovelli<-read.csv("data/scovelli_bootstrapped_I_partitions.csv")
```

#### Generating summary statistics from bootstraps

``` r
#Function to calculate confidence intervals from bootstrap results - employing percentile-based
#confidence intervals

calculate_CIs <- function(boot_results, conf_level = 0.95){
  
  alpha <- 1 - conf_level
  
  #Calculate percentile-based confidence intervals
  CIs <- apply(boot_results, 2, function(x){

    c(mean = mean(x),
      lower = quantile(x, alpha/2, na.rm = TRUE),
      upper = quantile(x, 1 - alpha/2, na.rm = TRUE))
    
  })
  
  return(t(CIs))
}
```

``` r
floridae_CIs <- calculate_CIs(boot_floridae)
fuscus_CIs <- calculate_CIs(boot_fuscus)
scovelli_CIs <- calculate_CIs(boot_scovelli)
```

``` r
female_CIs <- as.data.frame(rbind(floridae_CIs[grep("fem", rownames(floridae_CIs)), ],
                    fuscus_CIs[grep("fem", rownames(fuscus_CIs)), ],
                    scovelli_CIs[grep("fem", rownames(scovelli_CIs)), ]))

female_CIs$species <- c(rep("floridae", times = nrow(floridae_CIs[grep("fem", rownames(floridae_CIs)),])),
                        rep("fuscus", times = nrow(fuscus_CIs[grep("fem", rownames(fuscus_CIs)),])),
                        rep("scovelli", times = nrow(scovelli_CIs[grep("fem", rownames(scovelli_CIs)),])))

kable(round(female_CIs[,colnames(female_CIs) != "species"],6))
```

|                      |      mean | lower.2.5% | upper.97.5% |
|:---------------------|----------:|-----------:|------------:|
| I_1_fem              |  2.769904 |   1.317179 |    4.469052 |
| I_2_fem              |  0.142287 |   0.019819 |    0.284347 |
| coi1_2_fem           |  0.646908 |   0.448675 |    0.811199 |
| coi1_2given1_fem     | -0.040107 |  -0.171305 |    0.071605 |
| coi12_2_fem          |  0.789195 |   0.553450 |    0.996542 |
| coi12_2given1_fem    |  0.218287 |   0.026741 |    0.545730 |
| diff_12_fem          | -0.570908 |  -0.755247 |   -0.316249 |
| I_12_fem             |  2.948084 |   1.496197 |    4.706371 |
| I_3_fem              |  0.037796 |   0.001151 |    0.128762 |
| coi12_3_fem          |  0.664148 |   0.466648 |    0.815334 |
| coi12_3given2_fem    |  0.061706 |  -0.041057 |    0.201319 |
| coi123_3_fem         |  0.701944 |   0.484179 |    0.878165 |
| coi123_3given2_fem   |  0.155103 |  -0.007811 |    0.530414 |
| diff_123_fem         | -0.546841 |  -0.781240 |   -0.236131 |
| I_fem                |  3.164893 |   1.641990 |    4.903622 |
| I_1_fem.1            |  2.776615 |   1.365950 |    4.318333 |
| I_2_fem.1            |  0.036182 |   0.000204 |    0.125536 |
| coi1_2_fem.1         |  0.637906 |   0.403296 |    0.850000 |
| coi1_2given1_fem.1   |  0.007346 |   0.000000 |    0.023629 |
| coi12_2_fem.1        |  0.674088 |   0.431695 |    0.891692 |
| coi12_2given1_fem.1  |  0.087221 |   0.000408 |    0.252724 |
| diff_12_fem.1        | -0.586867 |  -0.826279 |   -0.342040 |
| I_12_fem.1           |  2.871181 |   1.497131 |    4.419234 |
| I_3_fem.1            |  0.000300 |   0.000000 |    0.000670 |
| coi12_3_fem.1        |  0.637013 |   0.402326 |    0.850646 |
| coi12_3given2_fem.1  |  0.002593 |  -0.001719 |    0.007432 |
| coi123_3_fem.1       |  0.637313 |   0.402999 |    0.851208 |
| coi123_3given2_fem.1 |  0.003178 |  -0.000726 |    0.008737 |
| diff_123_fem.1       | -0.634136 |  -0.849758 |   -0.399566 |
| I_fem.1              |  2.876952 |   1.510128 |    4.423255 |
| I_1_fem.2            |  1.793652 |   0.921824 |    3.099157 |
| I_2_fem.2            |  0.048545 |   0.022163 |    0.080368 |
| coi1_2_fem.2         |  0.509826 |   0.348427 |    0.664324 |
| coi1_2given1_fem.2   | -0.050495 |  -0.123410 |    0.005958 |
| coi12_2_fem.2        |  0.558371 |   0.400145 |    0.708964 |
| coi12_2given1_fem.2  |  0.037801 |  -0.022750 |    0.136803 |
| diff_12_fem.2        | -0.520570 |  -0.669581 |   -0.361070 |
| I_12_fem.2           |  1.780958 |   0.915341 |    3.048160 |
| I_3_fem.2            |  0.006277 |   0.001448 |    0.013424 |
| coi12_3_fem.2        |  0.525300 |   0.368940 |    0.668164 |
| coi12_3given2_fem.2  |  0.004739 |  -0.023838 |    0.039893 |
| coi123_3_fem.2       |  0.531576 |   0.374612 |    0.674015 |
| coi123_3given2_fem.2 |  0.018190 |  -0.015627 |    0.063624 |
| diff_123_fem.2       | -0.513387 |  -0.664373 |   -0.357792 |
| I_fem.2              |  1.803887 |   0.925266 |    3.061444 |

``` r
male_CIs <- as.data.frame(rbind(floridae_CIs[grep("mal", rownames(floridae_CIs)), ],
                  fuscus_CIs[grep("mal", rownames(fuscus_CIs)), ],
                  scovelli_CIs[grep("mal", rownames(scovelli_CIs)), ]))

male_CIs$species <- c(rep("floridae", times = nrow(floridae_CIs[grep("mal", rownames(floridae_CIs)),])),
                      rep("fuscus", times = nrow(fuscus_CIs[grep("mal", rownames(fuscus_CIs)),])),
                      rep("scovelli", times = nrow(scovelli_CIs[grep("mal", rownames(scovelli_CIs)),])))

male_CIs[,1:3]<-round(male_CIs[,1:3],6)
kable(male_CIs)
```

|                      |      mean | lower.2.5% | upper.97.5% | species  |
|:---------------------|----------:|-----------:|------------:|:---------|
| I_1_mal              |  2.462045 |   1.088524 |    4.241778 | floridae |
| I_2_mal              |  0.242932 |   0.042676 |    0.651593 | floridae |
| coi1_2_mal           |  0.574490 |   0.364829 |    0.734715 | floridae |
| coi1_2given1_mal     | -0.050284 |  -0.194790 |    0.088024 | floridae |
| coi12_2_mal          |  0.817423 |   0.549350 |    1.208123 | floridae |
| coi12_2given1_mal    |  0.411365 |   0.061422 |    1.038618 | floridae |
| diff_12_mal          | -0.406057 |  -0.670212 |   -0.017864 | floridae |
| I_12_mal             |  2.823127 |   1.354994 |    4.469692 | floridae |
| I_3_mal              |  0.036958 |   0.001456 |    0.120580 | floridae |
| coi12_3_mal          |  0.593693 |   0.391886 |    0.753157 | floridae |
| coi12_3given2_mal    |  0.040715 |  -0.052650 |    0.170490 | floridae |
| coi123_3_mal         |  0.630651 |   0.407503 |    0.803625 | floridae |
| coi123_3given2_mal   |  0.132065 |  -0.036535 |    0.486133 | floridae |
| diff_123_mal         | -0.498586 |  -0.735983 |   -0.187208 | floridae |
| I_mal                |  2.995907 |   1.425598 |    4.681028 | floridae |
| I_1_mal.1            |  2.758164 |   1.465155 |    4.135873 | fuscus   |
| I_2_mal.1            |  0.052594 |   0.000862 |    0.142914 | fuscus   |
| coi1_2_mal.1         |  0.640403 |   0.433333 |    0.833333 | fuscus   |
| coi1_2given1_mal.1   | -0.013746 |  -0.045992 |    0.000000 | fuscus   |
| coi12_2_mal.1        |  0.692997 |   0.475286 |    0.890739 | fuscus   |
| coi12_2given1_mal.1  |  0.120702 |   0.001653 |    0.336213 | fuscus   |
| diff_12_mal.1        | -0.572296 |  -0.802773 |   -0.308923 | fuscus   |
| I_12_mal.1           |  2.865119 |   1.601911 |    4.237304 | fuscus   |
| I_3_mal.1            |  0.001280 |   0.000001 |    0.003557 | fuscus   |
| coi12_3_mal.1        |  0.645317 |   0.438327 |    0.835956 | fuscus   |
| coi12_3given2_mal.1  |  0.007466 |  -0.001837 |    0.022238 | fuscus   |
| coi123_3_mal.1       |  0.646597 |   0.440519 |    0.838075 | fuscus   |
| coi123_3given2_mal.1 |  0.010018 |  -0.001044 |    0.029085 | fuscus   |
| diff_123_mal.1       | -0.636579 |  -0.831177 |   -0.429075 | fuscus   |
| I_mal.1              |  2.882603 |   1.633026 |    4.244581 | fuscus   |
| I_1_mal.2            |  0.176562 |   0.046647 |    0.513120 | scovelli |
| I_2_mal.2            |  0.087871 |   0.051752 |    0.134448 | scovelli |
| coi1_2_mal.2         |  0.103071 |   0.040816 |    0.187075 | scovelli |
| coi1_2given1_mal.2   |  0.000000 |   0.000000 |    0.000000 | scovelli |
| coi12_2_mal.2        |  0.190942 |   0.114479 |    0.279866 | scovelli |
| coi12_2given1_mal.2  |  0.096489 |   0.054936 |    0.148521 | scovelli |
| diff_12_mal.2        | -0.094454 |  -0.163478 |   -0.035954 | scovelli |
| I_12_mal.2           |  0.273051 |   0.124071 |    0.586807 | scovelli |
| I_3_mal.2            |  0.067074 |   0.013301 |    0.176254 | scovelli |
| coi12_3_mal.2        |  0.119365 |   0.047966 |    0.198556 | scovelli |
| coi12_3given2_mal.2  |  0.015298 |  -0.002205 |    0.039493 | scovelli |
| coi123_3_mal.2       |  0.186439 |   0.083490 |    0.323377 | scovelli |
| coi123_3given2_mal.2 |  0.080808 |   0.016995 |    0.205900 | scovelli |
| diff_123_mal.2       | -0.105631 |  -0.181201 |   -0.038511 | scovelli |
| I_mal.2              |  0.369157 |   0.181453 |    0.713790 | scovelli |

## Selection differentials

``` r
boot_floridae <- bootstrap_s(num_bootstraps,
                                       fem_succFL,
                                       mal_succFL)
write.csv(boot_floridae, "data/floridae_bootstrapped_s.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_floridae<-read.csv("data/floridae_bootstrapped_s.csv")
```

``` r
boot_fuscus <- bootstrap_s(num_bootstraps,
                                       fem_succFU,
                                       mal_succFU)
write.csv(boot_fuscus, "data/fuscus_bootstrapped_s.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_fuscus<-read.csv("data/fuscus_bootstrapped_s.csv")
```

``` r
boot_scovelli <- bootstrap_s(num_bootstraps,
                                       fem_succSC,
                                       mal_succSC)
write.csv(boot_scovelli, "data/scovelli_bootstrapped_s.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_scovelli<-read.csv("data/scovelli_bootstrapped_s.csv")
```

#### Generating summary statistics for selection differentials

``` r
floridae_CIs <- calculate_CIs(boot_floridae)
fuscus_CIs <- calculate_CIs(boot_fuscus)
scovelli_CIs <- calculate_CIs(boot_scovelli)
```

``` r
female_CIs <- as.data.frame(rbind(floridae_CIs[grep("fem", rownames(floridae_CIs)), ],
                    fuscus_CIs[grep("fem", rownames(fuscus_CIs)), ],
                    scovelli_CIs[grep("fem", rownames(scovelli_CIs)), ]))

female_CIs$species <- c(rep("floridae", times = nrow(floridae_CIs[grep("fem", rownames(floridae_CIs)),])),
                        rep("fuscus", times = nrow(fuscus_CIs[grep("fem", rownames(fuscus_CIs)),])),
                        rep("scovelli", times = nrow(scovelli_CIs[grep("fem", rownames(scovelli_CIs)),])))

kable(round(female_CIs[,colnames(female_CIs) != "species"],6))
```

|                  |      mean | lower.2.5% | upper.97.5% |
|:-----------------|----------:|-----------:|------------:|
| s1_fem           |  2.133775 |  -1.734951 |    6.334340 |
| s2_fem           | -2.143347 |  -5.303734 |    1.457109 |
| s3_fem           |  0.696840 |  -0.345598 |    1.813785 |
| s12_fem          | -0.009572 |  -1.585071 |    1.780581 |
| s123_fem         |  0.687267 |  -0.667753 |    2.067652 |
| s1_prime_fem     |  0.208629 |  -0.183208 |    0.599939 |
| s2_prime_fem     | -0.202862 |  -0.544145 |    0.135325 |
| s3_prime_fem     |  0.060289 |  -0.030323 |    0.150839 |
| s12_prime_fem    |  0.005767 |  -0.144067 |    0.186521 |
| s123_prime_fem   |  0.066056 |  -0.067314 |    0.196199 |
| s1_fem.1         |  2.139245 |   0.606776 |    3.683890 |
| s2_fem.1         | -1.580744 |  -2.877200 |   -0.344749 |
| s3_fem.1         | -0.012519 |  -0.144848 |    0.100207 |
| s12_fem.1        |  0.558501 |   0.042164 |    1.026939 |
| s123_fem.1       |  0.545982 |   0.084078 |    1.009678 |
| s1_prime_fem.1   |  0.434126 |   0.124474 |    0.736945 |
| s2_prime_fem.1   | -0.332289 |  -0.593602 |   -0.075584 |
| s3_prime_fem.1   |  0.002786 |  -0.019469 |    0.033067 |
| s12_prime_fem.1  |  0.101837 |   0.007846 |    0.187133 |
| s123_prime_fem.1 |  0.104623 |   0.021758 |    0.192363 |
| s1_fem.2         |  0.865463 |  -0.060077 |    1.870364 |
| s2_fem.2         | -0.634835 |  -1.443419 |    0.132140 |
| s3_fem.2         |  0.010415 |  -0.095005 |    0.116273 |
| s12_fem.2        |  0.230628 |  -0.104840 |    0.554641 |
| s123_fem.2       |  0.241043 |  -0.055714 |    0.552540 |
| s1_prime_fem.2   |  0.303439 |  -0.058973 |    0.676800 |
| s2_prime_fem.2   | -0.208775 |  -0.496868 |    0.069227 |
| s3_prime_fem.2   | -0.003703 |  -0.046062 |    0.036280 |
| s12_prime_fem.2  |  0.094664 |  -0.037590 |    0.225230 |
| s123_prime_fem.2 |  0.090960 |  -0.037613 |    0.214720 |

``` r
male_CIs <- as.data.frame(rbind(floridae_CIs[grep("mal", rownames(floridae_CIs)), ],
                  fuscus_CIs[grep("mal", rownames(fuscus_CIs)), ],
                  scovelli_CIs[grep("mal", rownames(scovelli_CIs)), ]))

male_CIs$species <- c(rep("floridae", times = nrow(floridae_CIs[grep("mal", rownames(floridae_CIs)),])),
                      rep("fuscus", times = nrow(fuscus_CIs[grep("mal", rownames(fuscus_CIs)),])),
                      rep("scovelli", times = nrow(scovelli_CIs[grep("mal", rownames(scovelli_CIs)),])))

male_CIs[,1:3]<-round(male_CIs[,1:3],6)
kable(male_CIs)
```

|                  |      mean | lower.2.5% | upper.97.5% | species  |
|:-----------------|----------:|-----------:|------------:|:---------|
| s1_mal           | -2.698896 |  -5.092572 |   -0.132139 | floridae |
| s2_mal           |  1.930512 |  -0.277298 |    4.071687 | floridae |
| s3_mal           | -0.418376 |  -1.072458 |    0.098118 | floridae |
| s12_mal          | -0.768384 |  -1.572385 |    0.134175 | floridae |
| s123_mal         | -1.186759 |  -1.894044 |   -0.468235 | floridae |
| s1_prime_mal     | -0.457204 |  -0.824185 |   -0.104706 | floridae |
| s2_prime_mal     |  0.329559 |   0.021811 |    0.671548 | floridae |
| s3_prime_mal     | -0.057883 |  -0.140687 |    0.010949 | floridae |
| s12_prime_mal    | -0.127645 |  -0.238950 |   -0.013372 | floridae |
| s123_prime_mal   | -0.185528 |  -0.287653 |   -0.086717 | floridae |
| s1_mal.1         | -0.614036 |  -2.209032 |    0.931892 | fuscus   |
| s2_mal.1         |  0.424356 |  -0.819361 |    1.648002 | fuscus   |
| s3_mal.1         |  0.010789 |  -0.052580 |    0.087109 | fuscus   |
| s12_mal.1        | -0.189680 |  -0.650162 |    0.287145 | fuscus   |
| s123_mal.1       | -0.178892 |  -0.631378 |    0.293832 | fuscus   |
| s1_prime_mal.1   | -0.189803 |  -0.529607 |    0.137317 | fuscus   |
| s2_prime_mal.1   |  0.138035 |  -0.099576 |    0.400359 | fuscus   |
| s3_prime_mal.1   |  0.000562 |  -0.023717 |    0.022673 | fuscus   |
| s12_prime_mal.1  | -0.051768 |  -0.153696 |    0.050421 | fuscus   |
| s123_prime_mal.1 | -0.051207 |  -0.152552 |    0.048431 | fuscus   |
| s1_mal.2         |  0.006749 |  -0.375404 |    0.321242 | scovelli |
| s2_mal.2         |  0.126983 |  -0.172755 |    0.450682 | scovelli |
| s3_mal.2         |  0.051408 |  -0.250024 |    0.374103 | scovelli |
| s12_mal.2        |  0.133733 |  -0.216121 |    0.485705 | scovelli |
| s123_mal.2       |  0.185140 |  -0.140124 |    0.533752 | scovelli |
| s1_prime_mal.2   |  0.047733 |  -0.083234 |    0.185796 | scovelli |
| s2_prime_mal.2   |  0.020855 |  -0.101119 |    0.132343 | scovelli |
| s3_prime_mal.2   |  0.033726 |  -0.063417 |    0.138662 | scovelli |
| s12_prime_mal.2  |  0.068587 |  -0.040649 |    0.183528 | scovelli |
| s123_prime_mal.2 |  0.102313 |  -0.003274 |    0.214213 | scovelli |
