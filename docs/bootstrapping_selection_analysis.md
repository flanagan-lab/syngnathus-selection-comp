Bootstrapping selection estimates
================



- [Read in the datasets](#read-in-the-datasets)
- [Opportunity for selection](#opportunity-for-selection)
  - [Partitioning *I* for all
    individuals](#partitioning-i-for-all-individuals)
  - [Partitioning *I* for only mated
    individuals](#partitioning-i-for-only-mated-individuals)
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
| I_1_mal              |  2.786641 |   1.379741 |    4.550087 | floridae |
| I_2_mal              |  0.138900 |   0.019109 |    0.274878 | floridae |
| coi1_2_mal           |  0.649312 |   0.456925 |    0.810739 | floridae |
| coi1_2given1_mal     | -0.039080 |  -0.167801 |    0.064899 | floridae |
| coi12_2_mal          |  0.788212 |   0.564657 |    0.994310 | floridae |
| coi12_2given1_mal    |  0.215568 |   0.026597 |    0.556546 | floridae |
| diff_12_mal          | -0.572643 |  -0.761159 |   -0.330327 | floridae |
| I_12_mal             |  2.963129 |   1.508897 |    4.769433 | floridae |
| I_3_mal              |  0.039895 |   0.001106 |    0.132617 | floridae |
| coi12_3_mal          |  0.665947 |   0.469292 |    0.821304 | floridae |
| coi12_3given2_mal    |  0.061680 |  -0.054640 |    0.194069 | floridae |
| coi123_3_mal         |  0.705841 |   0.498385 |    0.869893 | floridae |
| coi123_3given2_mal   |  0.161062 |  -0.015020 |    0.525905 | floridae |
| diff_123_mal         | -0.544780 |  -0.799817 |   -0.215633 | floridae |
| I_mal                |  3.185871 |   1.652233 |    4.981343 | floridae |
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

### Partitioning *I* for only mated individuals

``` r
boot_floridae <- bootstrap_partition_I(num_bootstraps,
                                       fem_succFL[fem_succFL$mated==1,],
                                       mal_succFL[mal_succFL$mated==1,])
write.csv(boot_floridae, "data/floridae_bootstrapped_I_partitions_mated.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_floridae<-read.csv("data/floridae_bootstrapped_I_partitions_mated.csv")
```

Because some of the S. fuscus trials only had one male and/or one female
mate (which cannot have a variance calculated on it for partitioning I,
if all non-mated individuals are removed), I first need to sub-set the
data to only include trials with 2 or more mated individuals.

``` r
fem_fu_mated<-fem_succFU[fem_succFU$mated==1,]
mal_fu_mated <- mal_succFU[mal_succFU$MatingSuccess>0,]
fem_tab<-table(fem_fu_mated$trial_num)
mal_tab<-table(mal_fu_mated$trial_num)

trials_to_use<-intersect(names(fem_tab[which(fem_tab>1)]),
                         names(mal_tab[which(mal_tab>1)]))

boot_fuscus <- bootstrap_partition_I(num_bootstraps,
                                       fem_fu_mated[which(fem_fu_mated$trial_num %in% trials_to_use),],
                                       mal_fu_mated[which(mal_fu_mated$trial_num %in% trials_to_use),])
write.csv(boot_fuscus, "data/fuscus_bootstrapped_I_partitions_mated.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_fuscus<-read.csv("data/fuscus_bootstrapped_I_partitions_mated.csv")
```

``` r
boot_scovelli <- bootstrap_partition_I(num_bootstraps,
                                       fem_succSC[fem_succSC$mated==1,],
                                       mal_succSC[mal_succSC$MatingSuccess>0,])
write.csv(boot_scovelli, "data/scovelli_bootstrapped_I_partitions_mated.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_scovelli<-read.csv("data/scovelli_bootstrapped_I_partitions_mated.csv")
```

#### Generating summary statistics from bootstraps (only mated)

``` r
floridae_CIs<-apply(boot_floridae,2,cim)
fuscus_CIs<-apply(boot_fuscus,2,cim)
scovelli_CIs<-apply(boot_scovelli,2,cim)
```

``` r
female_CIs<-cbind(t(floridae_CIs[,grep("fem",colnames(floridae_CIs))]),
      t(fuscus_CIs[,grep("fem",colnames(fuscus_CIs))]),
      t(scovelli_CIs[,grep("fem",colnames(scovelli_CIs))]))
colnames(female_CIs)<-c("floridae_low","floridae_upp",
                        "fuscus_low","fuscus_upp",
                        "scovelli_low","scovelli_upp")
kable(round(female_CIs,2))
```

|  | floridae_low | floridae_upp | fuscus_low | fuscus_upp | scovelli_low | scovelli_upp |
|:---|---:|---:|---:|---:|---:|---:|
| I_1_fem | 0.08 | 0.08 | 0.02 | 0.02 | 0.15 | 0.15 |
| I_2_fem | 0.19 | 0.20 | 0.13 | 0.14 | 0.06 | 0.06 |
| coi1_2_fem | -0.03 | -0.03 | 0.01 | 0.02 | -0.03 | -0.03 |
| coi1_2given1_fem | -0.03 | -0.02 | 0.01 | 0.01 | -0.03 | -0.03 |
| coi12_2_fem | 0.17 | 0.18 | 0.14 | 0.16 | 0.02 | 0.03 |
| coi12_2given1_fem | 0.14 | 0.15 | 0.14 | 0.15 | 0.02 | 0.03 |
| diff_12_fem | -0.03 | -0.02 | 0.00 | 0.00 | 0.00 | 0.00 |
| I_12_fem | 0.19 | 0.20 | 0.17 | 0.18 | 0.14 | 0.14 |
| I_3_fem | 0.06 | 0.07 | 0.00 | 0.00 | 0.01 | 0.01 |
| coi12_3_fem | 0.04 | 0.04 | 0.01 | 0.01 | 0.00 | 0.00 |
| coi12_3given2_fem | 0.04 | 0.04 | 0.00 | 0.00 | 0.00 | 0.00 |
| coi123_3_fem | 0.10 | 0.11 | 0.01 | 0.01 | 0.01 | 0.01 |
| coi123_3given2_fem | 0.09 | 0.10 | 0.00 | 0.01 | 0.01 | 0.01 |
| diff_123_fem | -0.01 | -0.01 | 0.00 | 0.00 | 0.00 | 0.00 |
| I_fem | 0.33 | 0.35 | 0.18 | 0.19 | 0.15 | 0.16 |

``` r
male_CIs<-cbind(t(floridae_CIs[,grep("mal",colnames(floridae_CIs))]),
      t(fuscus_CIs[,grep("mal",colnames(fuscus_CIs))]),
      t(scovelli_CIs[,grep("mal",colnames(scovelli_CIs))]))
colnames(male_CIs)<-c("floridae_low","floridae_upp",
                        "fuscus_low","fuscus_upp",
                        "scovelli_low","scovelli_upp")
kable(round(male_CIs,2))
```

|  | floridae_low | floridae_upp | fuscus_low | fuscus_upp | scovelli_low | scovelli_upp |
|:---|---:|---:|---:|---:|---:|---:|
| I_1_mal | 0.08 | 0.08 | 0.02 | 0.02 | 0.00 | 0.00 |
| I_2_mal | 0.18 | 0.20 | 0.16 | 0.17 | 0.09 | 0.09 |
| coi1_2_mal | -0.03 | -0.02 | -0.03 | -0.03 | 0.00 | 0.00 |
| coi1_2given1_mal | -0.03 | -0.02 | -0.03 | -0.02 | 0.00 | 0.00 |
| coi12_2_mal | 0.16 | 0.17 | 0.13 | 0.14 | 0.09 | 0.09 |
| coi12_2given1_mal | 0.14 | 0.14 | 0.13 | 0.14 | 0.09 | 0.09 |
| diff_12_mal | -0.03 | -0.02 | 0.00 | 0.00 | 0.00 | 0.00 |
| I_12_mal | 0.19 | 0.20 | 0.12 | 0.13 | 0.09 | 0.09 |
| I_3_mal | 0.06 | 0.06 | 0.00 | 0.00 | 0.06 | 0.07 |
| coi12_3_mal | 0.04 | 0.04 | 0.00 | 0.00 | 0.02 | 0.02 |
| coi12_3given2_mal | 0.04 | 0.04 | 0.00 | 0.00 | 0.01 | 0.01 |
| coi123_3_mal | 0.09 | 0.10 | 0.00 | 0.00 | 0.08 | 0.09 |
| coi123_3given2_mal | 0.09 | 0.10 | 0.00 | 0.00 | 0.07 | 0.08 |
| diff_123_mal | -0.01 | 0.00 | 0.00 | 0.00 | -0.01 | -0.01 |
| I_mal | 0.32 | 0.34 | 0.12 | 0.14 | 0.17 | 0.18 |

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
floridae_CIs<-apply(boot_floridae,2,cim)
fuscus_CIs<-apply(boot_fuscus,2,cim)
scovelli_CIs<-apply(boot_scovelli,2,cim)
```

``` r
female_CIs<-cbind(t(floridae_CIs[,grep("fem",colnames(floridae_CIs))]),
      t(fuscus_CIs[,grep("fem",colnames(fuscus_CIs))]),
      t(scovelli_CIs[,grep("fem",colnames(scovelli_CIs))]))
colnames(female_CIs)<-c("floridae_low","floridae_upp",
                        "fuscus_low","fuscus_upp",
                        "scovelli_low","scovelli_upp")
kable(round(female_CIs,2))
```

|  | floridae_low | floridae_upp | fuscus_low | fuscus_upp | scovelli_low | scovelli_upp |
|:---|---:|---:|---:|---:|---:|---:|
| s1_fem | 1.79 | 2.05 | 2.10 | 2.20 | 0.83 | 0.89 |
| s2_fem | -2.06 | -1.84 | -1.63 | -1.55 | -0.66 | -0.61 |
| s3_fem | 0.65 | 0.72 | -0.01 | -0.01 | 0.01 | 0.02 |
| s12_fem | -0.09 | 0.03 | 0.54 | 0.57 | 0.22 | 0.24 |
| s123_fem | 0.61 | 0.70 | 0.53 | 0.56 | 0.23 | 0.25 |
| s1_prime_fem | 0.18 | 0.20 | 0.43 | 0.45 | 0.29 | 0.31 |
| s2_prime_fem | -0.20 | -0.18 | -0.35 | -0.33 | -0.22 | -0.20 |
| s3_prime_fem | 0.06 | 0.06 | 0.00 | 0.00 | 0.00 | 0.00 |
| s12_prime_fem | 0.00 | 0.01 | 0.10 | 0.10 | 0.09 | 0.10 |
| s123_prime_fem | 0.06 | 0.07 | 0.10 | 0.11 | 0.08 | 0.09 |

``` r
male_CIs<-cbind(t(floridae_CIs[,grep("mal",colnames(floridae_CIs))]),
      t(fuscus_CIs[,grep("mal",colnames(fuscus_CIs))]),
      t(scovelli_CIs[,grep("mal",colnames(scovelli_CIs))]))
colnames(male_CIs)<-c("floridae_low","floridae_upp",
                        "fuscus_low","fuscus_upp",
                        "scovelli_low","scovelli_upp")
kable(round(male_CIs,2))
```

|  | floridae_low | floridae_upp | fuscus_low | fuscus_upp | scovelli_low | scovelli_upp |
|:---|---:|---:|---:|---:|---:|---:|
| s1_mal | 2.05 | 2.31 | -0.66 | -0.56 | 0.00 | 0.03 |
| s2_mal | -2.23 | -2.02 | 0.39 | 0.47 | 0.11 | 0.13 |
| s3_mal | 0.62 | 0.69 | 0.01 | 0.01 | 0.05 | 0.07 |
| s12_mal | 0.00 | 0.11 | -0.19 | -0.16 | 0.13 | 0.15 |
| s123_mal | 0.67 | 0.76 | -0.18 | -0.16 | 0.19 | 0.21 |
| s1_prime_mal | 0.20 | 0.23 | -0.19 | -0.17 | 0.05 | 0.06 |
| s2_prime_mal | -0.21 | -0.19 | 0.13 | 0.14 | 0.02 | 0.02 |
| s3_prime_mal | 0.05 | 0.06 | 0.00 | 0.00 | 0.03 | 0.04 |
| s12_prime_mal | 0.01 | 0.02 | -0.05 | -0.05 | 0.07 | 0.07 |
| s123_prime_mal | 0.06 | 0.07 | -0.05 | -0.05 | 0.10 | 0.11 |
