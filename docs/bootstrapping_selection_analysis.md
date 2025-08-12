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
| I_1_fem | 2.72 | 2.82 | 2.73 | 2.82 | 1.76 | 1.83 |
| I_2_fem | 0.14 | 0.15 | 0.03 | 0.04 | 0.05 | 0.05 |
| coi1_2_fem | 0.64 | 0.65 | 0.63 | 0.64 | 0.50 | 0.51 |
| coi1_2given1_fem | -0.04 | -0.04 | 0.01 | 0.01 | -0.05 | -0.05 |
| coi12_2_fem | 0.78 | 0.80 | 0.67 | 0.68 | 0.55 | 0.56 |
| coi12_2given1_fem | 0.21 | 0.23 | 0.08 | 0.09 | 0.04 | 0.04 |
| diff_12_fem | -0.58 | -0.56 | -0.59 | -0.58 | -0.53 | -0.52 |
| I_12_fem | 2.90 | 3.00 | 2.83 | 2.92 | 1.75 | 1.82 |
| I_3_fem | 0.04 | 0.04 | 0.00 | 0.00 | 0.01 | 0.01 |
| coi12_3_fem | 0.66 | 0.67 | 0.63 | 0.64 | 0.52 | 0.53 |
| coi12_3given2_fem | 0.06 | 0.07 | 0.00 | 0.00 | 0.00 | 0.01 |
| coi123_3_fem | 0.70 | 0.71 | 0.63 | 0.64 | 0.53 | 0.54 |
| coi123_3given2_fem | 0.15 | 0.17 | 0.00 | 0.00 | 0.02 | 0.02 |
| diff_123_fem | -0.56 | -0.54 | -0.64 | -0.63 | -0.52 | -0.51 |
| I_fem | 3.11 | 3.22 | 2.83 | 2.92 | 1.77 | 1.84 |

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
| I_1_mal | 2.73 | 2.84 | 2.72 | 2.80 | 0.17 | 0.18 |
| I_2_mal | 0.13 | 0.14 | 0.05 | 0.06 | 0.09 | 0.09 |
| coi1_2_mal | 0.64 | 0.66 | 0.63 | 0.65 | 0.10 | 0.11 |
| coi1_2given1_mal | -0.04 | -0.04 | -0.01 | -0.01 | 0.00 | 0.00 |
| coi12_2_mal | 0.78 | 0.80 | 0.69 | 0.70 | 0.19 | 0.19 |
| coi12_2given1_mal | 0.21 | 0.22 | 0.11 | 0.13 | 0.09 | 0.10 |
| diff_12_mal | -0.58 | -0.57 | -0.58 | -0.56 | -0.10 | -0.09 |
| I_12_mal | 2.91 | 3.02 | 2.82 | 2.91 | 0.27 | 0.28 |
| I_3_mal | 0.04 | 0.04 | 0.00 | 0.00 | 0.06 | 0.07 |
| coi12_3_mal | 0.66 | 0.67 | 0.64 | 0.65 | 0.12 | 0.12 |
| coi12_3given2_mal | 0.06 | 0.07 | 0.01 | 0.01 | 0.01 | 0.02 |
| coi123_3_mal | 0.70 | 0.71 | 0.64 | 0.65 | 0.18 | 0.19 |
| coi123_3given2_mal | 0.15 | 0.17 | 0.01 | 0.01 | 0.08 | 0.08 |
| diff_123_mal | -0.55 | -0.54 | -0.64 | -0.63 | -0.11 | -0.10 |
| I_mal | 3.13 | 3.24 | 2.84 | 2.92 | 0.36 | 0.38 |

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

Because some of the fuscus trials only had one male and/or one female
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
