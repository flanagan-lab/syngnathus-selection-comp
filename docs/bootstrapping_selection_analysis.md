Bootstrapping selection estimates
================



- [Read in the datasets](#read-in-the-datasets)
- [Opportunity for selection](#opportunity-for-selection)
- [Generating summary statistics from
  bootstraps](#generating-summary-statistics-from-bootstraps)

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

## Generating summary statistics from bootstraps

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
