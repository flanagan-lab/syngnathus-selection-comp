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
source("R/bootstrap_s.R")
source("R/calc_selection_diffs.R")
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
write.csv(boot_floridae, "data/floridae_bootstrapped_I_partitions_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_floridae<-read.csv("data/floridae_bootstrapped_I_partitions_R1.csv")
```

``` r
boot_fuscus <- bootstrap_partition_I(num_bootstraps,
                                       fem_succFU,
                                       mal_succFU)
write.csv(boot_fuscus, "data/fuscus_bootstrapped_I_partitions_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_fuscus<-read.csv("data/fuscus_bootstrapped_I_partitions_R1.csv")
```

``` r
boot_scovelli <- bootstrap_partition_I(num_bootstraps,
                                       fem_succSC,
                                       mal_succSC)
write.csv(boot_scovelli, "data/scovelli_bootstrapped_I_partitions_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_scovelli<-read.csv("data/scovelli_bootstrapped_I_partitions_R1.csv")
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
female_CIs[,colnames(female_CIs) != "species"]<-round(female_CIs[,colnames(female_CIs) != "species"],6)

kable(female_CIs)
```

|                      |      mean | lower.2.5% | upper.97.5% | species  |
|:---------------------|----------:|-----------:|------------:|:---------|
| I_1_fem              |  2.799501 |   0.000000 |    7.000000 | floridae |
| I_2_fem              |  0.194217 |   0.000000 |    0.716442 | floridae |
| coi1_2_fem           |  0.692449 |   0.000000 |    1.000000 | floridae |
| coi1_2given1_fem     | -0.024335 |  -0.682802 |    0.569276 | floridae |
| coi12_2_fem          |  0.886666 |   0.000000 |    1.366611 | floridae |
| coi12_2given1_fem    |  0.431135 |  -0.103038 |    1.769293 | floridae |
| diff_12_fem          | -0.455531 |  -1.000000 |    0.438153 | floridae |
| I_12_fem             |  3.206301 |   0.000000 |    7.000000 | floridae |
| I_3_fem              |  0.042162 |   0.000000 |    0.266712 | floridae |
| coi12_3_fem          |  0.715675 |   0.000000 |    1.000000 | floridae |
| coi12_3given2_fem    |  0.063701 |  -0.147706 |    0.531980 | floridae |
| coi123_3_fem         |  0.757837 |   0.000000 |    1.148851 | floridae |
| coi123_3given2_fem   |  0.159675 |  -0.074837 |    1.370745 | floridae |
| diff_123_fem         | -0.598161 |  -1.000000 |    0.223639 | floridae |
| I_fem                |  3.429677 |   0.000000 |    7.000000 | floridae |
| I_1_fem.1            |  4.144510 |   0.000000 |   12.000000 | fuscus   |
| I_2_fem.1            |  0.260180 |   0.000000 |    0.785593 | fuscus   |
| coi1_2_fem.1         |  0.772795 |   0.000000 |    1.000000 | fuscus   |
| coi1_2given1_fem.1   | -0.035363 |  -0.688368 |    0.415877 | fuscus   |
| coi12_2_fem.1        |  1.032975 |   0.000000 |    1.573512 | fuscus   |
| coi12_2given1_fem.1  |  0.876150 |   0.000000 |    3.130805 | fuscus   |
| diff_12_fem.1        | -0.156825 |  -1.000000 |    1.611298 | fuscus   |
| I_12_fem.1           |  4.985297 |   0.000000 |   12.000000 | fuscus   |
| I_3_fem.1            |  0.004042 |   0.000000 |    0.021873 | fuscus   |
| coi12_3_fem.1        |  0.774672 |   0.000000 |    1.000000 | fuscus   |
| coi12_3given2_fem.1  | -0.014267 |  -0.298341 |    0.154471 | fuscus   |
| coi123_3_fem.1       |  0.778714 |   0.000000 |    1.000000 | fuscus   |
| coi123_3given2_fem.1 |  0.000906 |  -0.215268 |    0.175778 | fuscus   |
| diff_123_fem.1       | -0.777808 |  -1.031653 |    0.000000 | fuscus   |
| I_fem.1              |  4.971936 |   0.000000 |   12.000000 | fuscus   |
| I_1_fem.2            |  1.750087 |   0.333284 |    7.000000 | scovelli |
| I_2_fem.2            |  0.081180 |   0.000000 |    0.271112 | scovelli |
| coi1_2_fem.2         |  0.525467 |   0.067651 |    1.000000 | scovelli |
| coi1_2given1_fem.2   | -0.085825 |  -0.456341 |    0.138017 | scovelli |
| coi12_2_fem.2        |  0.606647 |   0.169912 |    1.000000 | scovelli |
| coi12_2given1_fem.2  |  0.061369 |  -0.121554 |    0.386688 | scovelli |
| diff_12_fem.2        | -0.545278 |  -1.000000 |   -0.140626 | scovelli |
| I_12_fem.2           |  1.725631 |   0.275622 |    7.000000 | scovelli |
| I_3_fem.2            |  0.010618 |   0.000000 |    0.044410 | scovelli |
| coi12_3_fem.2        |  0.552329 |   0.152053 |    1.000000 | scovelli |
| coi12_3given2_fem.2  | -0.002354 |  -0.088010 |    0.110832 | scovelli |
| coi123_3_fem.2       |  0.562948 |   0.162297 |    1.000000 | scovelli |
| coi123_3given2_fem.2 |  0.017062 |  -0.067138 |    0.187331 | scovelli |
| diff_123_fem.2       | -0.545886 |  -1.000000 |   -0.159028 | scovelli |
| I_fem.2              |  1.740338 |   0.270934 |    7.000000 | scovelli |

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
| I_1_mal              |  2.711623 |   0.000000 |    7.000000 | floridae |
| I_2_mal              |  0.195575 |   0.000000 |    0.680822 | floridae |
| coi1_2_mal           |  0.677108 |   0.000000 |    1.000000 | floridae |
| coi1_2given1_mal     | -0.015699 |  -0.640554 |    0.509907 | floridae |
| coi12_2_mal          |  0.872683 |   0.000000 |    1.366619 | floridae |
| coi12_2given1_mal    |  0.446754 |  -0.080241 |    1.770821 | floridae |
| diff_12_mal          | -0.425929 |  -1.000000 |    0.438059 | floridae |
| I_12_mal             |  3.142678 |   0.000000 |    7.000000 | floridae |
| I_3_mal              |  0.035913 |   0.000000 |    0.189329 | floridae |
| coi12_3_mal          |  0.697727 |   0.000000 |    1.000000 | floridae |
| coi12_3given2_mal    |  0.050313 |  -0.199931 |    0.472162 | floridae |
| coi123_3_mal         |  0.733640 |   0.000000 |    1.036507 | floridae |
| coi123_3given2_mal   |  0.126406 |  -0.097080 |    0.890994 | floridae |
| diff_123_mal         | -0.607234 |  -1.000000 |    0.000000 | floridae |
| I_mal                |  3.319397 |   0.000000 |    7.000000 | floridae |
| I_1_mal.1            |  3.900907 |   0.779221 |   12.000000 | fuscus   |
| I_2_mal.1            |  0.262798 |   0.000000 |    0.794869 | fuscus   |
| coi1_2_mal.1         |  0.769311 |   0.454545 |    1.000000 | fuscus   |
| coi1_2given1_mal.1   | -0.106004 |  -1.021374 |    0.000000 | fuscus   |
| coi12_2_mal.1        |  1.032109 |   0.727130 |    1.495198 | fuscus   |
| coi12_2given1_mal.1  |  0.820117 |   0.000000 |    2.749957 | fuscus   |
| diff_12_mal.1        | -0.211992 |  -1.000000 |    1.337112 | fuscus   |
| I_12_mal.1           |  4.615020 |   1.355505 |   12.000000 | fuscus   |
| I_3_mal.1            |  0.008015 |   0.000000 |    0.056250 | fuscus   |
| coi12_3_mal.1        |  0.778969 |   0.456717 |    1.000000 | fuscus   |
| coi12_3given2_mal.1  |  0.018141 |  -0.187801 |    0.236136 | fuscus   |
| coi123_3_mal.1       |  0.786983 |   0.461273 |    1.000000 | fuscus   |
| coi123_3given2_mal.1 |  0.046435 |  -0.136112 |    0.478695 | fuscus   |
| diff_123_mal.1       | -0.740548 |  -1.000000 |   -0.303344 | fuscus   |
| I_mal.1              |  4.679596 |   1.397455 |   12.000000 | fuscus   |
| I_1_mal.2            |  0.131931 |   0.000000 |    0.466667 | scovelli |
| I_2_mal.2            |  0.123775 |   0.028925 |    0.288723 | scovelli |
| coi1_2_mal.2         |  0.105333 |   0.000000 |    0.333333 | scovelli |
| coi1_2given1_mal.2   |  0.000000 |   0.000000 |    0.000000 | scovelli |
| coi12_2_mal.2        |  0.229109 |   0.041024 |    0.508330 | scovelli |
| coi12_2given1_mal.2  |  0.137692 |   0.032969 |    0.325976 | scovelli |
| diff_12_mal.2        | -0.091416 |  -0.306314 |    0.000000 | scovelli |
| I_12_mal.2           |  0.269623 |   0.041024 |    0.692735 | scovelli |
| I_3_mal.2            |  0.065255 |   0.000307 |    0.276851 | scovelli |
| coi12_3_mal.2        |  0.131535 |  -0.011225 |    0.374708 | scovelli |
| coi12_3given2_mal.2  |  0.028709 |  -0.015748 |    0.124889 | scovelli |
| coi123_3_mal.2       |  0.196791 |  -0.001895 |    0.600231 | scovelli |
| coi123_3given2_mal.2 |  0.092790 |  -0.005811 |    0.438603 | scovelli |
| diff_123_mal.2       | -0.104000 |  -0.333082 |    0.013119 | scovelli |
| I_mal.2              |  0.391123 |   0.053025 |    1.103472 | scovelli |

### Partitioning *I* for only mated individuals

``` r
boot_floridae <- bootstrap_partition_I(num_bootstraps,
                                       fem_succFL[fem_succFL$mated==1,],
                                       mal_succFL[mal_succFL$mated==1,])
write.csv(boot_floridae, "data/floridae_bootstrapped_I_partitions_mated_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_floridae<-read.csv("data/floridae_bootstrapped_I_partitions_mated_R1.csv")
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
write.csv(boot_fuscus, "data/fuscus_bootstrapped_I_partitions_mated_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_fuscus<-read.csv("data/fuscus_bootstrapped_I_partitions_mated_R1.csv")
```

``` r
boot_scovelli <- bootstrap_partition_I(num_bootstraps,
                                       fem_succSC[fem_succSC$mated==1,],
                                       mal_succSC[mal_succSC$MatingSuccess>0,])
write.csv(boot_scovelli, "data/scovelli_bootstrapped_I_partitions_mated_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_scovelli<-read.csv("data/scovelli_bootstrapped_I_partitions_mated_R1.csv")
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
| I_1_fem | 0.18 | 0.19 | 0.06 | 0.06 | 0.19 | 0.20 |
| I_2_fem | 0.26 | 0.28 | 0.36 | 0.37 | 0.08 | 0.09 |
| coi1_2_fem | -0.01 | -0.01 | -0.04 | -0.03 | -0.05 | -0.05 |
| coi1_2given1_fem | 0.00 | 0.01 | -0.06 | -0.05 | -0.05 | -0.04 |
| coi12_2_fem | 0.26 | 0.27 | 0.32 | 0.33 | 0.03 | 0.03 |
| coi12_2given1_fem | 0.28 | 0.30 | 0.30 | 0.32 | 0.03 | 0.04 |
| diff_12_fem | 0.03 | 0.03 | -0.02 | -0.02 | 0.00 | 0.00 |
| I_12_fem | 0.46 | 0.49 | 0.31 | 0.32 | 0.18 | 0.19 |
| I_3_fem | 0.04 | 0.04 | 0.00 | 0.00 | 0.01 | 0.01 |
| coi12_3_fem | 0.07 | 0.08 | 0.01 | 0.02 | 0.00 | 0.00 |
| coi12_3given2_fem | 0.05 | 0.05 | 0.01 | 0.01 | 0.00 | 0.00 |
| coi123_3_fem | 0.11 | 0.12 | 0.02 | 0.02 | 0.01 | 0.01 |
| coi123_3given2_fem | 0.08 | 0.08 | 0.01 | 0.01 | 0.01 | 0.02 |
| diff_123_fem | -0.04 | -0.03 | -0.01 | 0.00 | 0.00 | 0.00 |
| I_fem | 0.59 | 0.63 | 0.33 | 0.35 | 0.19 | 0.21 |

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
| I_1_mal | 0.18 | 0.19 | 0.06 | 0.06 | 0.00 | 0.00 |
| I_2_mal | 0.26 | 0.27 | 0.36 | 0.38 | 0.12 | 0.12 |
| coi1_2_mal | -0.01 | 0.00 | -0.07 | -0.06 | 0.00 | 0.00 |
| coi1_2given1_mal | 0.00 | 0.01 | -0.10 | -0.08 | 0.00 | 0.00 |
| coi12_2_mal | 0.25 | 0.27 | 0.31 | 0.32 | 0.12 | 0.12 |
| coi12_2given1_mal | 0.29 | 0.30 | 0.30 | 0.31 | 0.12 | 0.12 |
| diff_12_mal | 0.03 | 0.04 | -0.01 | -0.01 | 0.00 | 0.00 |
| I_12_mal | 0.47 | 0.50 | 0.27 | 0.28 | 0.12 | 0.12 |
| I_3_mal | 0.04 | 0.05 | 0.00 | 0.00 | 0.05 | 0.06 |
| coi12_3_mal | 0.07 | 0.08 | 0.01 | 0.01 | 0.02 | 0.03 |
| coi12_3given2_mal | 0.05 | 0.05 | 0.01 | 0.01 | 0.02 | 0.02 |
| coi123_3_mal | 0.11 | 0.12 | 0.01 | 0.01 | 0.07 | 0.09 |
| coi123_3given2_mal | 0.08 | 0.09 | 0.01 | 0.01 | 0.06 | 0.08 |
| diff_123_mal | -0.04 | -0.03 | 0.00 | 0.00 | -0.01 | -0.01 |
| I_mal | 0.60 | 0.64 | 0.30 | 0.30 | 0.20 | 0.22 |

## Selection differentials

``` r
boot_floridae <- bootstrap_s(num_bootstraps,
                                       fem_succFL,
                                       mal_succFL)
write.csv(boot_floridae, "data/floridae_bootstrapped_s_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_floridae<-read.csv("data/floridae_bootstrapped_s_R1.csv")
```

``` r
boot_fuscus <- bootstrap_s(num_bootstraps,
                                       fem_succFU,
                                       mal_succFU)
write.csv(boot_fuscus, "data/fuscus_bootstrapped_s_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_fuscus<-read.csv("data/fuscus_bootstrapped_s_R1.csv")
```

``` r
boot_scovelli <- bootstrap_s(num_bootstraps,
                                       fem_succSC,
                                       mal_succSC)
write.csv(boot_scovelli, "data/scovelli_bootstrapped_s_R1.csv",
          quote=FALSE, row.names = FALSE)
```

``` r
boot_scovelli<-read.csv("data/scovelli_bootstrapped_s_R1.csv")
```

#### Generating summary statistics for selection differentials

``` r
floridae_CIs<-calculate_CIs(boot_floridae)
fuscus_CIs<-calculate_CIs(boot_fuscus)
scovelli_CIs<-calculate_CIs(boot_scovelli)
```

``` r
female_CIs<-cbind(floridae_CIs[grep("fem",rownames(floridae_CIs)),],
                  fuscus_CIs[grep("fem",rownames(fuscus_CIs)),],
                  scovelli_CIs[grep("fem",rownames(scovelli_CIs)),])
colnames(female_CIs)[1:3]<-paste0("floridae_",colnames(female_CIs)[1:3])
colnames(female_CIs)[4:6]<-paste0("fuscus",colnames(female_CIs)[4:6])
colnames(female_CIs)[7:9]<-paste0("scovelli_",colnames(female_CIs)[7:9])

kable(round(female_CIs,2))
```

|  | floridae_mean | floridae_lower.2.5% | floridae_upper.97.5% | fuscusmean | fuscuslower.2.5% | fuscusupper.97.5% | scovelli_mean | scovelli_lower.2.5% | scovelli_upper.97.5% |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| s1_fem | 2.19 | -11.59 | 14.78 | 0.94 | -5.88 | 7.16 | NA | -1.71 | 5.50 |
| s2_fem | -2.48 | -12.20 | 8.59 | -0.51 | -6.22 | 5.07 | NA | -4.07 | 1.38 |
| s3_fem | 0.50 | -1.25 | 3.62 | -0.18 | -1.15 | 0.48 | -0.11 | -0.56 | 0.23 |
| s12_fem | -0.29 | -5.09 | 4.01 | 0.43 | -1.11 | 2.01 | 0.56 | -0.65 | 1.80 |
| s123_fem | 0.20 | -4.12 | 4.06 | 0.25 | -1.08 | 1.55 | 0.45 | -0.65 | 1.50 |
| s1_prime_fem | 0.21 | -0.91 | 1.34 | 0.16 | -1.00 | 1.24 | NA | -0.60 | 1.41 |
| s2_prime_fem | -0.23 | -1.02 | 0.66 | -0.09 | -1.11 | 0.86 | NA | -1.11 | 0.55 |
| s3_prime_fem | 0.05 | -0.12 | 0.30 | -0.03 | -0.20 | 0.08 | -0.04 | -0.19 | 0.09 |
| s12_prime_fem | -0.02 | -0.44 | 0.34 | 0.07 | -0.19 | 0.36 | 0.19 | -0.23 | 0.49 |
| s123_prime_fem | 0.02 | -0.36 | 0.37 | 0.04 | -0.18 | 0.28 | 0.15 | -0.25 | 0.43 |

``` r
male_CIs<-cbind(floridae_CIs[grep("mal",rownames(floridae_CIs)),],
      fuscus_CIs[grep("mal",rownames(floridae_CIs)),],
      scovelli_CIs[grep("mal",rownames(floridae_CIs)),])
colnames(male_CIs)[1:3]<-paste0("floridae_",colnames(male_CIs)[1:3])
colnames(male_CIs)[4:6]<-paste0("fuscus_",colnames(male_CIs)[4:6])
colnames(male_CIs)[7:9]<-paste0("scovelli_",colnames(male_CIs)[7:9])
kable(round(male_CIs,2))
```

|  | floridae_mean | floridae_lower.2.5% | floridae_upper.97.5% | fuscus_mean | fuscus_lower.2.5% | fuscus_upper.97.5% | scovelli_mean | scovelli_lower.2.5% | scovelli_upper.97.5% |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| s1_mal | 1.89 | -12.48 | 14.40 | -1.51 | -7.89 | 4.69 | 0.08 | -1.36 | 1.03 |
| s2_mal | -2.28 | -11.41 | 9.55 | 1.40 | -3.94 | 6.57 | 0.46 | -0.53 | 1.61 |
| s3_mal | 0.46 | -1.13 | 3.45 | -0.23 | -1.00 | 0.42 | -0.15 | -1.31 | 1.23 |
| s12_mal | -0.38 | -5.30 | 3.51 | -0.11 | -1.70 | 1.48 | 0.54 | -0.81 | 1.95 |
| s123_mal | 0.08 | -4.09 | 3.80 | -0.34 | -1.76 | 1.03 | 0.40 | -1.04 | 2.00 |
| s1_prime_mal | 0.19 | -0.97 | 1.22 | -0.24 | -1.18 | 0.78 | 0.03 | -0.31 | 0.31 |
| s2_prime_mal | -0.22 | -1.01 | 0.72 | 0.22 | -0.64 | 1.01 | 0.11 | -0.15 | 0.36 |
| s3_prime_mal | 0.04 | -0.12 | 0.29 | -0.04 | -0.17 | 0.07 | -0.03 | -0.31 | 0.31 |
| s12_prime_mal | -0.03 | -0.46 | 0.33 | -0.02 | -0.26 | 0.27 | 0.14 | -0.19 | 0.45 |
| s123_prime_mal | 0.01 | -0.37 | 0.35 | -0.05 | -0.28 | 0.17 | 0.11 | -0.24 | 0.51 |
