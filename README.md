# Measuring Sexual Selection in _Syngnathus_ pipefish

This is a repository for the analysis of sexual selection pressues in three species of pipefish from the genus _Syngnathus_. This includes the dusky pipefish _Syngnathus floridae_, the Northern pipefish _Syngnathus fuscus_, and the Gulf pipefish _Syngnathus scovelli_. The goals of this study are the following:

  1. Investigate the contributions of pre- and post-copulatory sexual selection across all three species.
  2. Generate Bateman's gradient for male and females within each species.
  3. Attempt to identify traits that may be targets of pre- and post-copulatory selection using selection differentials.

These analyses support the findings in the manuscript "Comparing mating systems and sexual selection pressures across three congeneric species of pipefish that span the continuum of sexual dimorphism", currently submitted to a journal for consideration.

## Data
The RMarkdown documents refer to data that is found in the data/ directory at the top of the repo. 

### Data used to calculate estimates of selection
The datasets in this section were made from the original raw data (see "Data Availability" for how to access the raw data).

  - `all_fem_meso_*.csv` and `all_mal_meso_*.csv`: These files contain the data about the morphometrics of all male and female pipefish and various information relating to the trials they were in.
     
     - For the `_scovelli` datasets, the columns **svl**, **depth**, **bp_area**, and **bp_length** were not included in the original dataset. The original photos were obtained and those additional measurements were gathered. See the "Data Availability" section for more detail about the changes made to the publicly accessible _S. Scovelli_ dataset.
       
     - Variable description for these datasets:
       | Variable | Description | Species the variable is found |
       | -------- | ----------- | ----------------------------- |
       | trial_num & fishID | Indicates the trial the pipefish was used in and the unique identifier it was given (1 - 8 for _S. scovelli_ and _S. floridae_, 1 - 6 for _S. fuscus_) | All three |
       | col_location | Indication of which location the pipefish was captured, either TBML (27.8730210, -82.5397645), FD (27.6273507, -82.7002240), or DC (28.0516324, -82.7921436) | _S. floridae_ |
       | col_month | The month the pipefish was collected during | _S. floridae_ and _S. fuscus_ |
       | lat_to_trial | The latency to entering a trial (i.e., how long the pipefish was in the lab prior to entering a trial) measured as days | _S. floridae_ and _S. fuscus_ |
       | weight | Weight in grams of each pipefish | _S. floridae_ and _S. fuscus_ |
       | length, depth & svl | The standard length (measured from the tip of the snout to the end of the caudal fin), the torso depth, and snout-vent length (measured from the tip of the snout to the urogenital opening) in millimetres | All three |
       | bp_area & bp_length | Area and length of the male brood pouch in millimetres | MALE pipefish for all three |
       | snout_len | Length in millimetres of the snout | _S. fuscus_ |
       | preg_status | Indicates pregnant or non-pregnant as either a 1 or 0 (does NOT indicate number of mates) | _S. fuscus_ and _S. floridae_ MALE |
       | lat_to_preg | How long after starting the trial it took for the male to become pregnant in terms of days | MALE pipefish for all three |
       | MatingSuccess | Total number of mates obtained by that pipefish | _S. scovelli_ |
       | mate_cat | Mating category stating whether the pipefish mated once (s) or more than once (m) | _S. scovelli_ FEMALES |
       | totalEggs | Total number of eggs stored/transferred across all matings | _S. scovelli_ |
       | NumDeveloped | The number of eggs stored/transferred across all matings which showed signs of development | _S. scovelli_ |
       | NumUndeveloped | The number of eggs stored/transferred across all matings which did NOT show signs of development | _S. scovelli_ |
       | per_surviving | The proportion of total eggs stored/transferred which showed some signs of development | _S. scovelli_ |
       

  - `EmbryoParentage_*.csv`: This is the file that is used to calculate reproductive fitness and mating success for male and female pipefish. Each row represents a section of a male's pouch and contains information about who the confirmed mother was (`momID_XX`) for each of the genotyped embryos (`babyID_XX`) and the total number of developed (`num_embryos_dev`) and undeveloped (`num_embryos_non_dev`) embryos in that section.
    
      - _*_ There will be one male and female .csv file for only two of the three species (`_floridae` and `_fuscus`). _S. scovelli_ does not have a file in this format as the publically accesible dataset already contains information about reproductive fitness and mating success.

### Datasets containing fitness estimates generated using the all_XX_meso_XX.csv and EmbryoParentage_XX.csv datasets
The .csv files outlined here are ones which were generated in the `selection_analysis_*.Rmd` documents (see below) and then used in either the `cross_species_comp.Rmd` or `bootstrapping_selection_analysis.Rmd` documents.

  - `*_fem_bateman.csv`: There will be one file for _S. floridae_ (FL_XXX), _S. fuscus_ (FU_XXX), and _S. scovelli_ (SS_XXX). These files contain the relative mating success (`MatingSuccess`) and relative number of developed eggs (`rel_repo_fitness`) for female pipefish of the three species. This data was used to generate the Bateman gradients.
    
  - `episode_select_data.csv`: This dataset includes the opportunity for selection averaged across the trials (`average_cal`) for each episode of selection (`episode_sel`) for both sexes (`sex`) and all three species (`species`). Lower and upper 95% confidence intervals are also included (`lower` and `upper`). The means and CIs presented here are generated from the bootstrapping analysis.
    
  - `fem_fitness_*.csv` and `mal_fitness_*.csv`: These datasets are an extension of the `all_fem_meso_*.csv` and `all_mal_meso_*.csv` datasets described above with additional information about reproductive and mating success and subset to only include trials where at least one individual has mated. These datasets were used to generate information about the Bateman Gradient and for the bootstrapping analyses.
    
     - Variables not already described for these datasets:
       | Variable | Description | Species the variable is found |
       | -------- | ----------- | ----------------------------- |
       | momID/femID & maleID| Full ID for each pipefish including information about species (FU = _S. fuscus_, FL = _S. floridae_), trial number, sex (M or F) and fishID | All three |
       | depth_adj | Torso depth adjusted for by the snout-vent length | All three |
       | mated | Indicates a pipefish who has mated at least once (1) or not at all (0) | All three |
       | Sex | Denotes whether a pipefish is male (M) or female (F) | All three |

  - `select_diff_boot_aves.csv`: This dataset includes the selection differentials for snout-vent length averaged across the trials (`average_cal`) for each episode of selection (`episode_sel`) for both sexes (`sex`) and all three species (`species`). Lower and upper 95% confidence intervals are also included (`lower` and `upper`). The means and CIs presented here are generated from the bootstrapping analysis.

### Data Availability
The raw data files containing raw genotyping sequences generated from a 3130xl Genetic Analyzer for _Syngnathus fuscus_ (https://doi.org/10.5281/zenodo.14053624) and _Syngnathus floridae_ (https://doi.org/10.5281/zenodo.10558631) are archived on Zenodo. Original photographs for _S. floridae_, _S. fuscus_, and _Syngnathus scovelli_ that were used to generate the morphometric data (stored in the `all_fem_meso_XX.csv` and the `all_mal_meso_XX.csv` files) are archived on Zenodo (https://doi.org/10.5281/zenodo.22906552). 

The _S. scovelli_ data was pulled from a Dryad repository (http://dx.doi.org/10.5061/dryad.bk03m). Column names were adjusted to match the ones used for the _S. fuscus_ and _S. floridae_ datasets, and additional snout-vent length and torso depth measurements were taken from the photographs, alongside brood pouch area and length for male pipefish, and added to the _S. scovelli_ data (see `all_fem_meso_scovelli.csv` and `all_mal_meso_scovelli.csv`). Length was also converted from cm to mm to match the other species.

## Navigating this repository
The analysis is documented in a series of RMarkdown documents.

### RMarkdown documents
All RMarkdown documents used for the various analyses are located in the directory docs/. Every .Rmd file has been knit into a .md file which is readable on GitHub. Packages and data required to run each .Rmd file can be found at the top of the document alongside version numbers where appropriate. The .Rmd documents do the following things:

  - `selection_analysis_*.Rmd`: Read in the corresponding datasets, calculate reproductive and mating fitness based on the embryo parentage data, calculate summary statistics for males and females, and lastly quantify selection in terms of opportunity for selection, selection differentials, and Bateman's gradient for males and females. The decomposition of selection into pre- and post-copulatory episodes for the opportunity for selection and the selection differentials is also presented in these documents.
    
      - _*_ There will be one .Rmd file for each of three species (`_floridae`, `_fuscus`, `_scovelli`).
        
      - The document `selection_analysis_floridae.Rmd` contains the most detail with the code adapted from this .Rmd for the other two species. In any areas where there are major changes, sufficient detail is provided in the species's .Rmd file.
        
      - To run the .Rmd files for _S. floridae_ and _S. fuscus_ you will need the corresponding `all_meso_XXX.csv` and `EmbryoParentage_XXX.csv` files (stored in /data, see "Data" section above) as well as the `calc_fitness.R` script stored in the directory /R (see below). To run the .Rmd file for _S. scovelli_, all that is needed are the corresponding `all_meso_XXX.csv` files.

  - `cross_species_comp.Rmd`: Reads in datasets generated in the `selection_analysis_*.Rmd` documents and creates several figures which compares selection metrics across the three species. In this document, figures 2, 3, and 4 from the manuscript are generated.
    
      - To run this .Rmd you will need only the datasets outlined at the beginning of the document which includes the `episode_select_data.csv` file, the `XX_fem_bateman.csv` files, the `fem_` and `mal_fitnessXX.csv` files, and the `select_diff_boot_aves.csv` file.

  - `bootstrapping_selection_analysis.Rmd`: Performs the bootstrapping approach to estimate and decompose the opportunity for selection and the selection differentials.
    
      - To run this .Rmd you will need to load the functions from the R scripts `bootstrap_partition_I.R`, `bootstrap_s.R`, `calc_selection_diffs.R`, and `partition_I.R`. You will also need the datasets outlined at the beginning of the document which includes the six `fem_` and `mal_fitnessXX.csv` files located in the data/ directory.

### R files
The directory R/ contains several supporting scripts used in the RMarkdown documents outlined above. They do the following:

  - `calc_fitness.R`: Used to convert the results from the genetic parentage analysis into overall mating and reproductive success for _S. floridae_ and _S. fuscus_.

  - `opp_select_nozeros.R`: Outlines how we partitioned the opportunity for selection when unmated indivudals were only included in the first episode.

  - `bootstrap_partition_I.R`: Performs the bootstrapping procedure to decompose the opportunity for selection.

  - `bootstrap_s.R`: Performs the bootstrapping procedure to estimate and decompose the selection differentials and its confidence intervals.

  - `calc_selection_diffs.R`: Contains the code needed to calculate the selection differentials.

  - `partition_I.R`: Partitions the opportunity for selection into its various components. 

## Contributons
The data presented in this repository was collected by Nicole (Coley) M. Tosto (ORCID: 0000-0001-9858-7046). The code was written by NM Tosto and Sarah P. Flanagan (ORCID: 0000-0002-2226-4213). Please contact NM Tosto with any questions at coley.tosto@canterbury.ac.nz.
