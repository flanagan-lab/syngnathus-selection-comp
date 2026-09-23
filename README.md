# Measuring Sexual Selection in _Syngnathus_ pipefish

This is a repository for the analysis of sexual selection pressues in three species of pipefish from the genus _Syngnathus_. This includes the dusky pipefish _Syngnathus floridae_, the Northern pipefish _Syngnathus fuscus_, and the Gulf pipefish _Syngnathus scovelli_. The goals of this study are the following:

  1. Investigate the contributions of pre- and post-copulatory sexual selection across all three species.
  2. Generate Bateman's gradient for male and females within each species.
  3. Attempt to identify traits that may be targets of pre- and post-copulatory selection using selection differentials.

These analyses support the findings in the manuscript "Comparing mating systems and sexual selection pressures across three congeneric species of pipefish that span the continuum of sexual dimorphism", currently submitted to a journal for consideration.

## Data
The RMarkdown documents refer to data that is found in the data/ directory at the top of the repo. The datasets in this directory were made from the original raw data (see "Data Availability" for how to access the raw data).

  - `all_fem_meso_*.csv` and `all_mal_meso_*.csv`: These files contain the data about the morphometrics of all male and female pipefish and various information relating to the trials they were in.
     
     - For the `_scovelli` datasets, the columns **svl**, **depth**, **bp_area**, and **bp_length** were not included in the original dataset. The original photos were obtained and those additional measurements were gathered. See the "Data Availability" section for more detail about the changes made to the publicly accessible _S. Scovelli_ dataset.   

  - `EmbryoParentage_*.csv`: This is the file that is used to calculate reproductive fitness and mating success for male and female pipefish.
    
      - _*_ There will be one male and female .csv file for only two of the three species (`_floridae` and `_fuscus`). _S. scovelli_ does not have a file in this format as the publically accesible dataset already contains information about reproductive fitness and mating success. 

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
