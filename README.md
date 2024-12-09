# Measuring Sexual Selection in _Syngnathus_ pipefish

This is a repository for the analysis of sexual selection pressues in three species of pipefish from the genus _Syngnathus_. This includes the dusky pipefish _Syngnathus floridae_, the Northern pipefish _Syngnathus fuscus_, and the Gulf pipefish _Syngnathus scovelli_. The goals of this study are the following:

  1. Investigate the contributions of pre- and post-copulatory sexual selection across all three species.
  2. Generate Bateman's gradient for male and females within each species.
  3. Attempt to identify traits that may be targets of pre- and post-copulatory selection using selection differentials.

## Data
The code refers to data that is found in the data/ directory at the top of the repo. The datasets in this directory were made from the original raw data (see "Data Availability" for how to access the raw data).

  - `all_fem_meso_*.csv` and `all_mal_meso_*.csv`: These files contain the data about the morphometrics of all male and female pipefish and various information relating to the trials they were in.
     
     - For `S. scovelli` the columns **svl**, **depth**, **bp_area**, and **bp_length** were not included in the orignal dataset. The original photos were obtained and those additional measurements were gathered.   

  - `EmbryoParentage_*.csv`: This is the file that is used to calculate reproductive fitness and mating success for male and female pipefish.
    
      - _*_ There will be one male and female .csv file for only two of the three species (`_floridae` and `_fuscus`). _S. scovelli_ does not have a file in this format as the publically accesible dataset already contains information about reproductive fitness and mating success. 

### Data Availability
The raw data files containing raw genotyping sequences generated from a 3130xl Genetic Analyzer for _Synganthus fuscus_ are archived on Zenodo (https://zenodo.org/records/14053624).

## Navigating this repository
The analysis is documented in a series of RMarkdown documents.

### RMarkdown documents
All Rmarkdown documents used for the various analyses are located in the directory docs/. They do the following things:

  - `selection_analysis_*.Rmd`: Read in the corresponding datasets, calculate reproductive and mating fitness based on the embryo parentage data, calculate summary statistics for males and females, and lastly quantify selection in terms of opportunity for selection, selection differentials, and Bateman's gradient for males and females.
    
      - _*_ There will be one .Rmd file for each of three species (`_floridae`, `_fuscus`, `_scovelli`).
        
      - The document `selection_analysis_floridae.Rmd` contains the most detail with the code adapted from this .Rmd for the other two species. In any areas where there are major changes, sufficient detail is provided in the species's .Rmd file.

### R files
The directory R/ contains the script `calc_fitness.R` that was used to convert the results from the genetic parentage analysis into overall mating and reproductive success.

## Contributors
The contributors to this repository are Coley Tosto and Sarah Flanagan. Please contact Coley with any questions at coley.tosto@canterbury.ac.nz.
