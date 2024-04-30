General Outline:
- GISAID Epi-sets
- Package Versions
- Directory Organization

GISAID Epi-sets:
Training: 
EPI_SET ID: EPI_SET_240430wx
DOI: 10.55876/gis8.240430wx

Test set 01 (2023-on):
EPI_SET ID: EPI_SET_240430em 
DOI: 10.55876/gis8.240430em

Test set 02 (Mexico part 1):
EPI_SET ID: EPI_SET_240430xd 
DOI: 10.55876/gis8.240430xd

Test set 02 (Mexico part 2): 
EPI_SET ID: EPI_SET_240430aq
DOI: 10.55876/gis8.240430aq

Test set 03 (Italy): 
EPI_SET ID: EPI_SET_240430rt
DOI: 10.55876/gis8.240430rt


Package Versions: 
R Version 4.3.2
tidyverse - 2.0.0
janitor - 2.2.0
stringr - 1.5.1
phylotools - 0.2.2
Biostrings - 2.70.3
foreach - 1.5.2
doParallel - 1.0.17
Parallel - 4.3.2
Hmisc - 5.1-2
ggplot2 - 3.5.0
visdat - 0.6.0
bonsai - 0.2.1
lightgbm - 4.3.0
xgboost - 1.7.7.1
ranger - 0.16.0
yardstick - 1.3.1
tidymodels - 1.2.0


A breakdown of the folder directory for this project:
| - /Cleaned_Data/
	| - /aligned_files/
		| -  intermediate aligned files(.csv)
	| -  converted to amino acid aligned files(.csv)
	| -  metadata cleaned files (.csv)
	
| - /Merged_Dataframes/
	| - merged dataframes for inpute into modeling
	| - coefficeint, prediction files generated from modeling fits

| - /Raw_Data/ 
	| - raw .fasta and .tsv files from GISAID.org

| - /Scripts/
	| - Project scripts, 7 in total

| - /Tuning_Metric_Dfs/
	| - saved dataframes from hyperparameter tuning
