############# Training data merging files
# Data cleaning
library(tidyverse)
library(janitor)
default_wd <- getwd()
data_d <- str_c(default_wd, "/Raw_Data")
setwd(data_d)
Start_Aug_01 <- read_delim("gisaid_hcov-19_to_AUG_01.tsv", show_col_types = FALSE)
Aug_02_Nov_01 <- read_delim("gisaid_hcov-19_AUG_to_NOV_01.tsv", show_col_types = FALSE)
Nov_02_Dec_31 <- read_delim("gisaid_hcov-19_NOV_to_DEC_31.tsv", show_col_types = FALSE)
df <- rbind(Start_Aug_01, Aug_02_Nov_01, Nov_02_Dec_31)

df_train <- df |> 
  clean_names() |> 
  filter(host == "Human")

df_train |> 
  glimpse()

df_train |>
  filter(!is.na(last_vaccinated)) #only 891 rows, so cannot be used

df_train |> 
  filter(gender == "Male" | gender == "Female") |> 
  filter(patient_age != "unknown" | patient_age != "N/A") #13,778 entries

#Since making sure there are entries for gender, and for age leaves 13k entries
#I am going to go forward with this filtering.... however....

table(df_train$gender)
#why do some of them have numbers for gender?

df_train |> 
  filter(gender != "Male" & gender != "Female" & gender != "unknown") |> 
  head() |> 
  select(patient_age, gender)

df_train <- df_train |> 
  filter(gender != "unknown") |> 
  filter(patient_age != "unknown")

#so some of the data was entered wrongly, it has age and gender swapped. 
#I am going to write the merged file for now, as a checkpoint
#before mutating or filtering any further

write.csv(df_train, file = "merged_train_filtered_uncleaned.csv")

############# Test data merging files

#"2023 on"- a little over 2000 records for all regions
#"Mexico" - similar date range- 2021-01-01 to 2022-01-01 (roughly)
#         - randomly sampling 2000 records
#         - from 2021-01-01 to 2021-07-21 & 1000 from 2021-08-10 to 2022-01-01
#"Italy"  - similar date range
#         - 2000 randomly sampled from 2021-04-01 to 2022-01-01

test_set_01 <- read_delim("test_set_01_2023_on.tsv", show_col_types = FALSE) |> clean_names()
test_set_02_01 <- read_delim("test_set_02_mexico_part1.tsv", show_col_types = FALSE) |> clean_names()
test_set_02_02 <- read_delim("test_set_02_mexico_part2.tsv", show_col_types = FALSE) |> clean_names()
test_set_03 <- read_delim("test_set_03_italy.tsv", show_col_types = FALSE) |> clean_names()

#all of them have some missing info for gender age records, but not many
test_set_03 |> 
  filter(gender != "Male" & gender != "Female" & gender != "unknown") |> 
  head() |> 
  select(patient_age, gender)
#there are at least a few that have patient age and gender swapped like with the training set
#so I will do my best to save those but anything with gender unknown is going to be eliminated, 
#like with the training

df_test <- rbind(test_set_01, test_set_02_01, test_set_02_02, test_set_03)
#after filtering & sampling, will be left with a 6,000 record(ish) test set
#which can be sub-split into 3 different test sets

df_test |> select(gender) |> table()

df_test <- df_test |> 
  filter(gender != "unknown") |> 
  filter(patient_age != "unknown")

df_test <- df_test |> 
  separate_wider_delim(cols = collection_date, delim = "-", names = c("year", "month", "day")) |> 
  separate_wider_delim(cols = location, delim = " / ", names = c("continent", "country", "extra"),
                       too_few = "align_start", too_many = "merge") |> 
  mutate(test_set_id = case_when(year == 2023 | year == 2024 ~ "1",
                                 year != 2023 & year != 2024 & country == "Mexico" ~ "2",
                                 year != 2023 & year != 2024 & country == "Italy" ~ "3"))

set.seed(328)

mexico_list <- df_test |> 
  filter(test_set_id == 2) |> 
  slice_sample(n = 2000) |> 
  select(accession_id) |> 
  pull(accession_id)

italy_list <- df_test |> 
  filter(test_set_id == 3) |> 
  slice_sample(n = 2000) |> 
  select(accession_id) |> 
  pull(accession_id)

df_test_filtered <- df_test |> 
  filter(test_set_id == 01 | accession_id %in% mexico_list | accession_id %in% italy_list)

write.csv(df_test_filtered, file = "merged_test_filtered_uncleaned.csv")



