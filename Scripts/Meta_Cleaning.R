############# String repair on gender and age variables
# Data cleaning
library(tidyverse)
library(stringr)
default_wd <- getwd()
data_d <- str_c(default_wd, "/Raw_Data")
setwd(data_d)

test_df <- read_csv("merged_test_filtered_uncleaned.csv", col_names = TRUE)
train_df <- read_csv("merged_train_filtered_uncleaned.csv", col_names = TRUE)

#DOB instead
train_df |> 
  filter(str_detect(patient_age, "[0-9]{2}[-][A-Z]")) |> 
  glimpse() # 1 case
# i am going to exclude this as it would make this person 121 years old. I think this info is entered wrong
#but I wont assume what it should be i will just exclude it

#age and gender are wrong
train_df |> 
  filter(str_detect(patient_age, "[MF][0-9]")) |> 
  select(patient_age, gender) # 3 cases
#I am actually going to exclude any like this because I feel this data is conflicting and cannot be cleared up
#and there are none like this in the test df

test_df |> 
  select(patient_age) |> 
  table()
#<1 I think should get caught
#90+ should get caught
#ranges, months should be handled
#"ambulatory" "hospitalized", "ND" are weird

test_df |> 
  filter(patient_age == "Ambulatory" | patient_age == "Hospitalized"
         | patient_age == "ND" | patient_age == "unknow") |> 
  select(patient_age, gender, patient_status, passage) |> 
  print(n = 28)
#most seem to be wrongly entered into the database, I am going to exclude them. 
#since this doesnt appear to be a simple two columns were swapped, but rather entire frameshifts
#they will be excluded
#the unknow & ND ones will also be excluded as they are missing too much info


train_df_filtered <- train_df |> 
  filter(!str_detect(patient_age, "[0-9]{2}[-][A-Z]")) |> 
  filter(!str_detect(patient_age, "[MF][0-9]"))
  
test_df_filtered <- test_df |> 
  filter(patient_age != "Ambulatory" & patient_age != "Hospitalized" 
         & patient_age != "ND" & patient_age != "unknow")

age_correction <- function(df) {
  df |> 
    mutate(patient_age_corrected = case_when(
      patient_age == "Female" | patient_age == "Male" ~ gender,
      patient_age == "?" ~ NA,
      str_detect(patient_age, "[0-9]{2}[-][0-9]{2}") ~ str_extract(patient_age, "[0-9]+"),
      str_detect(patient_age, "[mM][oO]") & as.integer(str_extract(patient_age, "[0-9]+")) <= 6 ~ "0", #younger than 7 months
      str_detect(patient_age, "[mM][oO]") & as.integer(str_extract(patient_age, "[0-9]+")) > 6 
      & as.integer(str_extract(patient_age, "[0-9]+")) <= 15 ~ "1", #btwn 7 & 15 months
      str_detect(patient_age, "[mM][oO]") & as.integer(str_extract(patient_age, "[0-9]+")) > 15 ~ "2", #older than 15 months
      str_detect(patient_age, "[0-9][-][0-9]") ~ str_extract(patient_age, "[0-9]"), 
      str_detect(patient_age, "days") ~ "1", #266 days old calling it as 1 year
      str_detect(patient_age, "[0-9]{2}[+]") ~ "90",
      str_detect(patient_age, "[0-9]{2}[mMfF]") ~ str_extract(patient_age, "[0-9]+"),
      str_detect(patient_age, "[0-9]+ [-]") ~ str_extract(patient_age, "[0-9]+"),
      str_detect(patient_age, "[<][0-9]+") ~ str_extract(patient_age, "[0-9]+"),
      str_detect(patient_age, "[0-9]+[aA]") ~ str_extract(patient_age, "[0-9]+"), # no clue what A stands for
      .default = patient_age),
      
      gender = case_when(
        patient_age == "Female" | patient_age == "Male" ~ patient_age,
        .default = gender)) |> 
    mutate(patient_age_corrected = as.integer(patient_age_corrected))
}

train_df_corr <- train_df_filtered |> 
  age_correction()

test_df_corr <- test_df_filtered |> 
  age_correction()

train_df_corr |> 
  filter(is.na(patient_age_corrected)) |> 
  glimpse() #too much missing info, going to remove from training

test_df_corr |> 
  filter(is.na(patient_age_corrected)) |> 
  glimpse() #for the first, i am going to delete. 
#for the second, I think I can save it by assuming they meant 78

train_df_corr <- train_df_corr |> 
  filter(!is.na(patient_age_corrected))

test_df_corr <- test_df_corr |> 
  mutate(patient_age_corrected = case_when(
    accession_id == "EPI_ISL_2622160" ~ 78,
    .default = patient_age_corrected)) |> 
  filter(!is.na(patient_age_corrected))

#dropping variables I dont need
train_df_corr <- train_df_corr |> 
  select(accession_id, gender, patient_age_corrected, location, 
         patient_status, additional_host_information, lineage, clade) |> 
  rename(addnl_host_info = additional_host_information,
         patient_age = patient_age_corrected) |> 
  relocate(location, .after = clade) |> 
  relocate(patient_status, .before = gender)

test_df_corr <- test_df_corr |> 
  select(accession_id, patient_status, gender, patient_age_corrected, country, 
         additional_host_information, lineage, clade, continent) |> 
  rename(addnl_host_info = additional_host_information,
         patient_age = patient_age_corrected)

test_df_corr |> 
  select(gender) |> 
  table() #"missing"

test_df_corr |> 
  filter(gender == "Missing")

test_df_corr <- test_df_corr |> 
  filter(gender != "Missing")

############# Binning the patient_status variable

train_list <- train_df_corr |> 
  select(patient_status) |> 
  mutate(patient_status = tolower(patient_status)) |> 
  unique() |> 
  mutate(n = 1:n())

test_list <- test_df_corr |> 
  select(patient_status) |> 
  mutate(patient_status = tolower(patient_status)) |> 
  unique() |> 
  mutate(n = 1:n())

joined <- full_join(test_list, train_list, by = "patient_status")
#limited overlap..... TT_TT

#keeping to three categories as like the paper- unknown, mild, and severe
unknown_list <- c("not provided", "unown", "missing", "unknwon", "ou", "no", "yes", 
                  "not vaccinated", "vaccinated", "unkown", "unknow", "nasal swab", "Nasal")
mild_list <- c("asymptomatic", "symptomatic", "live", "outpatient", "non hospitalized", "non-hospitalized",
               "alive", "mild; recovered", "home isolation", "mild", "active", 
               "active, non-sentinel-surveillance (hospital)", "nurse house", "paucisymptomatic",
               "recovered after home treatment", "asymptomatic - ambulatory", "mild disease", 
               "paucisymptomatyc", "outpatient mild disease", "asymptomatyc", "oupatient", "moderate", 
               "symptomatic/mild illness/live", "asymtomatic", "mild infection", "outpatient, mild disease", 
               "symptomatic – mild disease", "not hospitalized", "live, recovered", "screening", "routine",
               "quarantine isolation")
severe_list <- c("hospitalized", "ambulatory", "recovered", "released", "severe acute respiratory syndrome", 
                 "hospital cases", "deceased", "hospitalised", "emergency room", "deceased (brought in dead)",
                 "hospitalized, live", "inpatient", "severe", "live, ambulatory care", "symptomatic ambulatory",
                 "symptomatic and ambulatory", "symptomatic - ambulatory", "fatal", 
                 "symptomatic - hospitalized", "severe, hospitalized", "hospitalized, deceased", 
                 "deceased, severe, hospitalized", "hospitalized, severe", "dead")

bin_patient_status <- function(df) {
  df |> 
    mutate(patient_status = tolower(patient_status)) |> 
    mutate(disease_status = case_when(
      patient_status %in% unknown_list ~ "unknown",
      patient_status %in% mild_list ~ "mild", 
      patient_status %in% severe_list ~ "severe")) |> 
    relocate(disease_status, .before = gender) |> 
    relocate(patient_status, addnl_host_info, .after = clade) }

train_df_corr <- train_df_corr |> 
  bin_patient_status() |> 
  filter(disease_status != "unknown")

test_df_corr <- test_df_corr |> 
  bin_patient_status() |> 
  filter(disease_status != "unknown")

setwd(str_c(default_wd, "/Cleaned_Data/"))

write.csv(train_df_corr, "training_metadata_cleaned.csv")
write.csv(test_df_corr, "testing_metadata_cleaned.csv")
