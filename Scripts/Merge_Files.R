#Merging the meta data files & the aligned amino acid sequence files
default_wd <- getwd()
cleaned_wd <- str_c(default_wd, "/Cleaned_Data/")
setwd(cleaned_wd)

library(tidyverse)

files <- as.list(list.files())

for (k in 2:length(files)) {
  files[[k]]
  df <- read_csv(files[[k]], col_types = cols(.default = "c")) |> 
    as_tibble()
  file_name <- str_split_i(files[[k]], ".csv", 1)
  
  assign(file_name, df)
}

seq_files <- list("test_set_01_2023_on_aligned_aa", "test_set_02_mexico_part1_aligned_aa", 
                     "test_set_02_mexico_part2_aligned_aa", "test_set_03_italy_aligned_aa",
                     "train_AUG_to_NOV_01_aligned_aa", "train_AUG_to_NOV_02_aligned_aa", 
                     "train_NOV_to_DEC_01_aligned_aa", "train_NOV_to_DEC_02_aligned_aa", 
                     "train_to_AUG_01_aligned_aa", "train_to_AUG_02_aligned_aa")

stop_codons_trimming <- function(df) {
  
  last_col <- df |> colnames() |> tail(1)
  
  df_mod <- df |> 
    select(-1) |> 
    unite("aa_united", aa_1:{{ last_col }}, sep = "", remove = TRUE, na.rm = TRUE) |> 
    mutate(aa_united = str_split_i(aa_united, "STOP", 1),
           length = nchar(aa_united)) |> 
    relocate(length, .before = aa_united) |> 
    separate_wider_position(cols = aa_united, 
                            widths = c(aa_ = rep(1, 1282)), 
                            too_few = "align_start")
  return(df_mod)
}


for (k in 1:10) {
  temp_df <- get(seq_files[[k]]) |> 
    stop_codons_trimming()
  
  assign(seq_files[[k]], temp_df)
  
}


for (i in 1:10) {
  get(seq_files[[i]]) %>% 
    keep(~all(is.na(.x))) %>% 
    names() |> 
    print()
}
#1276 and on are all NA for all columns


for (i in 1:10) {
  df_mod <- get(seq_files[[i]]) %>% 
    select(-c("aa_1276", "aa_1277", "aa_1278", "aa_1279", "aa_1280", "aa_1281", "aa_1282"))
  
  assign(seq_files[[i]], df_mod)
}

for (i in 1:4) {
  df_mod <- get(seq_files[[i]]) |> 
    mutate(test_set_id = parse_number(seq_files[[i]])) |> 
    relocate(test_set_id, .after = accession_id) |> 
    dplyr::rename(protein_length = 3) #length
  
  assign(seq_files[[i]], df_mod)
}

for (i in 5:10) {
  df_mod <- get(seq_files[[i]]) |> 
    dplyr::rename(protein_length = 2) #length
  
  assign(seq_files[[i]], df_mod)
}


testing_seq_data <- 
  rbind(test_set_01_2023_on_aligned_aa, test_set_02_mexico_part1_aligned_aa, 
        test_set_02_mexico_part2_aligned_aa, test_set_03_italy_aligned_aa)

training_seq_data <- 
  rbind(train_AUG_to_NOV_01_aligned_aa, train_AUG_to_NOV_02_aligned_aa, train_NOV_to_DEC_01_aligned_aa,
        train_NOV_to_DEC_02_aligned_aa, train_to_AUG_01_aligned_aa, train_to_AUG_02_aligned_aa)

#I suggest removing the intermediate data frames from the environment now

#joining time :)
#i wanted the meta data first so right join
testing_df <- right_join(testing_metadata_cleaned, 
                         testing_seq_data, 
                         by = "accession_id") |> 
  filter(protein_length > 800) #800 is a bit arbitrary but needed

training_df <- right_join(training_metadata_cleaned, 
                          training_seq_data,
                          by = "accession_id") |> 
  filter(protein_length > 800) #800 is a bit arbitrary but needed

library(visdat)
training_df |> 
  select(c(1280:1286)) |> 
  vis_miss(warn_large_data = FALSE) #last two columns as essentially NAs

testing_df |> 
  select(c(1280:1288)) |> 
  vis_miss(warn_large_data = FALSE) #last two columns are essentially NAs

training_df <- training_df |> 
  select(-c(1285, 1286))

testing_df <- testing_df |> 
  select(-c(1287, 1288))

#final cleaning, to prep for modeling
#toss out columns I wont be using
#disease_status can be left alone,
#gender is to be converted to 0 or 1
#and I am going to try out tokenization instead of dummy variable encoding...
#and thus anything with 'X' as the code is being converted to NA
training_df <- training_df |> 
  select(-c(1, 6, 7, 8, 9, 10, 11)) |> 
  mutate(patient_age = as.numeric(patient_age),
         gender = case_when(
           gender == "Male" ~ 0,
           gender == "Female" ~ 1)) |> 
  pivot_longer(cols = -c(1:4), 
               names_to = "aa_pos",
               values_to = "aa_code") |> 
  mutate(aa_code = case_when(
    aa_code == "X" ~ NA,
    .default = aa_code)) |> 
  pivot_wider(names_from = "aa_pos",
              values_from = "aa_code")


testing_df <- testing_df |> 
  select(-c(1, 6, 7, 8, 9, 10, 11, 13)) |> 
  mutate(patient_age = as.numeric(patient_age),
         gender = case_when(
           gender == "Male" ~ 0,
           gender == "Female" ~ 1)) |> 
  pivot_longer(cols = -c(1:5), 
               names_to = "aa_pos",
               values_to = "aa_code") |> 
  mutate(aa_code = case_when(
    aa_code == "X" ~ NA,
    .default = aa_code)) |> 
  pivot_wider(names_from = "aa_pos",
              values_from = "aa_code")

merged_dir <- str_c(default_wd, "/Merged_Dataframes/")
setwd(merged_dir)

write_csv(testing_df, "Testing_Data_Merged.csv")
write_csv(training_df, "Training_Data_Merged.csv")




