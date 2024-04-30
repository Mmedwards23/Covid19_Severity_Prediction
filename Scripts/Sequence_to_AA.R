#Converting the aligned sequences into amino acid residues
library(Hmisc)
library(tidyverse)
library(ggplot2)

default_wd <- getwd()
aligned_d <- str_c(default_wd, "/Cleaned_Data/aligned_files")
setwd(aligned_d)

aligned_list <- list(list.files())
df_list <- list()
#read the csvs back in
for (i in 1:length(aligned_list[[1]])) {
  df_list[[i]] <- str_split(aligned_list[[1]][i], ".csv")[[1]][1]
  temp <- as_tibble(read.csv(file = aligned_list[[1]][i])) |> select(-X)
  assign(df_list[[i]], temp)
}

#checking the ranges of sequence length
lengths <- tibble("min" = 0, "max" = 0)
for (i in 1:10) {
  df <- get(df_list[[i]]) |> 
    mutate(seq_length = nchar(aligned_seq))
  
  lengths[i, 1] <- min(df$seq_length)
  lengths[i, 2] <- max(df$seq_length)
}
lengths #something strange going on within the first dataset, we are missing almost 1000 basepairs

test_set_01_2023_on_aligned |> 
  mutate(seq_length = nchar(aligned_seq)) |> 
  filter(seq_length < 3000)
#these 8 need to be removed, there is no reason why their sequence length should be so truncated.... 
#im exporting their accession IDs to look into separately

setwd(default_wd)
truncated_df <- test_set_01_2023_on_aligned |> 
  mutate(seq_length = nchar(aligned_seq)) |> 
  filter(seq_length < 3000) |> 
  select(accession_id)

truncated <- as.list(truncated_df |> pull(accession_id))

truncated_df |> 
  write_csv(file = "Oddly_Truncated_Sequences.csv")

test_set_01_2023_on_aligned <- test_set_01_2023_on_aligned |> 
  filter(accession_id %nin% truncated)

#converting the BP to AA
#first I need to get rid of any deletions
#essentially here i am cutting out any of the -, which denote a deleted BP(or more than one)
#need to have the upper limit be 3835 BP, which ~ 1279 AA
for (i in 1:10) {
  df_mod <- get(df_list[[i]]) |> 
    mutate(aligned_seq = str_replace_all(aligned_seq, "[-]+", "")) |> 
    separate_wider_position(cols = aligned_seq, 
                            widths = c(aa_ = rep(3, 1279)), 
                            too_few = "align_start")
  
  assign(df_list[[i]], df_mod)
}

#if it doesnt start with ATG(the start codon), it is being removed
for (i in 1:10) {
  df_mod <- get(df_list[[i]]) |> 
    filter(aa_1 == "ATG") |> 
    pivot_longer(cols = -c(accession_id), 
                 names_to = "aa_position", 
                 values_to = "aa_seq")
  
  assign(df_list[[i]], df_mod)
}


basepairs_to_aa <- function(df) {
  df_mod <- df |> 
    mutate(aa_code = case_when(
      str_detect(aa_seq, "TT[TCY]") ~ "F",
      str_detect(aa_seq, "TT[AGR]") ~ "L",
      str_detect(aa_seq, "TC[ACGTRYKMSWBDHVN]") | str_detect(aa_seq, "AG[TCY]") ~ "S", 
      str_detect(aa_seq, "TA[TCY]") ~ "Y",
      str_detect(aa_seq, "TA[AGR]") | aa_seq == "TGA" ~ "STOP",
      str_detect(aa_seq, "TG[TCY]") ~ "C",
      aa_seq == "TGG" ~ "W",
      str_detect(aa_seq, "CT[ACGTRYKMSWBDHVN]") ~ "L", 
      str_detect(aa_seq, "CC[ACGTRYKMSWBDHVN]") ~ "P", 
      str_detect(aa_seq, "CA[TCY]") ~ "H",
      str_detect(aa_seq, "CA[AGR]") ~ "Q",
      str_detect(aa_seq, "CG[ACGTRYKMSWBDHVN]") | str_detect(aa_seq, "AG[AGR]") ~ "R", 
      str_detect(aa_seq, "AT[TCYA]") ~ "I",
      aa_seq == "ATG" ~ "M", 
      str_detect(aa_seq, "AC[ACGTRYKMSWBDHVN]") ~ "T", 
      str_detect(aa_seq, "AA[TCY]") ~ "N",
      str_detect(aa_seq, "AA[AGR]") ~ "K",
      str_detect(aa_seq, "GT[ACGTRYKMSWBDHVN]") ~ "V", 
      str_detect(aa_seq, "GC[ACGTRYKMSWBDHVN]") ~ "A", 
      str_detect(aa_seq, "GA[TCY]") ~ "D",
      str_detect(aa_seq, "GA[AGR]") ~ "E",
      str_detect(aa_seq, "GG[ACGTRYKMSWBDHVN]") ~ "G",
    
      #and if it doesnt fall into the above cases nicely... 
      str_detect(aa_seq, "R") | str_detect(aa_seq, "Y") | str_detect(aa_seq, "K") | 
        str_detect(aa_seq, "M") | str_detect(aa_seq, "S") |  str_detect(aa_seq, "W") | 
        str_detect(aa_seq, "B") | str_detect(aa_seq, "D") | str_detect(aa_seq, "H") | 
        str_detect(aa_seq, "V") | str_detect(aa_seq, "N")  ~ "X", 
      #indicating ambiguous nucleotide call cant be clarified
      #not much else I can do with this
    
      str_detect(aa_seq, "[ACGTRYKMSWBDHVN]{2}") ~ "", #if it is incomplete, and I cant use
      #the wobble effect to my advantage, then set it to blank
    
      .default = aa_seq #if anything else is going on, do not convert
    ))
  return(df_mod)
  }

cleaned_wd <- str_c(default_wd, "/Cleaned_Data/")
setwd(cleaned_wd)

library(foreach)
library(doParallel)

registerDoParallel(cores = 10)
foreach(i = 1:10) %dopar% {
  df_mod <- get(df_list[[i]]) |> 
    basepairs_to_aa() |> 
    select(-aa_seq) |> 
    filter(!is.na(aa_code)) |> #drop any rows with NA values
    #NAs will be introduced in the pivot wider as are needed
    pivot_wider(names_from = aa_position,
                values_from = aa_code)
  
  write.csv(df_mod, str_c(df_list[[i]], "_aa.csv"))
}
