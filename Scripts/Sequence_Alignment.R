#Grabbing the sequences for only the accession IDs i am going to use
#in both the training and test sets
library(stringr)
library(phylotools)
library(Biostrings)
library(tidyverse)

default_wd <- getwd()

train <- read.csv("Cleaned_Data/training_metadata_cleaned.csv")
test <- read.csv("Cleaned_Data/testing_metadata_cleaned.csv")

train_acc <- train |> 
  pull(accession_id)

test_acc <- test |> 
  pull(accession_id)

data_d <- str_c(default_wd, "/Raw_Data")
setwd(data_d)

fasta_list <- tibble(file_name = list.files())

fasta_list <- fasta_list |> 
  filter(str_detect(file_name, ".fasta")) |> 
  filter(!str_detect(file_name, "reference"))

reference_seq <- read.fasta(file = "EPI_ISL_402124-S_reference_genome.fasta")




read_fasta_files <- function(file_name) {
  fasta_tmp <- read.fasta(file = file_name_tmp, clean_name = TRUE)
  
  ifelse(str_detect(file_name_tmp, "test"), 
         fasta_mod <- fasta_tmp |> 
           mutate(accession_id = str_extract(seq.name, "EPI_ISL_[0-9]+")) |> 
           filter(accession_id %in% test_acc), #yes this is the test data
         
         fasta_mod <- fasta_tmp |> 
           mutate(accession_id = str_extract(seq.name, "EPI_ISL_[0-9]+")) |> 
           filter(accession_id %in% train_acc)) #no this is the train data
  
  return(fasta_mod) }

fasta_file_names <- list()


for (i in 1:nrow(fasta_list)) {
  file_name_tmp <- fasta_list$file_name[i]
  file_name_fixed <- str_split(file_name_tmp, ".fa")[[1]][1]
  
  fasta_file_names[[i]] <- file_name_fixed
  
  assign(file_name_fixed, read_fasta_files(file_name_tmp))
}

file_name_aligned <- list()

#splitting the training files up into 2 halves each to make computing a bit easier
fasta_file_names

train_to_AUG_01 <- gisaid_hcov_19_to_AUG_01[0:2500, ]
train_to_AUG_02 <- gisaid_hcov_19_to_AUG_01[2501:nrow(gisaid_hcov_19_to_AUG_01), ]
train_AUG_to_NOV_01 <- gisaid_hcov_19_AUG_to_NOV_01[0:2200, ]
train_AUG_to_NOV_02 <- gisaid_hcov_19_AUG_to_NOV_01[2201:nrow(gisaid_hcov_19_AUG_to_NOV_01), ]
train_NOV_to_DEC_01 <- gisaid_hcov_19_NOV_to_DEC_31[0:2000, ]
train_NOV_to_DEC_02 <- gisaid_hcov_19_NOV_to_DEC_31[2001:nrow(gisaid_hcov_19_NOV_to_DEC_31), ]

fasta_file_names <- list("train_to_AUG_01", "train_to_AUG_02", "train_AUG_to_NOV_01", 
                         "train_AUG_to_NOV_02", "train_NOV_to_DEC_01", "train_NOV_to_DEC_02",
                         fasta_file_names[[4]], fasta_file_names[[5]], fasta_file_names[[6]],
                         fasta_file_names[[7]])

aligning_fn <- function(i) {
  aligning <- get(fasta_file_names[[i]])
  aligned <- aligning |> 
    #head(n) |> 
    select(accession_id, seq.text) |> 
    mutate(aligned_seq = as.character(
      Biostrings::pattern(pairwiseAlignment(
        pattern = seq.text, subject = reference_seq$seq.text, type = "local")))) |> 
    select(accession_id, aligned_seq)
  return(aligned)}

aligned_d <- str_c(default_wd, "/Cleaned_Data/aligned_files")
setwd(aligned_d)

#library(parallel)
library(foreach)
library(doParallel)
detectCores() #10 cores
registerDoParallel(cores = 10)

file_name_aligned = list()

#off to the races. will take a while
foreach(i = 1:length(fasta_file_names)) %dopar% {
  file_name_aligned[[i]] <- str_c(fasta_file_names[[i]], "_aligned")
  file_name_csv <- str_c(file_name_aligned[[i]], ".csv")
  temp_aligned_obj <- aligning_fn(i)
  write.csv(temp_aligned_obj, file = file_name_csv)
  }


