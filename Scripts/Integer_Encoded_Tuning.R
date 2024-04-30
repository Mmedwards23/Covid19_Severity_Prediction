library(bonsai) #lightGBM
library(lightgbm)
library(xgboost)
library(ranger)
library(yardstick)
library(tidyverse)
library(tidymodels)
tidymodels_prefer()

default_wd <- getwd()

train_df <- read_csv("Merged_Dataframes/Training_Data_Merged.csv",
                     col_types = cols(.default = "c", gender = "i", patient_age = "i"))
test_df <- read_csv("Merged_Dataframes/Testing_Data_Merged.csv",
                    col_types = cols(.default = "c", gender = "i", patient_age = "i"))

train_df <- train_df |> 
  pivot_longer(cols = -c(1:4),
               names_to = "aa_pos",
               values_to = "aa_code") |> 
  mutate(aa_code = case_when(
    aa_code == "X" ~ "zunknown",
    aa_code = is.na(aa_code) ~ "zunknown",
    .default = aa_code)) |> 
  mutate(aa_code = as.integer(as.factor(aa_code))) |> 
  mutate(aa_code = case_when(
    aa_code == 21 ~ 0, #change the NA values to 0 instead
    .default = aa_code)) |> 
  pivot_wider(names_from = "aa_pos",
              values_from = "aa_code") |> 
  select(-1) #remove accession ID

test_df <- test_df |> 
  pivot_longer(cols = -c(1:5),
               names_to = "aa_pos",
               values_to = "aa_code") |> 
  mutate(aa_code = case_when(
    aa_code == "X" ~ "zunknown",
    aa_code = is.na(aa_code) ~ "zunknown",
    .default = aa_code)) |> 
  mutate(aa_code = as.integer(as.factor(aa_code))) |> 
  mutate(aa_code = case_when(
    aa_code == 21 ~ 0, #change the NA values to 0 instead
    .default = aa_code)) |> 
  pivot_wider(names_from = "aa_pos",
              values_from = "aa_code") |> 
  select(-1) #remove accession ID

write_csv(train_df, "Merged_Dataframes/integer_encoded_training_data.csv")
write_csv(test_df, "Merged_Dataframes/integer_encoded_testing_data.csv")

for (i in 1:3) {
  name <- str_c("test_set_", i)
  
  test_df_filt <- test_df |> 
    filter(test_set_id == i) |> 
    select(-test_set_id)
  
  assign(name, test_df_filt)
  
}

recipe_default <- recipe(disease_status ~ ., 
                         train_df) |>
  step_zv(all_predictors()) |>  # removes covariates with no variance
  step_other(all_nominal_predictors(), threshold = 0.03) |>  #anything with >3% occurrence is put under 'other'
  step_normalize(patient_age) |> #normalize age
  step_novel(all_nominal_predictors()) #allow for new factors in test data

#I have a lot of data, and while I could split up the training data into partitions or a validation set
#to use on tuning the parameters
#so i will mc split for 3, and then use the test_set_03 later after I have my 'results' 
#to make sure it doesnt look too wild

log <- logistic_reg() |>
  set_engine("glm") |> 
  set_mode("classification")

wf_log <- workflow() |>
  add_model(log) |>
  add_recipe(recipe_default)

#lasso
lasso <- logistic_reg(penalty = tune(), 
                      mixture = tune()) |> #elastic net with 0.5.....
  set_engine("glmnet") |> 
  set_mode("classification")

wf_lasso <- workflow() |> 
  add_model(lasso) |> 
  add_recipe(recipe_default)

#random forest
rf <- rand_forest(trees = tune(), min_n = 2) |> 
  #could tune mtry, trees, and min_n
  set_mode("classification") |>
  set_engine("ranger")

wf_rf <- workflow() |> 
  add_model(rf) |> 
  add_recipe(recipe_default)

#XGBoost
xgb <- boost_tree(trees = tune(), tree_depth = tune(), learn_rate = tune()) |> 
  set_mode("classification") |> 
  set_engine("xgboost")

wf_xgb <- workflow() |> 
  add_model(xgb) |> 
  add_recipe(recipe_default)

#LightGBM
lightgbm <- boost_tree(trees = tune(), tree_depth = tune(), learn_rate = tune()) |> 
  set_mode("classification") |> 
  set_engine("lightgbm")

wf_lightgbm <- workflow() |> 
  add_model(lightgbm) |> 
  add_recipe(recipe_default)

set.seed(3904)
covid_metrics <- metric_set(yardstick::roc_auc, yardstick::accuracy)

comparison_split <- train_df |>
  mc_cv(times = 3)

#models with tuning
# - elastic net/lasso: penalty, elastic net -ness
# - random forest: trees
# - xgboost: trees, tree_depth, learn_rate
# - lightgbm: trees, tree_depth, learn_rate

#elastic net/lasso: penalty
logistic_grid <- crossing(penalty = c(0.001, 0.01, 0.1),
                             mixture = seq(0, 1, by = 0.1))

library(doParallel)
cl <- makePSOCKcluster(6) 
registerDoParallel(cl) 

logistic_tuning <- wf_lasso  |>
  tune_grid(
    resamples = comparison_split,
    metrics = covid_metrics,
    grid = logistic_grid)

stopCluster(cl)

setwd(str_c(default_wd, "/Tuning_Metric_Dfs/"))

write_csv(logistic_tuning |> collect_metrics(), file = "logistic_tuning.csv")

#with this, accuracy favors mixture 0.6 or 0.2 with the penalty of 0.001
#roc_auc preferred any of the mixtures with a penalty of 0.001. so I am going to use
#mixture - 0.6 & penalty = 0.001
logistic_tuning |> 
  collect_metrics() |> 
  slice_max(mean)
#although the std_error was smaller with penalty of 0.001, I wanted to try elastic net which would be
#mixture of 0.5, as opposed to simply lasso

# - random forest: trees
rf_grid <- crossing(trees = seq(100, 2000, by = 50))

set.seed(7294)
cl <- makePSOCKcluster(6) 
registerDoParallel(cl) 
rf_tuning <- wf_rf |>
  tune_grid(resamples = comparison_split,
    metrics = covid_metrics,
    grid = rf_grid)

stopCluster(cl)

rf_tuning |> 
  collect_metrics() |> 
  slice_max(mean)

#300 trees is the winner!
write_csv(rf_tuning |> collect_metrics(), file = "rf_tuning.csv")

#xgboost
xgb_grid <- crossing(trees = seq(100, 500, by = 100),
                     tree_depth = seq(10, 30, by = 10),
                     learn_rate = c(0.001, 0.01, 0.1))

set.seed(662)
cl <- makePSOCKcluster(6) 
registerDoParallel(cl)
xgb_tuning <- wf_xgb |>
  tune_grid(resamples = comparison_split,
            metrics = covid_metrics,
            grid = xgb_grid)

stopCluster(cl)

xgb_tuning |> 
  collect_metrics() |> 
  slice_max(mean)

write_csv(xgb_tuning |> collect_metrics(), file = "xgb_tuning.csv")


#lightgbm
lightgbm_grid <- crossing(trees = seq(500, 2500, by = 500),
                     tree_depth = seq(10, 30, by = 10),
                     learn_rate = c(0.001, 0.01, 0.1))

set.seed(982)
cl <- makePSOCKcluster(6) 
registerDoParallel(cl)
lightgbm_tuning <- wf_lightgbm |>
  tune_grid(resamples = comparison_split,
            metrics = covid_metrics,
            grid = lightgbm_grid)

stopCluster(cl)

lightgbm_tuning |> 
  collect_metrics() |> 
  slice_max(mean, by = .metric, n = 3)

write_csv(lightgbm_tuning |> collect_metrics(), file = "lightgbm_tuning.csv")


