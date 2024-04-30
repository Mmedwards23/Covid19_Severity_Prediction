library(bonsai) #lightGBM
library(lightgbm)
library(xgboost)
library(ranger)
library(tidyverse)
library(tidymodels)
library(yardstick)
tidymodels_prefer()

default_wd <- getwd()

train_df <- read_csv("Merged_Dataframes/integer_encoded_training_data.csv")
test_df <- read_csv("Merged_Dataframes/integer_encoded_testing_data.csv")

test_df_list <- list()
for (i in 1:3) {
  name <- str_c("test_set_", i)
  test_df_list[[i]] <- name
  
  test_df_filt <- test_df |> 
    filter(test_set_id == i) |> 
    select(-test_set_id)
  
  assign(name, test_df_filt)
  
}

recipe_default <- recipe(disease_status ~ ., 
                         train_df) |>
  step_zv(all_predictors()) |>  # removes covariates with no variance
  step_unknown(all_nominal_predictors()) |> #set any NAs = unknown
  step_other(all_nominal_predictors(), threshold = 0.05) |>  #anything with >5% occurrence is put under 'other'
  step_normalize(patient_age) |> #normalize age
  step_dummy(all_nominal_predictors()) |> #dummy variables
  step_novel(all_nominal_predictors()) #allow for new factors in test data 

#default logistic
log <- logistic_reg() |>
  set_engine("glm") |> 
  set_mode("classification")

wf_log <- workflow() |>
  add_model(log) |>
  add_recipe(recipe_default)

#elastic net
log_elastic <- logistic_reg(penalty = 0.001, mixture = 0.6) |> 
  set_engine("glmnet") |> 
  set_mode("classification")

wf_elastic <- workflow() |> 
  add_model(log_elastic) |> 
  add_recipe(recipe_default)

#random forest
rf <- rand_forest(trees = 300) |> 
  set_mode("classification") |>
  set_engine("ranger")

wf_rf <- workflow() |> 
  add_model(rf) |> 
  add_recipe(recipe_default)

#XGBoost
xgb <- boost_tree(trees = 200, tree_depth = 10, learn_rate = 0.1) |> 
  set_mode("classification") |> 
  set_engine("xgboost")

wf_xgb <- workflow() |> 
  add_model(xgb) |> 
  add_recipe(recipe_default)

#LightGBM
lightgbm <- boost_tree(trees = 500, tree_depth = 30, learn_rate = 0.01) |> 
  set_mode("classification") |> 
  set_engine("lightgbm")

wf_lightgbm <- workflow() |> 
  add_model(lightgbm) |> 
  add_recipe(recipe_default)

#fit the models 
model_names <- c("log", "elastic", "rf",
                  "xgb", "lightgbm")
workflows <- list(wf_log, wf_elastic, wf_rf,
                       wf_xgb, wf_lightgbm)
workflows_tbl <- tibble(model_names = model_names,
                        workflow_objects = workflows)

library(doParallel)
cl <- makePSOCKcluster(6) 
registerDoParallel(cl) 

set.seed(5581)

workflows_tbl_fit <- workflows_tbl |>
  rowwise() |>
  mutate(fits = list(fit(workflow_objects, 
                         train_df)))

#test set 1
set.seed(8247)
workflows_tbl_fit <- workflows_tbl_fit |>
  mutate(pred_01_class = list(predict(fits, 
                                      test_set_1, 
                                      type = "class"))) |>
  mutate(pred_01_prob = list(predict(fits,
                                     test_set_1, 
                                     type = "prob")))

workflows_tbl_fit <- workflows_tbl_fit |>
  mutate(predictions = list(bind_cols(pred_01_class, pred_01_prob))) |>
  select(-c(pred_01_class, pred_01_prob))

predictions_tbl_1 <- workflows_tbl_fit |>
  select(model_names, 
         predictions) |>
  unnest(cols = c(predictions)) |>
  cbind(Disease_status = test_set_1 |>
          pull(disease_status)) |> 
  mutate(Disease_status = as.factor(Disease_status))

#test set 2
set.seed(7742)
workflows_tbl_fit <- workflows_tbl_fit |>
  mutate(pred_02_class = list(predict(fits, 
                                      test_set_2, 
                                      type = "class"))) |>
  mutate(pred_02_prob = list(predict(fits,
                                     test_set_2, 
                                     type = "prob")))

workflows_tbl_fit <- workflows_tbl_fit |>
  mutate(predictions = list(bind_cols(pred_02_class, pred_02_prob))) |>
  select(-c(pred_02_class, pred_02_prob))

predictions_tbl_2 <- workflows_tbl_fit |>
  select(model_names, 
         predictions) |>
  unnest(cols = c(predictions)) |>
  cbind(Disease_status = test_set_2 |>
          pull(disease_status)) |> 
  mutate(Disease_status = as.factor(Disease_status))

set.seed(209)
#test set 3
workflows_tbl_fit <- workflows_tbl_fit |>
  mutate(pred_03_class = list(predict(fits, 
                                      test_set_3, 
                                      type = "class"))) |>
  mutate(pred_03_prob = list(predict(fits,
                                     test_set_3, 
                                     type = "prob")))

workflows_tbl_fit <- workflows_tbl_fit |>
  mutate(predictions = list(bind_cols(pred_03_class, pred_03_prob))) |>
  select(-c(pred_03_class, pred_03_prob))

predictions_tbl_3 <- workflows_tbl_fit |>
  select(model_names, 
         predictions) |>
  unnest(cols = c(predictions)) |>
  cbind(Disease_status = test_set_3 |>
          pull(disease_status)) |> 
  mutate(Disease_status = as.factor(Disease_status))

stopCluster(cl)

predictions_tbl_1 <- predictions_tbl_1 |> 
  mutate(id = 1)

predictions_tbl_2 <- predictions_tbl_2 |> 
  mutate(id = 2)

predictions_tbl_3 <- predictions_tbl_3 |> 
  mutate(id = 3)

predictions_tbl <- rbind(predictions_tbl_1, predictions_tbl_2, predictions_tbl_3)

write_csv(predictions_tbl, "Merged_Dataframes/Predictions.csv")

#Bootstrapping coefficient estimates
#reuse recipe_default, wf_elastic

set.seed(5581)
#note, same seed as initial fit, although since being used with bootstrap
#it wont be the exact same... 
control_fit <- control_resamples(extract = tidy)

covid_bootstrap <- bootstraps(train_df, 
                             times = 20) #its a big dataset, so 20?

bootstrap_elastic <- wf_elastic |>
  fit_resamples(covid_bootstrap,
                control = control_fit)

bootstrap_e_coefs <- bootstrap_elastic |>
  select(id, .extracts) |>
  unnest(.extracts) |> 
  unnest(.extracts)

write_csv(bootstrap_e_coefs, "Merged_Dataframes/elastic_coefs.csv")


