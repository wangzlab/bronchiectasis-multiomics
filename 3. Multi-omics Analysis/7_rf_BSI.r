library(data.table)
library(dplyr)


library("tidymodels")
library(tidyverse)
library("workflows")
library("tune")
library(ranger)


metab <- fread("prediction_metab.txt", data.table = F)
metat <- fread("prediction_metat.txt", data.table = F)

groups <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
groups$new_grouping_info_cutoff10 %>% unique

colnames(metab)[!colnames(metab) %in% groups$Sample]
colnames(metat)[!colnames(metat) %in% groups$Sample]

metab <- metab %>% tibble::column_to_rownames("#") %>% t() %>% as.data.frame()
metat <- metat %>% tibble::column_to_rownames("#") %>% t() %>% as.data.frame()


meta.BSI <- fread("../_dataGZ/full_metadata_table.txt", data.table = F)
meta.BSI <- meta.BSI %>% dplyr::select(SampleID, Age, Gender, BMI,BSI)
meta.BSI$BSI <- as.integer(meta.BSI$BSI)
meta.BSI$SampleID[!meta.BSI$SampleID %in% groups$Sample]

files <- list.files("prediction/", full.names = T)

groups_vec <- setNames(c("Hi","Commensal","Pa","PPM"), nm = c("HID","NPD","PAD","PPMD"))
PredictionType = "regression"


RF.performance <- NULL

for(g in c("HID","NPD","PAD","PPMD")){
  # g = "HID"
  
  pf = files[grep(g, files)]
  pairs <- fread(pf, data.table = F, header = F, col.names = c("metat","metab"))
  
  sps.grp <- groups$Sample[groups$new_grouping_info_cutoff10 == groups_vec[g]]
  sps.grp.mb <- intersect(sps.grp, rownames(metab))
  sps.grp.mt <- intersect(sps.grp, rownames(metat))
  
  metab.sub <- metab[sps.grp.mb,]
  metat.sub <- metat[sps.grp.mt,]
  
  for(i in 1:nrow(pairs)){
    # i=1
    mb = pairs$metab[i]
    mt = pairs$metat[i]
    
    omicDat <- merge(metab.sub %>% dplyr::select(all_of(mb)),
          metat.sub %>% dplyr::select(all_of(mt)),
          by = 0)
    rfDat <- merge(omicDat, meta.BSI, by.x = "Row.names", by.y = "SampleID") %>%
      tibble::column_to_rownames("Row.names") %>%
      rename(Y = BSI)
    
    # create a cross validation version of the training set for parameter tuning
    my_cv <- vfold_cv(rfDat, v = 5)
    
    # define the recipe
    my_recipe <- 
      # which consists of the formula (outcome ~ predictors)
      recipe(Y ~ ., data = rfDat)  # %>%
    #step_knnimpute(all_predictors()) # missing value
    
    # Specify the model --------------------------------------
    rf_model <- 
      # specify that the model is a random forest
      rand_forest() %>%  # ?rand_forest shows the tunable parameters
      # specify that the `mtry` parameter needs to be tuned
      set_args(mtry = tune(), trees=tune(),min_n=10) %>%
      # select the engine/package that underlies the model
      set_engine("ranger", importance = "impurity") %>% # if to examine the variable importance of your final model need to set importance = 
      # choose either the continuous regression or binary classification mode
      set_mode(PredictionType) 
    
    #  Put it all together in a workflow -----------------------------------
    # set the workflow
    rf_workflow <- workflow() %>%
      # add the recipe
      add_recipe(my_recipe) %>%
      # add the model
      add_model(rf_model) 
    
    # Tune the parameters ----------------------------------
    # specify which values eant to try
    rf_grid <- expand.grid(mtry = c(1,2,3), trees=c(100, 500,1000,1500)) 
    # extract results
    if(PredictionType == "regression"){
      rf_tune_results <- rf_workflow %>%
        tune_grid(resamples = my_cv, #CV object
                  grid = rf_grid, # grid of values to try
                  metrics = metric_set( rmse, mae, rsq) # metrics we care about
        )
      resultSelect = "rmse"
    }else{
      rf_tune_results <- rf_workflow %>%
        tune_grid(resamples = my_cv, #CV object
                  grid = rf_grid, # grid of values to try
                  metrics = metric_set(accuracy, roc_auc) # metrics we care about
        )
      resultSelect = "roc_auc"
    }
    
    # print results
    results_df <- rf_tune_results %>%
      collect_metrics()
    
    # Finalize the workflow ------------------------------------------
    param_final <- rf_tune_results %>%
      select_best(metric = resultSelect) # 采用最高accuracy的参数
    param_final
    
    res_c <- results_df %>%
      filter(mtry == param_final$mtry,
             trees == param_final$trees,
             .config == param_final$.config) %>%
      reshape2::dcast( mtry + trees + `.config` ~ `.metric`, value.var = "mean") %>%
      mutate(Group = g, 
             MetaB = mb,
             MetaT = mt)
    
    RF.performance <- bind_rows(RF.performance, res_c)
  } # loop through pairs
}# loop through groups


write.csv(RF.performance, file = "rf_BSI.csv", quote = F, row.names = F)
