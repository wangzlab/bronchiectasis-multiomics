library(data.table)
library(dplyr)

# profiles ==================================
profile_metat_Pm <- fread("Pm_module.txt", data.table = F) %>% tibble::column_to_rownames("#SampleID") %>% t() %>% as.data.frame()
profile_metat_Rm <- fread("Rm_module.txt", data.table = F) %>% tibble::column_to_rownames("#SampleID") %>% t() %>% as.data.frame()
profile_metat_Vd <- fread("Vd_module.txt", data.table = F) %>% tibble::column_to_rownames("#SampleID") %>% t() %>% as.data.frame()
profile_metat_Vp <- fread("Vp_module.txt", data.table = F) %>% tibble::column_to_rownames("#SampleID") %>% t() %>% as.data.frame()

profile_metab <- fread("../_dataGZ/metaB_module.txt", data.table = F)
profile_metab <- profile_metab %>%
  tibble::column_to_rownames("V1") %>%
  t() %>% as.data.frame()

profile_hostt <- fread("../_dataGZ/hostT_pathway.txt", data.table = F)
profile_hostt <- profile_hostt %>%
  tibble::column_to_rownames("#") %>%
  t() %>% as.data.frame()

# predictors ================================
Pm_metat_metab <- fread("prediction_Pm_metat_metab.txt", data.table = F, header = F, sep = "-", col.names = c("metat","metab")); head(Pm_metat_metab)
Pm_metat_metab_hostt <- fread("prediction_Pm_metat_metab_hostt.txt", data.table = F, header = F, col.names = c("metat","metab","hostt")); head(Pm_metat_metab_hostt)
Rm_metat_metab <- fread("prediction_Rm_metat_metab.txt", data.table = F, header = F, sep = "-", col.names = c("metat","metab")); head(Rm_metat_metab)
Rm_metat_metab_hostt <- fread("prediction_Rm_metat_metab_hostt.txt", data.table = F, header = F, col.names = c("metat","metab","hostt")); head(Rm_metat_metab_hostt)
Vd_metat_metab <- fread("prediction_Vd_metat_metab.txt", data.table = F, header = F, sep = "-", col.names = c("metat","metab")); head(Vd_metat_metab)
Vd_metat_metab_hostt <- fread("prediction_Vd_metat_metab_hostt.txt", data.table = F, header = F, col.names = c("metat","metab","hostt")); head(Vd_metat_metab_hostt)
Vp_metat_metab <- fread("prediction_Vp_metat_metab.txt", data.table = F, header = F, sep = "-", col.names = c("metat","metab")); head(Vp_metat_metab)
Vp_metat_metab_hostt <- fread("prediction_Vp_metat_metab_hostt.txt", data.table = F, header = F, col.names = c("metat","metab","hostt")); head(Vp_metat_metab_hostt)

# Y: Pa Hi =====================================
Pa.Hi <- fread("../_dataGZ/Pa_Hi.txt", data.table = F)
Pa.Hi <- Pa.Hi %>% tibble::column_to_rownames("#") %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("SampleID")
head(Pa.Hi)

# Covars =====================================
meta <- fread("../_dataGZ/full_metadata_table.txt", data.table = F)

# groups to extract Stable ===================
groups <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
meta$SampleID %in% groups$Sample

Pa.Hi$SampleID  %in% groups$Sample
Pa.Hi <- Pa.Hi %>%
  mutate(Cluster = sapply(SampleID, function(x) groups$Group[which(groups$Sample == x)])) %>%
  select(SampleID, Cluster, Pa_kraken, Hi_kraken)

Pa.Hi <- merge(Pa.Hi, meta %>% select(SampleID, Age, Gender, BMI), by = "SampleID")

meta.Stable <- Pa.Hi %>% filter(Cluster == "Stable")


# check IDs ===============================
rownames(profile_hostt) %in% groups$Sample
rownames(profile_metab) %in% groups$Sample
rownames(profile_metat_Pm) %in% groups$Sample
rownames(profile_metat_Rm) %in% groups$Sample
rownames(profile_metat_Vd) %in% groups$Sample
rownames(profile_metat_Vp) %in% groups$Sample

metatMatch <- fread("../_dataGZ/metat_match.txt", data.table = F,header = F)
rownames(profile_metat_Pm) <- sapply(rownames(profile_metat_Pm), function(x) metatMatch$V2[which(metatMatch$V1 == x)])
rownames(profile_metat_Rm) <- sapply(rownames(profile_metat_Rm), function(x) metatMatch$V2[which(metatMatch$V1 == x)])
rownames(profile_metat_Vd) <- sapply(rownames(profile_metat_Vd), function(x) metatMatch$V2[which(metatMatch$V1 == x)])
rownames(profile_metat_Vp) <- sapply(rownames(profile_metat_Vp), function(x) metatMatch$V2[which(metatMatch$V1 == x)])


# prediction ===============================

library("tidymodels")
library(tidyverse)
library("workflows")
library("tune")
library(ranger)


PredictionType = "regression"


for(spc in c("Pm","Rm","Vd","Vp")){  # "Pm","Rm","Vd","Vp"
  # spc = "Pm"
  
  if(grepl("Rm",spc)){
    Y = "Hi_kraken"
    col.rm = "Pa_kraken"
  }else if(grepl("(Vd|Vp|Pm)",spc)){
    Y = "Pa_kraken"
    col.rm = "Hi_kraken"
  }else{
    print("Wrong file name")
    break
  }
  
  profile_metat <- eval(parse(text = paste0("profile_metat_",spc)))
  omics <- merge(merge(profile_metat, profile_metab, by=0),
                 profile_hostt, by.x = "Row.names", by.y = 0)  %>%
    tibble::column_to_rownames("Row.names")
  
  metat_metab <- eval(parse(text = paste0(spc, "_metat_metab")))
  metat_metab_hostt <- eval(parse(text = paste0(spc, "_metat_metab_hostt")))
  pairs <- bind_rows(metat_metab_hostt, metat_metab)
  pairs$metab <- sub("MetaB_","", pairs$metab)
  
  RF.performance <- NULL
  
  for(i_prs in 1:nrow(pairs)){
    # i_prs = 1
    metat <- pairs$metat[i_prs]
    metab <- pairs$metab[i_prs]
    hostt <- pairs$hostt[i_prs]
    
    if(!metat %in% colnames(omics) | !metab %in% colnames(omics)) next
    if(is.na(hostt)){
      rfDat <- merge(meta.Stable %>% dplyr::select(-all_of(col.rm), -Cluster), 
                     omics %>% dplyr::select(all_of(c(metat, metab))), 
                     by.x = "SampleID", by.y=0)
        
    }else{
      rfDat <- merge(meta.Stable %>% dplyr::select(-all_of(col.rm), -Cluster), 
                     omics %>% dplyr::select(all_of(c(metat, metab, hostt))),
                     by.x = "SampleID", by.y=0)
    }
    
    colnames(rfDat)[colnames(rfDat) == Y] = "Y"
    head(rfDat)
    
    # create a cross validation version of the training set for parameter tuning
    my_cv <- vfold_cv(rfDat, v = 5, repeats = 2)
    rfDat$SampleID <- NULL
    
    # define the recipe
    my_recipe <- 
      recipe(Y ~ ., data = rfDat) %>%
      step_impute_median(all_numeric_predictors()) %>% # 处理数值变量的缺失值：用中位数填充
      step_impute_mode(all_nominal_predictors()) # 处理分类变量的缺失值：用众数填充
    
    # Specify the model --------------------------------------
    rf_model <- 
      # specify that the model is a random forest
      rand_forest() %>%  # ?rand_forest shows the tunable parameters
      # specify that the `mtry` parameter needs to be tuned
      set_args(mtry = tune(), trees=tune(),min_n=5) %>%
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
    # if(is.na(hostt)) mtrys = c(1,2) else mtrys = c(1,2,3) 
    rf_grid <- expand.grid(mtry =  c(1,2,3) , trees=c(100, 200, 500,1000)) 
    # extract results
    if(PredictionType == "regression"){
      rf_tune_results <- rf_workflow %>%
        tune_grid(resamples = my_cv, #CV object
                  grid = rf_grid, # grid of values to try
                  metrics = metric_set(rmse, mae, rsq)# metrics we care about
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
      mutate(Species = spc, 
             Y = Y,
             MetaT = metat,
             MetaB = metab,
             HostT = hostt)
    
    RF.performance <- bind_rows(RF.performance, res_c)
  }# loop through pairs
  
  write.csv(RF.performance, file = paste0("RF.performance_PaHi.inStable_",spc,".csv"), quote = F, row.names = F)
  
}#loop through species
