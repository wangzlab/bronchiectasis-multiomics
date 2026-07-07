list.files()

library(data.table)
library(dplyr)
library(pROC)
groups <- fread("NPD_Health_EMBARC_CAMEB2_China.txt", data.table = F)
groups$Cohort %>% unique

CIfiles <- list.files("BDI_diff_cutoffs/", full.names = T)

dataSets <- list(c("EMBARC", "China", "CAMEB2"), c("EMBARC"), c("China"), c("CAMEB2"))


AUC.res <- NULL
for(ds in dataSets){
  # ds = c("EMBARC", "China")
  
  groups_sub <- groups %>% filter(Cohort %in% ds)
  
  for(f in CIfiles){
    # f = CIfiles[1]
    
    dat.CI <- fread(f, data.table = F) %>% rename(CI = V5)
    
    testD <- merge(groups_sub, dat.CI, by.x = "SampleID", by.y = "V1")
    
    testD$is_positive <- ifelse(testD$Group == "Commensal", 1, 0)
    
    # 计算ROC曲线和AUC
    roc_obj <- roc(response = testD$is_positive, predictor = testD$CI)
    auc(roc_obj)
    
    res_c <- data.frame(cohort = paste(ds, collapse = " & "),
                        file = basename(f),
                        AUC = as.numeric(auc(roc_obj)))
    AUC.res <- bind_rows(AUC.res, res_c)
  }
}

write.csv(AUC.res, file = 'AUC_different_BDI.csv', quote = F, row.names = F)
