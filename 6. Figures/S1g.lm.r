library(data.table)
library(dplyr)


meta.cohorts <- fread("../_data/combine_metadata(1).txt", data.table = F) # groups and cohorts

# clinical meta --------------
vars <- c("BSI","AE_1year","FEV1pred_perc","MMRC")
covars <- c("Age","Gender")

meta.617 <- fread("../_data/combinedData_617sps.csv", data.table = F) 
meta.617 <- meta.617 %>% 
  rename(AE_1year = Exacerbation_frequency, FEV1pred_perc = FEV1_percent_predicted, SampleID=Sequencing.ID) %>% 
  dplyr::select(SampleID, Study, all_of(vars), all_of(covars))
sapply(meta.617, class)

meta.EMBARC <- meta.617 %>% filter(Study == "EMBARC")
meta.CAMEB2 <- meta.617 %>% filter(Study == "CAMEB2")

# meta.118 <- fread("../metadata_118sps.txt", data.table = F)
meta.118 <- readxl::read_excel("../_data/metadata_118sps_more.xlsx")  # 比 "../metadata_118sps.txt"要全，数据一致
meta.118$Group %>% table()
meta.118.stable <- meta.118 %>% filter(Group == "Stable")  %>% dplyr::select(SampleID, all_of(vars), all_of(covars))
sapply(meta.118.stable, class)
meta.118.stable$BSI <- as.integer(meta.118.stable$BSI)

# species --------------

species <- fread("../_data/combine_metagenome_rel(1).txt", data.table = F)
species$species <- sub(".*;s__","",species$`#NAME`)
test <- species %>% dplyr::select(species, `#NAME`)

species$`#NAME` <- NULL
species.t <- species %>% tibble::column_to_rownames("species") %>% t() %>% as.data.frame()

all(meta.118.stable$SampleID %in% rownames(species.t))
all(meta.CAMEB2$SampleID %in% rownames(species.t))

all(meta.EMBARC$SampleID %in% rownames(species.t))
meta.EMBARC$SampleID[!meta.EMBARC$SampleID %in% rownames(species.t)] # checked 

rownames(species.t)[!rownames(species.t) %in% c(meta.118$SampleID,meta.CAMEB2$SampleID,meta.EMBARC$SampleID)]


library(vegan)
hellinger_data <- decostand(species.t, "hellinger") # 直接对含零的数据操作

library(MASS)

res_BSI <- NULL
res_BSI.cat <- NULL
res_AE_1year <- NULL
res_FEV1pred <- NULL
res_mMRC <- NULL

for(dn in c("meta.118.stable", "meta.CAMEB2","meta.EMBARC")){  # "meta.118.stable","meta.CAMEB2","meta.EMBARC"
  # dn = "meta.118.stable"
  
  meta <- eval(parse(text = dn))
  meta$BSI.cat <- cut(meta$BSI,
                      breaks = c(-Inf,5,9,Inf),
                      labels = c("mild","moderate","severe"),
                      include.lowest = T, right = F)
  meta$BSI.cat <- factor(meta$BSI.cat, levels = c("mild","moderate","severe"))
  meta$MMRC <- factor(paste0("level", meta$MMRC), levels = c("level0","level1","level2","level3","level4","level5"))
  fullTestDat <- merge(meta, hellinger_data, by.x = "SampleID", by.y = 0)
  
  for(spc in colnames(hellinger_data)[1:ncol(hellinger_data)]){
    # spc = colnames(hellinger_data)[1]
    which(colnames(hellinger_data) == spc)
    
    testD <- fullTestDat %>% rename(Species = !!spc) %>% dplyr::select(SampleID, all_of(c(vars, "BSI.cat", covars, "Species")))
    if(sum(testD$Species > 0) <=1 ) next
    
    # species vs. BSI ----------------------------------
    robust_model_bsi <- rlm(BSI ~ Species + Age + Gender, data = testD)
    robust_summary <- summary(robust_model_bsi) # 获取详细结果
    coefficients <- robust_summary$coefficients
    
    n <- nrow(testD)  
    p <- length(coefficients[, "Value"])  
    df_residual <- n - p  
    
    p_values <- 2 * pt(abs(coefficients[, "t value"]), 
                       df = df_residual, 
                       lower.tail = FALSE)
    
    
    res_c_BSI <- data.frame(
      dataSet = sub("meta.", "", dn),
      var = "BSI",
      Species = spc, 
      Method = "MASS::rlm",
      Estimate = coefficients["Species", "Value"],
      Std.Error = coefficients["Species", "Std. Error"],
      t.value = coefficients["Species", "t value"],
      p.value = p_values["Species"]
    )
    res_BSI <- bind_rows(res_BSI, res_c_BSI)
    
    # species vs. BSI.cat ----------------------------------
    ordinal_model <- 
      try(polr(BSI.cat ~ Species + Age + Gender, data = testD, Hess = TRUE, 
               control = list(maxit = 10000)), # 注意：Hess = TRUE 是为了后续计算p值存储海森矩阵
          silent = TRUE
      )
    
    if(inherits(ordinal_model, "try-error")) {
      res_c_BSI.cat <- data.frame(
        dataSet = sub("meta.", "", dn),
        var = "BSI.cat",
        Species = spc, 
        Method = "MASS::polr",
        Estimate = NA,
        Std.Error = NA,
        t_value = NA,
        p_value = NA,
        OR = NA,
        CI_low = NA,
        CI_high = NA
      )
    }else{
      # 计算p值（需要手动完成）
      ctable <- coef(summary(ordinal_model)) # 1. 获取系数摘要
      p_values <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2 # 2. 计算p值 
      ctable_complete <- cbind(ctable, "p value" = p_values)   # 3. 合并结果
      
      # 使用基于正态近似的置信区间（而不是轮廓似然区间）
      ci_default <- confint.default(ordinal_model)  # 使用confint.default而不是confint
      odds_ratios <- exp(coef(ordinal_model))
      ci_odds <- exp(ci_default)
      
      res_c_BSI.cat <- data.frame(
        dataSet = sub("meta.", "", dn),
        var = "BSI.cat",
        Species = spc, 
        Method = "MASS::polr",
        Estimate = ctable_complete["Species", "Value"],
        Std.Error = ctable_complete["Species", "Std. Error"],
        t_value = ctable_complete["Species", "t value"],
        p_value = p_values["Species"],
        OR = odds_ratios["Species"],
        CI_low = ci_odds["Species",1],
        CI_high = ci_odds["Species",2]
      )
      print(res_c_BSI.cat)
    }
    
    res_BSI.cat <- bind_rows(res_BSI.cat, res_c_BSI.cat)
    
    # species vs. AE_1year ----------------------------------
    nb_model <- glm.nb(AE_1year ~ Species + Age + Gender, data = testD)  
    model_summary <- summary(nb_model)
    coefficients <- model_summary$coefficients
    species_results <- coefficients["Species", ]
    irr <- exp(species_results["Estimate"]) 
    ci_default <- confint.default(nb_model)["Species", ]
    ci_irr <- exp(ci_default)
    
    res_c_AE_1year <- data.frame(
      dataSet = sub("meta.", "", dn),
      var = "AE_1year",
      Species = spc, 
      Method = "MASS::glm.nb",
      Estimate = species_results["Estimate"],
      Std.Error = species_results["Std. Error"],
      z_value = species_results["z value"],
      IRR = irr,
      CI_low = ci_irr[1],    
      CI_high = ci_irr[2],   
      p_value = species_results["Pr(>|z|)"]
    )
    res_AE_1year <- bind_rows(res_AE_1year, res_c_AE_1year)
    
    # species vs. FEV1%pred  ----------------------------------
    robust_model_bsi <- rlm(FEV1pred_perc ~ Species + Age + Gender, data = testD)
    robust_summary <- summary(robust_model_bsi) # 获取详细结果
    coefficients <- robust_summary$coefficients
    
    # 计算p值
    n <- nrow(testD) 
    p <- length(coefficients[, "Value"]) 
    df_residual <- n - p  
    p_values <- 2 * pt(abs(coefficients[, "t value"]), 
                       df = df_residual, 
                       lower.tail = FALSE)
    
    res_c_FEV1pred <- data.frame(
      dataSet = sub("meta.", "", dn),
      var = "FEV1pred_perc",
      Species = spc, 
      Method = "MASS::rlm",
      Estimate = coefficients["Species", "Value"],
      Std.Error = coefficients["Species", "Std. Error"],
      t.value = coefficients["Species", "t value"],
      p.value = p_values["Species"]
    )
    res_FEV1pred <- bind_rows(res_FEV1pred, res_c_FEV1pred)
    
    # species vs. mMRC ------------------------------------------
    ordinal_model <-
      try(polr(MMRC ~ Species + Age + Gender, data = testD, Hess = TRUE, 
               control = list(maxit = 10000)),  # 注意：Hess = TRUE 是为了后续计算p值存储海森矩阵
          silent = TRUE
      )
    
    if(inherits(ordinal_model, "try-error")){
      res_c_mMRC <- data.frame(
        dataSet = sub("meta.", "", dn),
        var = "MMRC",
        Species = spc, 
        Method = "MASS::polr",
        Estimate = NA,
        Std.Error = NA,
        t_value = NA,
        p_value = NA,
        OR = NA,
        CI_low = NA,
        CI_high = NA
      )
    }else{
      ctable <- coef(summary(ordinal_model)) 
      p_values <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2 
      ctable_complete <- cbind(ctable, "p value" = p_values)  
      
      ci_default <- confint.default(ordinal_model) 
      odds_ratios <- exp(coef(ordinal_model))
      ci_odds <- exp(ci_default)
      
      res_c_mMRC <- data.frame(
        dataSet = sub("meta.", "", dn),
        var = "MMRC",
        Species = spc, 
        Method = "MASS::polr",
        Estimate = ctable_complete["Species", "Value"],
        Std.Error = ctable_complete["Species", "Std. Error"],
        t_value = ctable_complete["Species", "t value"],
        p_value = p_values["Species"] ,
        OR = odds_ratios["Species"],
        CI_low = ci_odds["Species",1],
        CI_high = ci_odds["Species",2]
      )
    }
    
    res_mMRC <- bind_rows(res_mMRC, res_c_mMRC)
  }
}


save(res_BSI, res_BSI.cat, res_FEV1pred, res_mMRC, res_AE_1year,
     file = "S1g.lmCovars.Res.RData")