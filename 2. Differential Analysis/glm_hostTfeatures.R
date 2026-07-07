library(data.table)
library(dplyr)

tmp1 <- fread("../_dataGZ/meta_confounder.txt", data.table = F)
meta_origin <- fread("../_dataGZ/group_all.txt", data.table = F)
meta_origin <- merge(meta_origin %>% select(-MetaT), 
                     tmp1 %>% select(-Group),
                     by.x="#NAME", by.y="SampleID")

grouping <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
colnames(grouping)[2] <- "Grouping_new"

sample_to_group <- setNames(grouping$Grouping_new, grouping$Sample)
meta_origin$MicrobGrouping <- sample_to_group[meta_origin$`#NAME`]

print(table(meta_origin$MicrobGrouping, useNA = "always"))

dat <- fread("../_dataGZ/hostT_feature.txt", data.table = F)
colnames(dat)[1] <- "#NAME"
colnames(dat) %in% meta_origin$`#NAME`


feature_abb_df <- 
  cbind.data.frame(
    feature = dat$`#NAME`,
    abb = paste0("feature", seq(1,nrow(dat), 1)),
    stringsAsFactors=F
  )
write.table(feature_abb_df, file = "hostT_feature.abb.txt", sep = "\t",quote = F, row.names = F)


feature_df.full <-  
  dat %>%
  mutate(feature = sapply(`#NAME`, function(x) feature_abb_df$abb[which(feature_abb_df$feature == x)])) %>%
  select(-`#NAME`) %>%
  tibble::column_to_rownames("feature") %>%
  t() %>%
  as.data.frame(stringsAsFactors = F)


sps.PPM <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "PPM")]
sps.Pa <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "Pa")]
sps.Hi <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "Hi")]
sps.Commensal <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "Commensal")]

sps.Health <- meta_origin$`#NAME`[which(meta_origin$Type == "Health")]
sps.Stable <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable")]
sps.Exacerbation <- meta_origin$`#NAME`[which(meta_origin$Type == "Exacerbation")]

comparison_list <- list(
  c("Pa", "PPM"),     
  c("Hi", "PPM"),      
  c("PPM", "Health"),  
  c("PPM", "Commensal"),
  c("Pa", "Commensal"),  
  c("Hi", "Commensal"),
  c("Commensal","Health"),
  c("Stable","Health"),
  c("Pa","Health"),
  c("Hi","Health")
)



for(cp in comparison_list){  
  # cp <- comparison_list[[1]]
  Groups_df <- bind_rows(
    cbind.data.frame(
      Sample=eval(parse(text = paste0("sps.",cp[1]))),
      Group=cp[1],
      stringsAsFactors=F
    ),
    cbind.data.frame(
      Sample=eval(parse(text = paste0("sps.",cp[2]))),
      Group=cp[2],
      stringsAsFactors=F
    )
  )
  sps.avail <- intersect(Groups_df$Sample, rownames(feature_df.full))
  
  # match sample sequences:
  feature_df <- feature_df.full[match(sps.avail, rownames(feature_df.full)),]
  meta <- meta_origin[match(sps.avail, meta_origin$`#NAME`),]
  
  # loop through features to calculate glm
  pval_res <- NULL
  for(ft in colnames(feature_df)){
    # ft = colnames(feature_df)[1]
    
    testDat <- 
      cbind.data.frame(
        feature =  feature_df[,ft],
        meta ,
        stringsAsFactors=F
      ) %>%
      mutate(Group = sapply(`#NAME`, function(x) Groups_df$Group[which(Groups_df$Sample == x)]))
    testDat$Group <- factor(testDat$Group, levels = c(cp[2], cp[1]))
    
    g <- glm(feature ~ Group + Current_smoking + Gender + BMI + Age, data = testDat)
    g.s <- summary(g)
    p.glm <- g.s$coefficients[2,4]
    # p.glm.name <- rownames( g.s$coefficients)[2]
    # if(p.glm.name != "GroupStable") print(paste0("p glm name is ", p.glm.name))
    
    w <- wilcox.test(feature ~ Group, data = testDat) 
    p.wilcox <- w$p.value
    
    avg <- 
      testDat %>%
      group_by(Group) %>%
      summarise(avg = mean(feature))  %>%
      reshape2::dcast(.~Group, value.var = "avg") %>%
      select(-1)
    
    colnames(avg) <- paste0("avg_", colnames(avg))
    
    res_c <- 
      cbind.data.frame(
        Feature = ft, 
        glm.p = p.glm,
        wilcox.p = p.wilcox,
        avg )
    
    pval_res <- bind_rows(pval_res, res_c)
  }
  
  pval_res$DataSource = "hostT_feature"
  pval_res$Feature.full = sapply(pval_res$Feature, 
                                 function(x) feature_abb_df$feature[which(feature_abb_df$abb == x)])
  pval_res$Feature <- NULL
  pval_res$glm.fdr = p.adjust(pval_res$glm.p, method = "fdr")
  pval_res$wilcox.fdr = p.adjust(pval_res$wilcox.p, method = "fdr")
  write.table(pval_res,
              file = paste0("hostT_feature_", paste(cp[1],cp[2], sep = ".vs."), ".txt"),
              sep="\t", quote = F, row.names = F)
}


# stable vs Exacerbations -----------------
# paired wilcoxon  and lmer.p

cp =  c("Stable","Exacerbation")

Groups_df <- bind_rows(
  cbind.data.frame(
    Sample=eval(parse(text = paste0("sps.",cp[1]))),
    Group=cp[1],
    stringsAsFactors=F
  ),
  cbind.data.frame(
    Sample=eval(parse(text = paste0("sps.",cp[2]))),
    Group=cp[2],
    stringsAsFactors=F
  )
)

sps.avail <- intersect(Groups_df$Sample, rownames(feature_df.full))

# match sample sequences:
feature_df <- feature_df.full[match(sps.avail, rownames(feature_df.full)),]
meta <- meta_origin[match(sps.avail, meta_origin$`#NAME`),]

library(lmerTest)

Stable.Exacerb_pval <- NULL
for(ft in colnames(feature_df)[1:ncol(feature_df)]){
  # ft = colnames(feature_df)[1]
  
  print(which(colnames(feature_df) == ft))
  testDat <- 
    cbind.data.frame(
      feature =  feature_df[,ft],
      meta,
      stringsAsFactors=F
    ) %>%
    mutate(Group = sapply(`#NAME`, function(x) Groups_df$Group[which(Groups_df$Sample == x)]))
  
  testDat.w <- testDat %>% reshape2::dcast(SubjectID~Group, value.var = "feature")
  testDat.w <- testDat.w[complete.cases(testDat.w),]
  
  avg <- 
    testDat.w %>%
    reshape2::melt(id.var = "SubjectID", variable.name="Group") %>%
    group_by(Group) %>%
    summarise(avg = mean(value))  %>%
    reshape2::dcast(.~Group, value.var = "avg") %>%
    select(-1)
  
  colnames(avg) <- paste0("avg_", colnames(avg))
  
  # the paired wilcox p
  w <- wilcox.test(x = testDat.w$Exacerbation, y = testDat.w$Stable, paired = TRUE) 
  p.wilcox <- w$p.value
  
  
  # the lmr p-value
  testDat <- testDat %>% 
    filter(SubjectID %in% testDat.w$SubjectID) 
  
  ml1 <- tryCatch({
    lmer(feature ~ Group + (1|SubjectID) , data=testDat)
  }, error = function(e) {
    NULL
  })
  
  if(is.null(ml1)) next
  
  ml1.s <- summary(ml1, ddf="Kenward-Roger") 
  p.lmer <- ml1.s$coefficients[2,"Pr(>|t|)"] 
  
  
  res_c <- 
    cbind.data.frame(
      Feature = ft, 
      # glm.p = p.glm,
      wilcox.p = p.wilcox,
      lmer.p = p.lmer,
      avg )
  
  Stable.Exacerb_pval <- bind_rows(Stable.Exacerb_pval, res_c)
}

Stable.Exacerb_pval$DataSource = "hostT_feature"
Stable.Exacerb_pval$Feature.full = sapply(Stable.Exacerb_pval$Feature,  function(x) feature_abb_df$feature[which(feature_abb_df$abb == x)])
Stable.Exacerb_pval$Feature <- NULL

Stable.Exacerb_pval$lmer.fdr = p.adjust(Stable.Exacerb_pval$lmer.p, method = "fdr")
Stable.Exacerb_pval$wilcox.fdr = p.adjust(Stable.Exacerb_pval$wilcox.p, method = "fdr")

write.table(Stable.Exacerb_pval,
            file = "hostT_feature_Stable.vs.Exacerbations.txt",
            sep="\t", quote = F, row.names = F)
