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

dat <- fread("../_dataGZ/metaB_module.txt", data.table = F)
colnames(dat)[1] <- "#NAME" 


feature_df.full <-  
  dat %>%
  tibble::column_to_rownames("#NAME") %>%
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
  c("Stable","Health")
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
  
  # 对每一个Feature,算glm
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
  
  pval_res$DataSource = "metaB_module"
  write.table(pval_res,
              file = paste0("metaB_module_", paste(cp[1],cp[2], sep = ".vs."), ".txt"),
              sep="\t", quote = F, row.names = F)
}
