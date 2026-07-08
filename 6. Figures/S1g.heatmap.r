load("S1g.lmCovars.Res.RData") 
load("../Figure1/1a.topSpecies.1e-2.RData")

# keep only topSpecies
res_BSI <- res_BSI %>% filter(Species %in% speciesTop) 
res_AE_1year <- res_AE_1year %>% filter(Species %in% speciesTop) 
res_FEV1pred <- res_FEV1pred %>% filter(Species %in% speciesTop)
res_mMRC <- res_mMRC %>% filter(Species %in% speciesTop)

res_BSI <- res_BSI %>% 
  mutate(Zscore = -qnorm(p.value/2) ) %>%
  mutate(direction = sign(Estimate)) %>%
  mutate(Zscore = Zscore * direction)

max.limit <- 10 
p_color.co <- 0.25
p_label.co <- 0.05

res_BSI$PlotValue <- sapply(1:nrow(res_BSI),
                            function(i){
                              if(res_BSI$p.value[i] > p_color.co){
                                NA
                              }else{
                                #limits: 
                                if( abs(res_BSI$Zscore[i]) > max.limit) {
                                  max.limit*res_BSI$direction[i]
                                }else{
                                  res_BSI$Zscore[i]
                                }
                              }
                            })

res_BSI$Text <- sapply(1:nrow(res_BSI),
                       function(i){
                         if( res_BSI$p.value[i] > p_label.co){
                           NA
                         } else if(res_BSI$direction[[i]] == 1) {
                           "+"
                         }else if(res_BSI$direction[[i]] == -1) {
                           "-"
                         }
                       }) 

if(T){
  # cluster as a whole
  res_BSI.w <- res_BSI %>% reshape2::dcast(Species~dataSet, value.var = "Zscore") %>% tibble::column_to_rownames("Species")
  #res_BSI.w[is.na(res_BSI.w)] <- 0
  
  if(T){
    library(ggdendro)
    df <- res_BSI.w
    x <- as.matrix(scale(df))
    dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
    col.ord <- order.dendrogram(dd.col)
    
    dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
    row.ord <- order.dendrogram(dd.row)
    
    xx <- scale(df)[col.ord, row.ord] 
    xx_names <- attr(xx, "dimnames") 
    #df <- as.data.frame(xx)
    ddata_x <- dendro_data(dd.row) 
    ddata_y <- dendro_data(dd.col) 
  }
  
  xx_names[[1]]
  xx_names[[2]]
}

species.lvls <- (res_BSI %>%
                   group_by(Species) %>%
                   summarise(avgEstimate = mean(Estimate)) %>%
                   arrange(desc(avgEstimate)))$Species 

res_BSI$Species <- factor(res_BSI$Species, levels = species.lvls)
res_BSI$dataSet <- factor(res_BSI$dataSet, levels = rev(c("CAMEB2","EMBARC","118.stable")))

library(ggplot2)
p_BSI <- ggplot(res_BSI) +
  geom_tile(aes(x=Species, y=dataSet, fill = PlotValue), color="#E3E3E3") +
  scale_fill_gradient2(low = '#dd9196', high = '#4dabaf', mid = '#F5F5F4', midpoint = 0,na.value='white') +
  geom_text(aes(x=Species, y=dataSet, label=Text, color=Text)) +
  scale_color_manual(values = c("black","black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))
p_BSI


# AE_1year ------
res_AE_1year <- res_AE_1year %>% 
  mutate(Zscore = -qnorm(p_value/2) ) %>%
  mutate(direction = sign(Estimate)) %>%
  mutate(Zscore = Zscore * direction)

res_AE_1year$PlotValue <- sapply(1:nrow(res_AE_1year),
                                 function(i){
                                   if(res_AE_1year$p_value[i] > p_color.co){
                                     NA
                                   }else{
                                     #limits: 
                                     if( abs(res_AE_1year$Zscore[i]) > max.limit) {
                                       max.limit*res_AE_1year$direction[i]
                                     }else{
                                       res_AE_1year$Zscore[i]
                                     }
                                   }
                                 })

res_AE_1year$Text <- sapply(1:nrow(res_AE_1year),
                            function(i){
                              if( res_AE_1year$p_value[i] > p_label.co){
                                NA
                              } else if(res_AE_1year$direction[[i]] == 1) {
                                "+"
                              }else if(res_AE_1year$direction[[i]] == -1) {
                                "-"
                              }
                            }) 

res_AE_1year$Species <- factor(res_AE_1year$Species, levels = species.lvls)
res_AE_1year$dataSet <- factor(res_AE_1year$dataSet, levels = rev(c("CAMEB2","EMBARC","118.stable")))


p_AE_1year <- ggplot(res_AE_1year) +
  geom_tile(aes(x=Species, y=dataSet, fill = PlotValue), color="#E3E3E3") +
  scale_fill_gradient2(low = '#dd9196', high = '#4dabaf', mid = '#F5F5F4', midpoint = 0,na.value='white') +
  geom_text(aes(x=Species, y=dataSet, label=Text, color=Text)) +
  scale_color_manual(values = c("black","black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))

p_AE_1year

# FEV1pred ---------
res_FEV1pred <- res_FEV1pred %>% 
  mutate(Zscore = -qnorm(p.value/2) ) %>%
  mutate(direction = sign(Estimate)) %>%
  mutate(Zscore = Zscore * direction)

res_FEV1pred$PlotValue <- sapply(1:nrow(res_FEV1pred),
                                 function(i){
                                   if(res_FEV1pred$p.value[i] > p_color.co){
                                     NA
                                   }else{
                                     #limits: 
                                     if( abs(res_FEV1pred$Zscore[i]) > max.limit) {
                                       max.limit*res_FEV1pred$direction[i]
                                     }else{
                                       res_FEV1pred$Zscore[i]
                                     }
                                   }
                                 })

res_FEV1pred$Text <- sapply(1:nrow(res_FEV1pred),
                            function(i){
                              if( res_FEV1pred$p.value[i] > p_label.co){
                                NA
                              } else if(res_FEV1pred$direction[[i]] == 1) {
                                "+"
                              }else if(res_FEV1pred$direction[[i]] == -1) {
                                "-"
                              }
                            }) 

res_FEV1pred$Species <- factor(res_FEV1pred$Species, levels = species.lvls)
res_FEV1pred$dataSet <- factor(res_FEV1pred$dataSet, levels = rev(c("CAMEB2","EMBARC","118.stable")))

library(ggplot2)
p_FEV1pred <- ggplot(res_FEV1pred) +
  geom_tile(aes(x=Species, y=dataSet, fill = PlotValue), color="#E3E3E3") +
  scale_fill_gradient2(low = '#dd9196', high = '#4dabaf', mid = '#F5F5F4', midpoint = 0,na.value='white') +
  geom_text(aes(x=Species, y=dataSet, label=Text, color=Text)) +
  scale_color_manual(values = c("black","black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))
p_FEV1pred

# mMRC ---------
res_mMRC <- res_mMRC %>% 
  mutate(Zscore = -qnorm(p_value/2) ) %>%
  mutate(direction = sign(Estimate)) %>%
  mutate(Zscore = Zscore * direction)

res_mMRC$PlotValue <- sapply(1:nrow(res_mMRC),
                             function(i){
                               if( is.na(res_mMRC$p_value[i]) ){
                                 NA
                               } else if(res_mMRC$p_value[i] > p_color.co){
                                 NA
                               }else{
                                 #limits: 
                                 if( abs(res_mMRC$Zscore[i]) > max.limit) {
                                   max.limit*res_mMRC$direction[i]
                                 }else{
                                   res_mMRC$Zscore[i]
                                 }
                               }
                             })

res_mMRC$Text <- sapply(1:nrow(res_mMRC),
                        function(i){
                          if( is.na(res_mMRC$direction[i]) ){
                            NA
                          }else if( res_mMRC$p_value[i] > p_label.co){
                            NA
                          }else if(res_mMRC$direction[[i]] == 1) {
                            "+"
                          }else if(res_mMRC$direction[[i]] == -1) {
                            "-"
                          }
                        }) 

res_mMRC$Species <- factor(res_mMRC$Species, levels = species.lvls)
res_mMRC$dataSet <- factor(res_mMRC$dataSet, levels = rev(c("CAMEB2","EMBARC","118.stable")))

library(ggplot2)
p_MMRC <- ggplot(res_mMRC) +
  geom_tile(aes(x=Species, y=dataSet, fill = PlotValue), color="#E3E3E3") +
  scale_fill_gradient2(low = '#dd9196', high = '#4dabaf', mid = '#F5F5F4', midpoint = 0,na.value='white') +
  geom_text(aes(x=Species, y=dataSet, label=Text, color=Text)) +
  scale_color_manual(values = c("black","black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))
p_MMRC

library(ggpubr)

ggarrange(
  p_BSI + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
  p_FEV1pred + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
  p_AE_1year + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
  p_MMRC,
  ncol = 1, heights = c(0.15,0.15,0.15,0.55)
)
ggsave(filename = "S1g.heatmap_lm_topSpecies.pdf", width = 7, height = 5.5) 
