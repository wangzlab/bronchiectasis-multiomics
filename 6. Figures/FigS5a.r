library(ggdendro)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)

# plotting settings and necessary data match settings-------

# hm settings 
padj_co = 0.25
Colors <- c("#ADDCDE","#F3CFD2", "#47B1B6","#E6949A","#cc0202")  # "#cc0202"
names(Colors) <- c("Not-significant +","Not-significant -", 
                   "Significant positive","Significant negative",
                   "Significant + linked")

curveColors = c("#47B1B6","#E6949A")
names(curveColors) <- c(1, -1)  

# curve settings
#assocLimits <- c(0, 1)  # association的最小值和最大值
assocLimits_hostt <- c(0,0.8)
assocLimits_metat <- c(0,1)
lineWid.range <- c(0.5, 1.5)  # association最小值和最大值 分别对应的线的最细和最粗值
alpha.range <- c(0.2, 1)  # association最小值和最大值 分别对应的线的最浅和最深的透明度

# match list
sp_cp_matchlist <- list(
  Pm=c("Prevotella_melaninogenica","Pa_kraken"),
  Vp=c("Veillonella_parvula","Pa_kraken"),
  Rm = c("Rothia_mucilaginosa","Hi_kraken"),
  Vd = c("Vd","Pa_kraken")
)

limits_df_all <- NULL
CurveLimits_df_all <- NULL

sp_cps <- setNames(c("Pa","Pa","Hi","Pa"),
                   c("Pm", "Vp", "Rm","Vd"))

# read data =======================

for(sp in names(sp_cps)){
  # sp = "Pm"
  sp_cp <- sp_cps[sp]
  
  # files ----
  Files <- list.files("FigS5a", full.names = T)
  
  selF_metat <- Files[grepl(sp, Files) & grepl("metat_select", Files)]  #"multi-omic_heatmap_update_20251231/Pm_metat_select"
  selF_metab <- Files[grepl(sp, Files) & grepl("metab_select", Files)]  #"multi-omic_heatmap_update_20251231/Pm_metab_select"
  selF_hostt <- Files[grepl(sp, Files) & grepl("hostt_select", Files)]  #"multi-omic_heatmap_update_20251231/Pm_hostt_select"
  
  
  # -----
  if(F){
    curveF_exac_sp <- Files[grepl(sp, Files) & grepl("metat_exac_asso.txt", Files)]  # "data/Pi_metat_exac_asso.txt" 
    curveF_exac_Pa <- "data/Pa_Hi_gsva_exac_asso.txt" # same for all species
    curveF_exac_Hi <- "data/Pa_Hi_gsva_exac_asso.txt" # same for all species
    
    assoF_exac_metab_metat <- Files[grepl(sp, Files) & grepl("metab_exac_asso.txt", Files)]  #"data/Pi_metab_exac_asso.txt"
    linkF_exac_metab_metat <- Files[grepl(sp, Files) & grepl("metab_exac_link.txt", Files)]  #"data/Pi_metab_exac_link.txt"
    
    assoF_exac_metab_hostt <- "data/metab_hostt_exac_asso.txt"  # same for all species
    linkF_exac_metab_hostt <- "data/metab_hostt_exac_link.txt"  # same for all species
    
    
  }
  
  # -----
  if(F){
    curveF_health_sp <- Files[grepl(sp, Files) & grepl("metat_health_asso.txt", Files)]  # "data/Pi_metat_health_asso.txt"  
    curveF_health_Pa <- "data/Pa_Hi_gsva_health_asso.txt"   # same for all species
    curveF_health_Hi <- "data/Pa_Hi_gsva_health_asso.txt"   # same for all species
    
    assoF_health_metab_metat <- Files[grepl(sp, Files) & grepl("metab_health_asso.txt", Files)]  # "data/Pi_metab_health_asso.txt"
    linkF_health_metab_metat <- Files[grepl(sp, Files) & grepl("metab_health_link.txt", Files)]  # "data/Pi_metab_health_link.txt" 
    
    assoF_health_metab_hostt <- "data/metab_hostt_health_asso.txt"  # same for all species
    linkF_health_metab_hostt <- "data/metab_hostt_health_link.txt"   # same for all species
    
    
  }
  
  # -----
  
  curveF_stable_sp <- Files[grepl(sp, Files) & grepl("metat_assoc.txt", Files)]  #"data/Pi_metat_stable_asso.txt" 
  curveF_stable_Pa <- "data/Pa_Hi_gsva_stable_asso.txt"  # same for all species
  curveF_stable_Hi <- "data/Pa_Hi_gsva_stable_asso.txt"  # same for all species
  
  assoF_stable_metab_metat <- Files[grepl(sp, Files) & grepl("metat_metab_assoc", Files)]  #"multi-omic_heatmap_update_20251231/Pm_metat_metab_assoc.txt"
  linkF_stable_metab_metat <- Files[grepl(sp, Files) & grepl("metat_metab_link", Files)]  #"multi-omic_heatmap_update_20251231/Pm_metat_metab_link.txt"
  
  assoF_stable_metab_hostt <- "multi-omic_heatmap_update_20251231/metab_hostt_assoc.txt"  # same for all species
  linkF_stable_metab_hostt <- "multi-omic_heatmap_update_20251231/metab_hostt_link.txt"  # same for all species
  
  # read data -------
  sel.metat <- fread(selF_metat, data.table = F, header = F)
  sel.metab <- fread(selF_metab, data.table = F, header = F)
  sel.hostt <- fread(selF_hostt, data.table = F, header = F)
  
  
  for(stg in c("stable")){
    # stg = "stable"
    
    assoF_left <- eval(parse(text = paste0("assoF_",stg,"_metab_metat")))
    asso_left <- fread(assoF_left, data.table = F)
    linkF_left <- eval(parse(text = paste0("linkF_", stg, "_metab_metat"))) 
    link_left <- fread(linkF_left, data.table = F)
    
    assoF_right <- eval(parse(text = paste0("assoF_", stg, "_metab_hostt")))
    asso_right <- fread(assoF_right, data.table = F)
    linkF_right <- eval(parse(text = paste0("linkF_",stg,"_metab_hostt")))
    link_right <-  fread(linkF_right, data.table = F)
    
    curveF_left <- eval(parse(text = paste0("curveF_",stg,"_sp")))
    curve_left <- fread(curveF_left, data.table = F)
    curveF_right <- eval(parse(text = paste0("curveF_",stg,"_",sp_cp)))
    curve_right <- fread(curveF_right, data.table = F)
    
    
    colnames(asso_left)[1:2] <- c("metat","metab")
    colnames(link_left)[1:2] <- c("metab","metat")
    
    colnames(asso_right)[1:2] <- c("metab","hostt")
    colnames(link_right)[1:2] <- c("metab","hostt")
    
    colnames(curve_left)[1:2] <- c("species","metat")
    colnames(curve_right)[1:2] <- c("species","hostt")
    
    
    # left panel: plotDat1 -----
    plotDat1 <- asso_left %>%
      filter(metat %in% sel.metat$V1) %>%
      filter(metab %in% sel.metab$V1) 
    
    
    plotDat1 <- merge(plotDat1,
                      link_left %>% mutate(linked="linked") %>% select(metab, metat, linked) %>% unique(),
                      by=c("metab","metat"), all.x = T)
    plotDat1$linked[is.na(plotDat1$linked)] <- "not linked"
    
    
    plotDat1$color <- 
      sapply(1:nrow(plotDat1),
             function(i){
               if(plotDat1$`p-values`[i] <= padj_co & plotDat1$linked[i] =="linked"){
                 "Significant + linked"
               }else if(plotDat1$`p-values`[i] <= padj_co & plotDat1$association[i] > 0){
                 "Significant positive"
               }else if(plotDat1$`p-values`[i] <= padj_co & plotDat1$association[i] < 0){
                 "Significant negative"
               }else if(plotDat1$`p-values`[i] > padj_co & plotDat1$association[i] > 0) {
                 "Not-significant +"
               }else{
                 "Not-significant -"
               }
             })
    
    # clustering 
    
    if( stg == "stable") {
      dat.w <- asso_left %>% 
        reshape2::dcast(metab~metat, value.var = "association") %>% 
        filter(metab %in% sel.metab$V1) %>%
        select(metab, all_of(sel.metat$V1))
      rownames(dat.w) <- dat.w$metab; dat.w <- dat.w[-1]
      
      if(T){
        df <- t(dat.w) 
        x <- as.matrix(scale(df))
        dd.col <- as.dendrogram(hclust(dist(x), method = "average"))
        col.ord <- order.dendrogram(dd.col)
        
        dd.row <- as.dendrogram(hclust(dist(t(x)), method = "average"))
        row.ord <- order.dendrogram(dd.row)
        
        xx <- scale(df)[col.ord, row.ord] 
        xx_names <- attr(xx, "dimnames") 
        #df <- as.data.frame(xx)
        ddata_x <- dendro_data(dd.row) 
        ddata_y <- dendro_data(dd.col) 
      }   
      order_metab <- xx_names[[2]]
      order_metat <- xx_names[[1]]
    }
    
    plotDat1$metat <- factor(plotDat1$metat, levels = order_metat) 
    plotDat1$metab <- factor(plotDat1$metab, levels = order_metab)  
    
    P1 <- ggplot(data = plotDat1, aes(x=metat,y=metab)) +
      geom_tile(aes(fill=color, alpha=abs(association)),color="white") +
      theme(axis.text.x = element_text(angle = 90))+
      #scale_fill_manual(values=c("white","#e2e2e2","#cc0202"))+
      scale_fill_manual(values = Colors) +
      scale_alpha(limits = c(0, 1), range = c(0.5,1))+
      theme(  panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())
    # axis.title.x = element_text(colour=NA),
    # axis.title.y = element_blank())
    
    P1
    
    limits_df <- 
      data.frame(species=sp,
                 stage=stg,
                 FtrPair = paste(colnames(asso_left)[1:2], collapse = "_vs._"),
                 min.abs.assoc = min(abs(plotDat1$association)),
                 max.abs.assoc = max(abs(plotDat1$association)))
    limits_df_all <- bind_rows(limits_df_all, limits_df)
    
    # right panel: plotDat2 -----
    plotDat2 <- asso_right %>%
      filter(hostt %in% sel.hostt$V1) %>%
      filter(metab %in% sel.metab$V1) 
    
    plotDat2 <- merge(plotDat2,
                      link_right %>% mutate(linked="linked") %>% select(metab, hostt, linked) %>% unique(),
                      by=c("metab","hostt"), all.x = T)
    plotDat2$linked[is.na(plotDat2$linked)] <- "not linked"
    
    
    plotDat2$color <- 
      sapply(1:nrow(plotDat2),
             function(i){
               if(plotDat2$`p-values`[i] <= padj_co & plotDat2$linked[i] =="linked"){
                 "Significant + linked"
               }else if(plotDat2$`p-values`[i] <= padj_co & plotDat2$association[i] > 0){
                 "Significant positive"
               }else if(plotDat2$`p-values`[i] <= padj_co & plotDat2$association[i] < 0){
                 "Significant negative"
               }else if(plotDat2$`p-values`[i] > padj_co & plotDat2$association[i] > 0) {
                 "Not-significant +"
               }else{
                 "Not-significant -"
               }
             })
    
    
    
    # clustering 
    if( stg == "stable") {
      dat.w <- asso_right %>% 
        reshape2::dcast(metab~hostt, value.var = "association") %>% 
        filter(metab %in% sel.metab$V1) %>%
        select(metab, all_of(sel.hostt$V1))
      rownames(dat.w) <- dat.w$metab; dat.w <- dat.w[-1]
      
      
      if(T){
        df <- t(dat.w) 
        x <- as.matrix(scale(df))
        dd.col <- as.dendrogram(hclust(dist(x), method = "average"))
        col.ord <- order.dendrogram(dd.col)
        
        dd.row <- as.dendrogram(hclust(dist(t(x)), method = "average"))
        row.ord <- order.dendrogram(dd.row)
        
        xx <- scale(df)[col.ord, row.ord] 
        xx_names <- attr(xx, "dimnames") 
        #df <- as.data.frame(xx)
        ddata_x <- dendro_data(dd.row) 
        ddata_y <- dendro_data(dd.col) 
      } 
      order_hostt <- xx_names[[1]]
    }
    
    plotDat2$hostt <- factor(plotDat2$hostt, levels = order_hostt) 
    plotDat2$metab <- factor(plotDat2$metab, levels = order_metab)  
    
    
    P2 <- ggplot(data = plotDat2, aes(x=hostt,y=metab))+
      geom_tile(aes(fill=color, alpha=abs(association)),color="white") +
      theme(axis.text.x = element_text(angle = 90))+
      #scale_fill_manual(values=c("white","#e2e2e2","#cc0202"))+
      scale_fill_manual(values = Colors) +
      scale_alpha(limits = c(0, 1), range = c(0.5,1))+
      theme(  panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())
    # axis.title.x = element_text(colour=NA),
    # axis.title.y = element_blank())
    
    P2
    
    limits_df <- 
      data.frame(species=sp_cp,
                 stage=stg,
                 FtrPair = paste(colnames(asso_right)[1:2], collapse = "_vs._"),
                 min.abs.assoc = min(abs(plotDat2$association)),
                 max.abs.assoc = max(abs(plotDat2$association)))
    limits_df_all <- bind_rows(limits_df_all, limits_df)
    
    
    # curve left : P3 ====================
    
    plotDat3 <- curve_left %>% 
      filter(species == sp_cp_matchlist[[sp]][1]) %>%
      filter(metat %in% order_metat) %>%
      mutate(metat = factor(metat, levels=order_metat)) %>%
      arrange(metat) %>%
      mutate(x2 = 1:length(order_metat)) %>% # 20个metat module
      mutate(x1 = (1+length(order_metat))/2) %>%
      mutate(y1= 1, y2 = 0) %>%
      mutate(signAssoci = as.character(sign(association))) %>%
      mutate(curvation =  seq(0.2,-0.2, by=(-0.4)/(length(order_metat)-1)))  # curvature要从0.2平滑过渡到-0.2
    
    limits_df <- 
      data.frame(species=sp,
                 stage=stg,
                 FtrPair = paste(colnames(curve_left)[1:2], collapse = "_vs._"),
                 min.abs.assoc = min(abs(plotDat3$association)),
                 max.abs.assoc = max(abs(plotDat3$association)))
    CurveLimits_df_all <- bind_rows(CurveLimits_df_all, limits_df)
    
    # 逐个绘制曲线，并为每条曲线指定不同的曲率值
    curves_metat <- lapply(1:nrow(plotDat3), function(i) {
      geom_curve(aes(x = plotDat3$x1[i], y = plotDat3$y1[i], xend = plotDat3$x2[i], yend = plotDat3$y2[i],  
                     color = plotDat3$signAssoci[i],
                     linewidth = abs(plotDat3$association[i]), 
                     alpha=abs(plotDat3$association[i])),
                 #arrow = arrow(length = unit(0.3, "inches")),
                 curvature = plotDat3$curvation[i]) 
    })
    
    
    # 添加所有曲线到图表中
    p3_metat <- ggplot() + xlim(0,length(order_metat)) + ylim(0, 1.5) +
      curves_metat +
      scale_x_continuous(breaks = 1:length(order_metat), labels = plotDat3$metat) + # 修改x轴文本
      theme_void() +
      theme(axis.text.x = element_text(angle = 90))  +
      scale_linewidth(limits = assocLimits_metat, range = lineWid.range) +  #手动调节curve粗细
      scale_alpha(limits = assocLimits_metat, range = alpha.range) +  #调节curve的透明度
      scale_color_manual(values = curveColors)
    p3_metat
    
    
    # P4: hostt ================================
    
    plotDat4 <- curve_right %>% 
      filter(species == sp_cp_matchlist[[sp]][2]) %>%
      filter(hostt %in% order_hostt) %>%
      mutate(hostt = factor(hostt, levels=order_hostt)) %>%
      arrange(hostt) %>%
      mutate(x2 = 1:length(order_hostt)) %>% # 20个hostt module
      mutate(x1 = (1+length(order_hostt))/2) %>%
      mutate(y1= 1, y2 = 0) %>%
      mutate(signAssoci = as.character(sign(association))) %>%  # as.character是因为要aes(color), scale_color_manual的时候定义的是分类变量
      mutate(curvation =  seq(0.2,-0.2, by=(-0.4)/(length(order_hostt)-1)))  # curvature要从0.2平滑过渡到-0.2
    
    limits_df <- 
      data.frame(species=sp_cp_matchlist[[sp]][2],
                 stage=stg,
                 FtrPair = paste(colnames(curve_right)[1:2], collapse = "_vs._"),
                 min.abs.assoc = min(abs(plotDat4$association)),
                 max.abs.assoc = max(abs(plotDat4$association)))
    CurveLimits_df_all <- bind_rows(CurveLimits_df_all, limits_df)
    
    
    
    # 逐个绘制曲线，并为每条曲线指定不同的曲率值
    curves_hostt <- lapply(1:nrow(plotDat4), function(i) {
      geom_curve(aes(x = plotDat4$x1[i], y = plotDat4$y1[i], xend = plotDat4$x2[i], yend = plotDat4$y2[i],  
                     color = plotDat4$signAssoci[i],
                     linewidth = abs(plotDat4$association[i]), 
                     alpha=abs(plotDat4$association[i])),
                 #arrow = arrow(length = unit(0.3, "inches")),
                 curvature = plotDat4$curvation[i]) 
    })
    
    
    # 添加所有曲线到图表中
    p4_hostt <- ggplot() + xlim(0,length(order_hostt)) + ylim(0, 1.5) +
      curves_hostt +
      scale_x_continuous(breaks = 1:length(order_hostt), labels = plotDat4$hostt) + # 修改x轴文本
      theme_void() +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_linewidth(limits = assocLimits_hostt, range = lineWid.range) +  #手动调节curve粗细
      scale_alpha(limits = assocLimits_hostt, range = alpha.range) +  #调节curve的透明度
      scale_color_manual(values = curveColors)
    p4_hostt
    
    
    P_corrHms <- ggarrange(P1 + theme(legend.position = "none"), 
                           P2+theme(axis.text.y = element_blank()),nrow = 1 )
    #ggsave(P_corrHms, filename = paste0(sp, "_", stg, "_corrHm_select2.pdf"), width = 10, height = 4)
    ggsave(P_corrHms, filename = paste0("select_20251231/", sp, "_corrHm.pdf"), width = 10, height = 4)
    
    P_curves <- ggarrange(p3_metat + theme(legend.position = "none"),
                          p4_hostt, heights = c(0.8,1), widths = c(0.385,0.615))
    #ggsave(P_curves,  filename = paste0(sp, "_", stg, "_curves_select2.pdf"), width = 10, height = 3)
    ggsave(P_curves,  filename = paste0("select_20251231/", sp, "_curves.pdf"), width = 10, height = 3)
    
  }# loop through health, stable and exac
} # loop through species

#save(limits_df_all, CurveLimits_df_all, 
#     file = "limits.txt")