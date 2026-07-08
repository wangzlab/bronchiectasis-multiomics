library(data.table)
library(dplyr)

splist <- fread("sp_list.txt", data.table = F, header = F)

meta <- fread("../_dataGZ/group_all.txt",data.table = F)
sps.stable <- meta$`#NAME`[meta$Type =="Stable"]

dat.metaT <- fread("../_dataGZ/metaT_taxa.txt", data.table = F)
dat.metaT$`#NAME`

dat.metaT$Species <- sapply(strsplit(dat.metaT$`#NAME`, ";", fixed=T), "[[", 7)
dat.metaT$`#NAME` <- NULL

dat.metaT.t <- dat.metaT %>%
  tibble::column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()

rownames(dat.metaT.t) %in% meta$`#NAME`

dat.metaT.sub <- dat.metaT.t %>%
  select(all_of(splist$V7)) 
dat.metaT.sub <- dat.metaT.sub[sps.stable,]

library(Hmisc)

corRes <- rcorr(as.matrix(dat.metaT.sub), type = "spearman")
rho <- corRes$r
P <- corRes$P

#rho[upper.tri(rho)] <- NA
#P[upper.tri(P)] <- NA

rho.l <- 
  rho %>%
  reshape2::melt(value.name = "rho")
P.l <-
  P %>% 
  reshape2::melt(value.name = "P")

all(rho.l$Var1 == P.l$Var1)
all(rho.l$Var2 == P.l$Var2)

res <- cbind.data.frame(rho.l, "spearman.p" = P.l$P)
res.hi.pa <- res %>%
  mutate(Var2 = as.character(Var2)) %>%
  filter(Var2 %in% c("s__Haemophilus_influenzae","s__Pseudomonas_aeruginosa"))
write.csv(res.hi.pa, file = "4a.spearmanStable_HI.PA_speciesCorr.csv", quote = F, row.names = F)

# plotting ===================
dat.curves <- fread("4a.spearmanStable_HI.PA_speciesCorr.csv", data.table=F)

Hi.sps <- fread("Hi_assoc_sp_list.txt", data.table = F, header = F)
Hi.sps$V7 %in% dat.curves$Var1

Pa.sps <- fread("Pa_assoc_sp_list.txt", data.table = F, header = F)
Pa.sps$V7 %in% dat.curves$Var1

# data of metat taxonomy 
meta <- fread("../_dataGZ/group_all.txt",data.table = F)
sps.stable <- meta$`#NAME`[meta$Type =="Stable"]

dat.metaT <- fread("../_dataGZ/metaT_taxa.txt", data.table = F)
dat.metaT$`#NAME`

dat.metaT$Species <- sapply(strsplit(dat.metaT$`#NAME`, ";", fixed=T), "[[", 7)
dat.metaT$`#NAME` <- NULL

dat.metaT.t <- dat.metaT %>%
  tibble::column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()
dat.metat_stable <- dat.metaT.t[sps.stable, ]


library(ggdendro)
library(ggplot2)

hp.fullname <- setNames(c("s__Haemophilus_influenzae","s__Pseudomonas_aeruginosa"),
                        nm = c("Hi","Pa"))

for(hp in c("Hi","Pa")){
  # hp = "Pa"
  
  tmp <- eval(parse(text = paste0(hp,".sps")))
  species <- tmp$V7
  
  dat.metat_sub <- dat.metat_stable[,species]
  res <- rcorr(as.matrix(dat.metat_sub), type="spearman")
  r.full = res$r
  p.full = res$P
  
  if(T){
    df <- t(r.full) 
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
  
  assign(paste0("order_",hp),xx_names[[1]], envir = .GlobalEnv)
  
  r.l <- reshape2::melt(r.full, value.name = "rho")
  r.l$Var1 <- factor(r.l$Var1, levels = xx_names[[1]])
  r.l$Var2 <- factor(r.l$Var2, levels = xx_names[[1]])
  
  r.leveled <- r.l %>% reshape2::dcast(Var1~Var2) %>% tibble::column_to_rownames("Var1")
  
  r.leveled[upper.tri(r.leveled, diag =T)] <- NA  #先把顺序搞好了以后在去upper.tri赋值NA
  colnames(r.leveled) <- paste0(colnames(r.leveled),".x")  #为了区分画图时哪个是x，哪个是y
  rownames(r.leveled) <- paste0(rownames(r.leveled),".y")
  
  plotDat <- as.matrix(r.leveled) %>% reshape2::melt(value.name = "rho")
  colnames(plotDat)[1:2] <- c("Y","X")
  plotDat$Y <- factor(plotDat$Y, levels = rev(paste0(xx_names[[1]],".y"))) #y轴的levels要倒个方向否则是倒三角
  
   
  p.l <- p.full %>% reshape2::melt(value.name = "pval") %>% mutate(Var1=as.character(Var1), Var2=as.character(Var2))
  
  plotDat <- plotDat %>%
    mutate(species1 = sub(".x$","",X),
           species2 = sub(".y$", "", Y))
  
  plotDat <- merge(plotDat, p.l, by.x = c("species1","species2")  , by.y=c("Var1","Var2"))
  plotDat <- plotDat %>% mutate(rho2 = ifelse(pval<0.1, rho, NA))  #当pvalue < 0.1时才显示rho的颜色
  
  species.names <- levels(plotDat$Y)
  
  P1 <- ggplot(data = plotDat, aes(x=X,y=Y))+
    geom_tile(aes(fill=rho2),color="white") +
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_gradient2(low = "#E6949A", 
                         high = "#47B1B6", 
                         mid = "white",
                         midpoint = 0,
                         #   limits = c(-1000, 1000),
                         na.value="white") +
    theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) 
  
  P1
  
  # export results：
  write.table(plotDat, file = paste0("4a.corrHmPlotDat_",hp,".txt"), 
              sep = "\t", quote = F, row.names = F)
  
  plotDat_curve <-
    merge(plotDat %>% 
            filter(species1 == species2) %>%
            mutate(x.coord = as.integer(X),
                   y.coord = 0) %>%
            mutate(HP.x = length(species)/2, HP.y=1)  %>%
            select(-rho, -pval,-rho2),
          dat.curves %>% 
            filter(Var2 == hp.fullname[hp]),
          by.x = "species1", by.y="Var1") %>%
    arrange(x.coord) %>%
    mutate(curvation =  seq(0.2,-0.2, by=(-0.4)/(length(species)-1))) %>% # curvature要从0.2平滑过渡到-0.2
    mutate(signAssoci = as.character(sign(rho)))
  
  # export results：
  write.table(plotDat_curve, 
              file = paste0("4a.corrCurvePlotDat_",hp,".txt"), 
              sep = "\t", quote = F, row.names = F)
  
  
  # plotting
  curves_metab  <- lapply(1:nrow(plotDat_curve), function(i) {
    geom_curve(aes(x = plotDat_curve$HP.x[i], y = plotDat_curve$HP.y[i], 
                   xend = plotDat_curve$x.coord[i], yend = plotDat_curve$y.coord[i],  
                   color = plotDat_curve$signAssoci[i],
                   linewidth = abs(plotDat_curve$rho[i]), 
                   alpha=abs(plotDat_curve$rho[i])),
               #arrow = arrow(length = unit(0.3, "inches")),
               curvature = plotDat_curve$curvation[i]) 
  })
  
  assocLimits <- c(0.15, 0.61)  # association的最小值和最大值
  lineWid.range <- c(0.4, 1.2)  # association最小值和最大值 分别对应的线的最细和最粗值
  alpha.range <- c(0.2, 1)  # association最小值和最大值 分别对应的线的最浅和最深的透明度
  
  
  Colors = c("#47B1B6","#E6949A")
  names(Colors) <- c(1, -1)  
  
  
  P2 <- ggplot() + xlim(0,length(species)) + 
    ylim(0, 1) +   # 向下的curves
    # ylim(0, length(species)) +  # 向对角线的curves
    curves_metab +
    scale_x_continuous(breaks = 1:length(species), labels = plotDat_curve$species1) + # 修改x轴文本
    theme_void() +
    theme(legend.position = "bottom") +
    scale_linewidth(limits = assocLimits, range = lineWid.range) +  #手动调节curve粗细
    scale_alpha(limits = assocLimits, range = alpha.range) +  #调节curve的透明度
    scale_color_manual(values = Colors)
  P2
  
  
  #ggsave(P2, device = 'pdf', filename=paste0("4a.",hp,"_curves.pdf"), width = 8, height = 1.5)
  
  
  library(ggpubr)
  # pdf(paste0("1.",hp,".pdf"), width = 9, height = 4.5)
  P <- ggarrange(P1, ggarrange(P2, ggplot(), nrow=2,heights = c(0.35,0.65)), widths = c(0.6,0.4))
  # dev.off()
  
  ggsave(P, device = 'pdf', filename=paste0("4a.",hp,".pdf"), width = 13, height = 7)
}




# Bar plots =====================================================================
FC.Hi <- fread("Hi_assoc_sp_metagt.txt", data.table = F)
FC.Pa <- fread("Pa_assoc_sp_metagt.txt", data.table = F)

for(hp in c("Hi","Pa")){
  #hp = "Hi"
  FC = eval(parse(text = paste0("FC.",hp)))
  order_species <- eval(parse(text = paste0("order_",hp)))
  
  plotDat.bar <- 
    FC %>%
    select(V1, Commensal.metag.fc, PPM.metag.fc, Commensal.metag.p, PPM.metag.p) %>%
    mutate(species = sapply(strsplit(V1, ";", fixed = T),"[[", 7))  %>%
    mutate(logP.Commensal = (-log10(Commensal.metag.p))*sign(Commensal.metag.fc),
           logP.PPM = (-log10(PPM.metag.p))*sign(PPM.metag.fc)) 
  #  mutate(log2FC.Commensal = log2(Commensal.metag.fc), log2FC.PPM = log2(PPM.metag.fc))
  
  plotDat.bar$Taxa.lvled <- factor(plotDat.bar$species, levels = order_species)
  plotDat.bar$logP.Commensal <- ifelse(plotDat.bar$Commensal.metag.p < 0.1,plotDat.bar$logP.Commensal, NA )
  plotDat.bar$logP.PPM <- ifelse(plotDat.bar$PPM.metag.p < 0.1,plotDat.bar$logP.PPM, NA )
  
  npd.values <- plotDat.bar$logP.Commensal[!is.na(plotDat.bar$logP.Commensal)]
  assign(paste0("npd.vals_",hp), npd.values, envir = .GlobalEnv)
  
  ppm.values <- plotDat.bar$logP.PPM[!is.na(plotDat.bar$logP.PPM)]
  assign(paste0("ppm.vals_",hp), ppm.values, envir = .GlobalEnv)
}

# 跑完Hi和Pa后得到了 npd.vals_Hi, npd.vals_Pa, ppm.vals_Hi, ppm.vals_Pa
npd.values_full <- c(npd.vals_Hi, npd.vals_Pa)
ppm.values_full <- c(ppm.vals_Hi, ppm.vals_Pa)

values.full <- c(npd.values_full, ppm.values_full)

# 再各自画图：
for(hp in c("Hi","Pa")){
  #hp = "Hi"
  FC = eval(parse(text = paste0("FC.",hp)))
  order_species <- eval(parse(text = paste0("order_",hp)))
  
  plotDat.bar <- 
    FC %>%
    select(V1, Commensal.metag.fc, PPM.metag.fc, Commensal.metag.p, PPM.metag.p) %>%
    mutate(species = sapply(strsplit(V1, ";", fixed = T),"[[", 7))  %>%
    mutate(logP.Commensal = (-log10(Commensal.metag.p))*sign(Commensal.metag.fc),
           logP.PPM = (-log10(PPM.metag.p))*sign(PPM.metag.fc)) 
  #  mutate(log2FC.Commensal = log2(Commensal.metag.fc), log2FC.PPM = log2(PPM.metag.fc))
  
  plotDat.bar$Taxa.lvled <- factor(plotDat.bar$species, levels = order_species)
  plotDat.bar$logP.Commensal <- ifelse(plotDat.bar$Commensal.metag.p < 0.1,plotDat.bar$logP.Commensal, NA )
  plotDat.bar$logP.PPM <- ifelse(plotDat.bar$PPM.metag.p < 0.1,plotDat.bar$logP.PPM, NA )
  
  P_NPD <- ggplot(plotDat.bar) +
    geom_tile(aes(x=Taxa.lvled, y=0, fill=logP.Commensal), color="#E9E9E9") +
    scale_fill_gradientn(
      colors = c("#E6949A", "#E9E9E9","#47B1B6" ),  # 自定义颜色
      values = scales::rescale(c(min(min(values.full), 0)-0.1,
                                 0,  
                                 max(0,max(values.full))+0.1 )),  
      limits = c(min(min(values.full), 0)-0.1,   max(0,max(values.full))+0.1 ),
      na.value = "white"  # 设置 NA 值的颜色
    ) +
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    ggtitle(paste0(hp,"_Commensal_-logP*direction"))
  
  
  
  
  P_PPM <- ggplot(plotDat.bar) +
    geom_tile(aes(x=Taxa.lvled, y=0, fill=logP.PPM), color="#E9E9E9") +
    scale_fill_gradientn(
      colors = c("#E6949A", "#E9E9E9","#47B1B6" ),  # 自定义颜色
      values = scales::rescale(c(min(min(values.full), 0)-0.1,
                                 0,  
                                 max(0,max(values.full))+0.1 )),  
      limits = c(min(min(values.full), 0)-0.1,   max(0,max(values.full))+0.1 ),
      na.value = "white"  # 设置 NA 值的颜色
    ) +
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    ggtitle(paste0(hp,"_PPM_-logP*direction"))
  
  ggarrange(P_NPD, P_PPM)
  ggsave(filename = paste0("4a.bar_-logp.dire_",hp,".pdf"), width = 13, height = 4)
  
}


