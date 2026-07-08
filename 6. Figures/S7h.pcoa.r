library(data.table)
library(dplyr)
dat <- fread("all_combined_sp.txt", data.table = F)
dat.t <- dat %>% 
  tibble::column_to_rownames("#SampleID") %>%
  t() %>%
  as.data.frame()

apply(dat.t, 1, sum)

meta <- fread("all_combined_meta.txt", data.table = F)
meta <- meta %>% filter(!Group %in% c("Asthma2" , "Exacerbations"))

# make sure data and groups have the same amount of samples 
common.sps <- intersect(rownames(dat.t), meta$`#SampleID`)

dat.t <- dat.t[match(common.sps, row.names(dat.t)),]
group_df <- meta[match(common.sps, meta$`#SampleID`),]



# pcoa 
library(vegan)
Dist <- vegdist(dat.t, method = "bray", binary = F)

pcoa <- cmdscale(Dist, k = (nrow(dat.t) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 


# groups
site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)
site$group <- sapply(site$name, function(x) group_df$Group[which(group_df$`#SampleID` == x)])


#前 2 轴解释量
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')


# adonis
library("vegan")
arg.adonis <- adonis2(as.matrix(dat.t) ~ Group, data = meta, permutations = 999)
arg.adonis$R2[1] ; arg.adonis$`Pr(>F)`[1]




#ggplot2 作图
library(grid)
grob1 <- grobTree(textGrob(paste("Adonis.r2=",round(arg.adonis$R2[1], 3),
                                 ", pvalue=",round(arg.adonis$`Pr(>F)`[1],3),
                                 sep = ""), 
                           x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=13)))

library(RColorBrewer)
library(ggplot2)

site$group <- factor(site$group, levels = c("Healthy","BC","ABO","CP","Asthma","Stable"))

Colors <- setNames(c("#836daf","#cea8a6","#ceb8c3","#a2d1b5","#abcfdd","#c95462"),
                   nm=c("Healthy","BC","ABO","CP","Asthma","Stable"))


p <- ggplot(data = site, aes(X1, X2)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.2, show.legend = TRUE) +    #添加置信椭圆，注意不是聚类
  geom_point(aes(color = group),size=2.5) +
  geom_point(size=2.5,shape=21) +
  scale_color_manual(values = Colors) +
  scale_fill_manual(values = Colors) +
  annotation_custom(grob1)+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pcoa1, y = pcoa2, title = "PCoA")

p
ggsave(filename = "S7h.pcoa.pdf", width = 5.5, height = 4)
