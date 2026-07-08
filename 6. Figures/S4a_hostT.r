library(data.table)
library(dplyr)


meta.origin <- fread("../_dataGZ/group_all.txt", data.table = F)
#matches <- fread("metat_match.txt", data.table = F, header = F)
dat <- fread("../_dataGZ/hostt.txt",data.table = F)
dat <- dat %>% select(-EntrezID, -Name)
sapply(dat, class)
head(dat[,1:6]) 
colnames(dat) %in% meta.origin$`#NAME`
dat$`#NAME`

feature_df <-  
  dat %>%
  tibble::column_to_rownames("Symbols") %>%
  t() %>%
  as.data.frame(stringsAsFactors = F) 


meta <- meta.origin[match(rownames(feature_df), meta.origin$`#NAME`),]
newMicrobGrouping <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
meta$metaG.newGroup <- sapply(meta$`#NAME`,function(x) newMicrobGrouping$new_grouping_info_cutoff10[which(newMicrobGrouping$Sample == x)])
meta$grp.color <- 
  sapply(1:nrow(meta),
         function(i){
           if(meta$Type[i] == "Stable"){
             paste0(meta$Type[i],"-", meta$metaG.newGroup[i])
           }else{
             meta$Type[i]
           }
         })


sps.exacer <- meta$`#NAME`[meta$Type == "Exacerbation"]
sps.keep <- rownames(feature_df)[!rownames(feature_df) %in% sps.exacer]
feature_df <- feature_df[sps.keep,]
meta <- meta[match(rownames(feature_df), meta$`#NAME` ),]
meta$`#NAME` == rownames(feature_df)

if(any(sapply(feature_df,sum) == 0)) dat_positive <- feature_df[,-which(sapply(feature_df,sum) == 0)] else dat_positive <- feature_df

arg.pca <- prcomp(dat_positive, scale. = TRUE, center = T) 
summary(arg.pca,loadings = T)

arg.pca$rotation %>% head() # 每个otu的coordinates
arg.pca$x %>% head() # 每个sample的coordinates
sp_coordinates <-  as.data.frame(arg.pca$x[,1:2])
head(sp_coordinates)
sp_coordinates$Group <- sapply(rownames(sp_coordinates),
                               function(x) meta$Type[which(meta$`#NAME` == x)])
sp_coordinates$group_metaG <- sapply(rownames(sp_coordinates), function(x)meta$metaG.newGroup[which(meta$`#NAME` == x)] )
sp_coordinates$grp.color <- 
  sapply(1:nrow(sp_coordinates),
         function(i){
           if(sp_coordinates$Group[i] == "Stable"){
             paste0(sp_coordinates$Group[i],"-", sp_coordinates$group_metaG[i])
           }else{
             sp_coordinates$Group[i]
           }
         })


pca_exp <- summary(arg.pca,loadings = T)$importance[2,]
pca1 <- paste('PC1 :', round(100*pca_exp[1], 2), '%')
pca2 <- paste('PC2 :', round(100*pca_exp[2], 2), '%')


# 使用vegan包中的adonis函数进行PERMANOVA分析
library("vegan")
arg.adonis <- adonis2(dat_positive ~ grp.color, data = meta, permutations = 999)
arg.adonis$R2[1] ; arg.adonis$`Pr(>F)`[1]


# 之后在绘图代码中将PERMANVOA结果在PCA图中进行展示,同时修改图例
library(grid)
grob <- grobTree(textGrob(paste("adonis: R2 = ",round(arg.adonis$R2[1],3), " p-value = ", arg.adonis$`Pr(>F)`[1], sep = ""), 
                          x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=10)))

Colors <- setNames(c("#c9c7c7","#8abbd0","#75b18f","#d2898a","#deb56f"),
                   nm = c("Health","Stable-PPM","Stable-Pa","Stable-Hi","Stable-Commensal"))

P0 <- ggplot(sp_coordinates,aes(PC1,PC2))+
  geom_point(size=3, aes(fill = grp.color), shape=21) +
  stat_ellipse(aes(color = grp.color),  # 如果要加 polygon:  ", fill=group"
               #geom = 'polygon', alpha=0.1,
               level = 0.9) +   
  scale_color_manual(values=Colors)+
  scale_fill_manual(values=Colors)+
  #  scale_x_continuous(limits = c(-230,400)) + 
  #  scale_y_continuous(limits = c(-210,130)) + 
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
  annotation_custom(grob)+
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pca1, y = pca2)
P0
ggsave( filename = "S4a.hostT.pdf", P0, device = "pdf" , width = 5.5, height = 3.5)
