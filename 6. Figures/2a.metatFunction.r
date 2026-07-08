
library(data.table)
library(dplyr)

meta.origin <- fread("../_dataGZ/group_all.txt", data.table = F)
matches <- fread("../_dataGZ/metat_match.txt", data.table = F, header = F)
dat <- fread("../_dataGZ/metat_ko.txt",data.table = F)
head(dat[,1:6]) 
colnames(dat) %in% matches$V1
dat$`#SampleID`


# calculate relative abundance
colSums(dat[,-1]) #每个colsum不一样

test <- sapply(dat[,-1], 
               function(x){
                 # print(class(x))
                 # print(length(x))
                 x/sum(x)
               }) %>% as.data.frame(stringsAsFactors=F) 
rownames(test) <- dat$`#SampleID`


feature_df <-  
  test %>%
  # tibble::column_to_rownames("#SampleID") %>%
  t() %>%
  as.data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column("Sample.id2") %>%
  mutate(sample = sapply(Sample.id2, function(x) matches$V2[which(matches$V1 == x)])) %>%
  select(-Sample.id2) %>%
  tibble::column_to_rownames("sample")

meta <- meta.origin[match(rownames(feature_df), meta.origin$`#NAME`),]
newMicrobGrouping <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)


meta$`#NAME` %in% newMicrobGrouping$Sample
meta$newGroup <- sapply(meta$`#NAME`,function(x) newMicrobGrouping$new_grouping_info_cutoff10[which(newMicrobGrouping$Sample == x)])
meta$grp.color <- 
  sapply(1:nrow(meta),
         function(i){
           if(meta$Type[i] == "Stable"){
             paste0(meta$Type[i],"-", meta$newGroup[i])
           }else{
             meta$Type[i]
           }
         })

sps.exacer <- meta$`#NAME`[meta$Type == "Exacerbation"]
sps.keep <- rownames(feature_df)[!rownames(feature_df) %in% sps.exacer]
feature_df <- feature_df[sps.keep,]
meta <- meta[match(rownames(feature_df), meta$`#NAME` ),]
meta$`#NAME` == rownames(feature_df)


# pcoa 
library(vegan)
Dist <- vegdist(feature_df, method = "bray", binary = F)

pcoa <- cmdscale(Dist, k = (nrow(feature_df) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 

# groups
site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)
site$group <- sapply(site$name, function(x) meta$Type[which(meta$`#NAME` == x)])
site$group_metaG <- sapply(site$name, function(x)meta$newGroup[which(meta$`#NAME` == x)] )
site$grp.color <- 
  sapply(1:nrow(site),
         function(i){
           if(site$group[i] == "Stable"){
             paste0(site$group[i],"-", site$group_metaG[i])
           }else{
             site$group[i]
           }
         })


head(site)
sapply(site, class)
table(site$grp.color)


pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

library("vegan")
arg.adonis <- adonis2(feature_df ~ grp.color, data = meta, permutations = 999)
arg.adonis$R2[1] ; arg.adonis$`Pr(>F)`[1]
'[1] 0.1712094
[1] 0.001'


library(grid)
grob1 <- grobTree(textGrob(paste("Adonis R2=",round(arg.adonis$R2[1], 3),
                                 ", p-value=",round(arg.adonis$`Pr(>F)`[1],3),
                                 sep = ""), 
                           x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=13)))


Colors <- setNames(c("#c9c7c7","#8cc2d8","#7abf98","#e29192","#ecbf71"),
                   nm = c("Health","Stable-PPM","Stable-Pa","Stable-Hi","Stable-Commensal"))


P0 <- ggplot(site,aes(X1,X2))+
  geom_point(size=3, aes(fill = grp.color), shape=21) +
  stat_ellipse(aes(color = grp.color),  # 如果要加 polygon:  ", fill=group"
               #geom = 'polygon', alpha=0.1,
               level = 0.9) +    
  scale_color_manual(values=Colors)+
  scale_fill_manual(values=Colors)+
  #scale_x_continuous(limits = c(-0.5,0.45)) +   # limits看了pcoa和pc1.density图后再回来改
  #scale_y_continuous(limits = c(-0.38,0.37)) +   # limits看了pcoa和pc2.density图后再回来改
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
  annotation_custom(grob1)+
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pcoa1, y = pcoa2)
P0
ggsave(filename = "2a.metatFunction.pdf", width = 5.5, height = 3.5)
