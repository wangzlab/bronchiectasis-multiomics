library(data.table)
library(dplyr)
species <- fread("../_data/combine_metagenome_rel(1).txt", data.table = F)
species$species <- sub(".*;s__","",species$`#NAME`)
test <- species %>% dplyr::select(species, `#NAME`)

species$`#NAME` <- NULL
species.t <- species %>% tibble::column_to_rownames("species") %>% t() %>% as.data.frame()


meta <- fread("../_data/combine_metadata(1).txt", data.table = F) 

# 由于要加上region，用之前的meta信息
meta.617 <- fread("../_data/combinedData_617sps.csv", data.table = F) 
meta$region <- 
  sapply(1:nrow(meta),
         function(i){
           if(meta$Study[i] == "China"){
             ""
           }else{
             meta.617$European_Region[which(meta.617$Sequencing.ID == meta$SampleID[i])]
           }
         })
table(meta$region )

# china的只取stable的
meta$Disease %>% unique()
table(meta$Study, meta$Disease)

meta.stable <- meta %>% filter(Disease == "Stable")

library(dplyr)
species.stable <- species.t[meta.stable$SampleID,]

res <- apply(species.stable, 2, mean) %>% as.data.frame()
speciesTop <-  rownames(res)[res$. > 1e-3]
#save(speciesTop, file = "1.topSpecies.1e-3.RData")

topSpecies.stable <- species.stable %>% dplyr::select(all_of(speciesTop))

groups <- meta.stable %>%
  rename(Cohort = Study, Region = region) %>%
  mutate(Group = paste0(Cohort,"-",Region))
groups <- groups[match(rownames(species.stable), groups$SampleID),]

# pcoa 
library(vegan)

speciesDat <- topSpecies.stable

Dist <- vegdist(speciesDat, method = "bray", binary = F)

pcoa <- cmdscale(Dist, k = (nrow(speciesDat) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 

# groups
site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)
site <- merge(site, groups, by.x = "name", by.y = "SampleID")

pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

library(ggplot2)


set1_colors <- setNames(c("#de9295","#d8585c","#ec8b4f","#88cad6","#3da1bf","#367583"),
                        nm = c("CAMEB2-Asia","CAMEB2-UK","China-","EMBARC-Northern_And_Western_Europe","EMBARC-Southern_Europe","CAMEB2-UK"))

ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(fill= Group), shape=21, size=3) +
  stat_ellipse(aes(color = Group, fill = Group),geom = 'polygon', alpha=0.1,
               level = 0.9) +   
  scale_color_manual(values=set1_colors)+
  scale_fill_manual(values=set1_colors)+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pcoa1, y = pcoa2, title = "PCoA")
ggsave(filename = "S1c.PCoA_stable_colorByCohort.region.pdf", width = 8, height = 4.5)




# 计算几个region之间的centroid distance -------------------
site %>% head()
# 计算每组质心
centroids <- site %>%
  group_by(Group) %>%
  summarise(
    Centroid1 = mean(X1),
    Centroid2 = mean(X2)
  )
# 计算质心间距离
centroid_coords <- as.matrix(centroids[, -1])
rownames(centroid_coords) <- centroids$Group
centroid_dist <- dist(centroid_coords)
# 转换为矩阵格式
centroid_dist_matrix <- as.matrix(centroid_dist)
print(centroid_dist_matrix)

centroid_dist_matrix[upper.tri(centroid_dist_matrix)] <- NA

plotD.dist <- as.data.frame(as.table(centroid_dist_matrix))
summary(plotD.dist$Freq)
ggplot(plotD.dist, 
       aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 3)), color = "white", size = 4) +
  scale_fill_gradient(low = "#deefe3",
                      high = "#005b64",
                      na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs( title = "Centroid Distances Between Groups")
ggsave(filename = "S1d.region_centroid.distance.pdf", width = 6, height = 4.5)
