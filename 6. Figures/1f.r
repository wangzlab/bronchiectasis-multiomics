library(data.table)
library(dplyr)

list.files("../")

# meta --------------
vars <- c("FEV1pred_perc")
covars <- c("Age","Gender")


meta.617 <- fread("../_data/combinedData_617sps.csv", data.table = F) 
meta.617 <- meta.617 %>% 
  rename( FEV1pred_perc = FEV1_percent_predicted, SampleID=Sequencing.ID) %>% 
  dplyr::select(SampleID, Study, all_of(vars), all_of(covars))
sapply(meta.617, class)


meta.118 <- readxl::read_excel("../_data/metadata_118sps_more.xlsx")  
meta.118$Group %>% table()
meta.118.stable <- meta.118 %>% filter(Group == "Stable")  %>% dplyr::select(SampleID, all_of(vars), all_of(covars))
sapply(meta.118.stable, class)

meta <- bind_rows(
  meta.118.stable %>% mutate(Study = "China"),
  meta.617
)


# groups ------------------------------
ClusterID <- data.table::fread("../_data/new_grouping_cutoff10_v2.txt", data.table=F)
head(ClusterID)
colnames(ClusterID)[1:2] <- c("Sample", "Clusters")
ClusterID$Clusters %>% unique

# FEV1pred_perc ----------------------
pD.num <- sourceDat %>%
  dplyr::select(SampleID, Study, Clusters, FEV1pred_perc) %>%
  reshape2::melt(id.vars=c("SampleID","Study","Clusters"), variable.name="Variable") 
head(pD.num)

pD.num$Study <- factor(pD.num$Study, levels = c("China","CAMEB2","EMBARC"))
pD.num$Clusters <- factor(pD.num$Clusters, levels = c("Pa", "Hi", "PPM", "Commensal"))


Cluster.colors <- setNames(c("#7abf98","#8cc2d8","#ecbf71","#e29192"),
                           nm=c("Pa","PPM","Commensal","Hi"))

library(ggpubr)

mylist <- combn(as.character(unique(pD.num$Clusters)), m=2, simplify = F)

ggboxplot(pD.num, x="Clusters", y="value", outlier.shape = NA, fill = "Clusters" ) +
  stat_compare_means(comparisons = mylist, hide.ns = T, label="p.signif") +
  geom_jitter(color="gray", width = 0.2, alpha=0.5,  size=1) +
  scale_fill_manual(values = Cluster.colors) +
  facet_grid(Variable~Study,  scales = "free") +
  theme_bw()
ggsave(filename = "1f.FEV1_separateCohorts.pdf", width = 9, height = 5)

# fulldata
ggboxplot(pD.num, x="Clusters", y="value", outlier.shape = NA, fill = "Clusters" ) +
  stat_compare_means(comparisons = mylist, hide.ns = T, label="p.signif") +
  geom_jitter(color="gray", width = 0.2, alpha=0.5,  size=1) +
  scale_fill_manual(values = Cluster.colors) +
  facet_grid(Variable~.,  scales = "free") +
  theme_bw()
ggsave(filename = "1f.FEV1_fullData.pdf", width = 4.5, height = 5)
