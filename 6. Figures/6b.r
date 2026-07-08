library(data.table)
dat <- fread("All_Sample_BDI.txt", data.table = F)  
head(dat)

meta <- fread("../_data/combine_metadata(1).txt", data.table = F)
head(meta)

plotD.full <- merge(dat, meta, by = "SampleID")

# cluSters
meta.clusters <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
plotD.full$Cluster <- sapply(plotD.full$SampleID, function(x) meta.clusters$new_grouping_info_cutoff10[which(meta.clusters$Sample == x)])

library(dplyr)
library(ggplot2)
library(ggpubr)


# China Guangzhou -----------------------------------

pD2 <- plotD.full %>% filter(Study == "China") %>% filter(Disease %in% c("Stable","Healthy"))
pD2$Cluster2 <- ifelse(pD2$Disease == "Healthy", "Healthy", as.character(pD2$Cluster))
pD2$Cluster2 <- factor(pD2$Cluster2, levels = c("Pa","Hi","PPM","Commensal","Healthy"))
pairs2 <- combn(as.character(unique(pD2$Cluster2)),m=2,simplify = F)
ggboxplot(data = pD2, x="Cluster2", y="BDI", fill="Cluster2", outlier.shape = NA) +
  geom_jitter(color="#e4e4e4", alpha=0.5, width = 0.2)+
  scale_fill_manual(values = c("#8dc0a3","#e29192","#8cc2d8","#ecbf71","#c9c7c7")) +
  stat_compare_means(comparisons = pairs2,
                     label = "p.signif")
ggsave(filename = "6b.China.pdf", width = 3.3, height = 4.5)




# EMBARC-BRIDGE + Healthy --------------------------
pD4 <- plotD.full %>% filter(Study == "EMBARC")
pD4$Cluster2 <- ifelse(pD4$Disease == "Healthy", "Healthy", as.character(pD4$Cluster))
# healthy 
sp.healthy <- fread("88_German_healthy.txt", data.table = F, header = F)
dat$Sp_dd <- sub("Meta_(Sp_\\d+)_S.*$","\\1", dat$SampleID)
tmp <- dat %>% filter(Sp_dd %in% sp.healthy$V1) %>% mutate(Cluster2 = "Healthy", Disease = "Healthy") 
pD4 <- bind_rows(pD4, tmp)
pairs <- combn(as.character(unique(pD4$Cluster2)),m=2,simplify = F)
pD4$Cluster2 <- factor(pD4$Cluster2, levels = c("Pa","Hi","PPM","Commensal","Healthy"))
ggboxplot(data = pD4, x="Cluster2", y="BDI", fill="Cluster2", outlier.shape = NA) +
  geom_jitter(color="#e4e4e4", alpha=0.5, width = 0.2)+
  scale_fill_manual(values = c("#8dc0a3","#e29192","#8cc2d8","#ecbf71","#c9c7c7")) +
  stat_compare_means(comparisons = pairs,
                     label = "p.signif")
ggsave(filename = "6b.EMBARC_withHealthy.pdf", width = 3.3, height = 4.5)


# CAMEB2  -----------------------------------

pD3 <- plotD.full %>% filter(Study == "CAMEB2")
#pD3$Cluster <- factor(pD3$Cluster, levels = c("Pa","Hi","PPM","Commensal"))

CAMEB2.healthy <- fread("CI_CAMEB2_Healthy.txt", data.table = F)

pD3 <- bind_rows(
  pD3 %>% select(SampleID, BDI, Study, Disease, Cluster) %>% mutate(Cluster2 = Cluster),
  CAMEB2.healthy %>% rename(BDI = CI, Study = Cohort, Disease = Group) %>% mutate( Cluster2 = "Healthy") 
)

pairs3 <- combn(as.character(unique(pD3$Cluster2)),m=2,simplify = F)
pD3$Cluster2 <- factor(pD3$Cluster2, levels = c("Pa","Hi","PPM","Commensal","Healthy"))

ggboxplot(data = pD3, x="Cluster2", y="BDI", fill="Cluster2", outlier.shape = NA) +
  geom_jitter(color="#e4e4e4", alpha=0.5, width = 0.2)+
  scale_fill_manual(values = c("#8dc0a3","#e29192","#8cc2d8","#ecbf71","#c9c7c7")) +
  stat_compare_means(comparisons = pairs3,
                     label = "p.signif")
ggsave(filename = "6b.CAMEB2.pdf", width = 3.3, height = 4.5)


# full cohort with healthy -------------------------
pD1 <- plotD.full %>% filter(Disease == "Stable") %>% mutate(Cluster2 = Cluster)
# healthy 
tmp2 <- pD2 %>% filter(Disease == "Healthy") %>% select(SampleID, BDI, Study, Disease, Cluster, Cluster2)
tmp3 <- pD3 %>% filter(Disease == "Healthy") %>% select(SampleID, BDI, Study, Disease, Cluster,Cluster2) 
tmp4 <- pD4 %>% filter(Disease == "Healthy") %>% select(SampleID, BDI, Study, Disease, Cluster,Cluster2) 
pD1 <- bind_rows(pD1, tmp2, tmp3, tmp4) 
pairs <- combn(as.character(unique(pD1$Cluster2)),m=2,simplify = F)
pD1$Cluster2 <- factor(pD1$Cluster2, levels = c("Pa","Hi","PPM","Commensal","Healthy"))
table(pD1$Cluster2)
ggboxplot(data = pD1, x="Cluster2", y="BDI", fill="Cluster2", outlier.shape = NA) +
  geom_jitter(color="#e4e4e4", alpha=0.5, width = 0.2)+
  scale_fill_manual(values = c("#8dc0a3","#e29192","#8cc2d8","#ecbf71","#c9c7c7")) +
  stat_compare_means(comparisons = pairs,
                     label = "p.signif")
ggsave(filename = "6b.fullCohort_w.healthy.pdf", width = 3.3, height = 4.5)

