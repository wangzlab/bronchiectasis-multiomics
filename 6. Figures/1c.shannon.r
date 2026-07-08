
library(data.table)
library(dplyr)

species <- fread("../_data/combine_metagenome_rel(1).txt", data.table = F)
species$species <- sub(".*;s__","",species$`#NAME`)
test <- species %>% dplyr::select(species, `#NAME`)


species$`#NAME` <- NULL
species.t <- species %>% tibble::column_to_rownames("species") %>% t() %>% as.data.frame()

meta <- fread("../_data/combine_metadata(1).txt", data.table = F) 



# microbe Clusters
microbGrouping <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
meta$microbGrp <- sapply(meta$SampleID, function(x) microbGrouping$new_grouping_info_cutoff10[which(microbGrouping$Sample == x)])


# stable only from China
meta$Disease %>% unique()
table(meta$Study, meta$Disease)

meta.stable <- meta %>% filter(Disease == "Stable")

# calculate Shannon
library(vegan)
shannon <- data.frame(Shannon= diversity(species.t, index = "shannon") )
write.csv(shannon, file = "1c.shannon_718sps.csv", quote = F, row.names = T)


all(meta.stable$SampleID %in% rownames(shannon))

plotD <- merge(meta.stable, shannon, by.x = 'SampleID', by.y=0)
plotD$microbGrp <- factor(plotD$microbGrp, levels = c("Pa","Hi","PPM","Commensal"))
head(plotD)

Cluster.colors <- setNames(
  c("#7abf98","#e29192","#8cc2d8", "#ebbe71"),
  nm = c("Pa","Hi","PPM","Commensal") 
)

library(ggplot2)
library(ggpubr)

myCps <- combn(as.character(unique(plotD$microbGrp)),m = 2, simplify = F)

ggboxplot(data = plotD, x="microbGrp", y = "Shannon", fill = "microbGrp", outlier.shape = NA) +
  geom_jitter(aes(x=microbGrp, y=Shannon),  color="#c9c9c9", width = 0.25, alpha=0.5) +
  scale_fill_manual(values = Cluster.colors) +
  stat_compare_means(comparisons= myCps)
ggsave(filename = "1c.Shannon.pdf", width = 5, height = 5)
