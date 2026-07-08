library(dplyr)
library(data.table)

species <- fread("../_data/combine_metagenome_rel(1).txt", data.table = F)
species$species <- sub(".*;s__","",species$`#NAME`)
test <- species %>% dplyr::select(species, `#NAME`)

species$`#NAME` <- NULL
species.t <- species %>% tibble::column_to_rownames("species") %>% t() %>% as.data.frame()

meta <- fread("../_data/combine_metadata(1).txt", data.table = F) 



meta.617 <- fread("../_data/combinedData_617sps.csv", data.table = F) 
meta$region <- 
  sapply(1:nrow(meta),
         function(i){
           if(meta$Study[i] == "China"){
             "China"
           }else{
             r = meta.617$European_Region[which(meta.617$Sequencing.ID == meta$SampleID[i])]
             paste0(meta$Study[i],"-", r)
           }
         })
table(meta$region )

microbGrouping <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
meta$microbGrp <- sapply(meta$SampleID, function(x) microbGrouping$new_grouping_info_cutoff10[which(microbGrouping$Sample == x)])

# china : stable only 
meta$Disease %>% unique()
table(meta$Study, meta$Disease)

meta.stable <- meta %>% filter(Disease == "Stable")


library(vegan)
shannon <- data.frame(Shannon= diversity(species.t, index = "shannon") )
#write.csv(shannon, file = "5.shannon_718sps.csv", quote = F, row.names = T)


all(meta.stable$SampleID %in% rownames(shannon))

plotD <- merge(meta.stable, shannon, by.x = 'SampleID', by.y=0)
plotD$Study <- factor(plotD$Study, levels = c("EMBARC", "CAMEB2", "China"))
head(plotD)


cohort.colors <- setNames(
  c("#efb7b7","#c4bed8","#88cbc0"),
  nm = c("CAMEB2","China","EMBARC") 
)


library(ggplot2)
library(ggpubr)

myCps <- combn(as.character(unique(plotD$Study)),m = 2, simplify = F)

ggboxplot(data = plotD, x="Study", y = "Shannon", fill = "Study", outlier.shape = NA) +
  geom_jitter(aes(x=Study, y=Shannon),  color="#c9c9c9", width = 0.25, alpha=0.5) +
  scale_fill_manual(values = cohort.colors) +
  stat_compare_means(comparisons= myCps) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = "S1b.Shannon_cohort.pdf", width = 5, height = 4)


# plot by region --------
unique(plotD$region)

region.colors <- setNames(
  c("#de9696","#d65351","#e89043","#96ccda","#40a0c1","#24707f"),
  nm = c("CAMEB2-Asia", "CAMEB2-UK", "China","EMBARC-Northern_And_Western_Europe","EMBARC-Southern_Europe","EMBARC-UK") 
)

plotD$region <- factor(plotD$region, levels = names(region.colors))

myCps2 <- combn(as.character(unique(plotD$region)),m = 2, simplify = F)

ggboxplot(data = plotD, x="region", y = "Shannon", fill = "region", outlier.shape = NA) +
  geom_jitter(aes(x=region, y=Shannon),  color="#c9c9c9", width = 0.25, alpha=0.5) +
  scale_fill_manual(values = region.colors) +
  stat_compare_means(comparisons= myCps2) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = "S1c.Shannon_region.pdf", width = 7, height = 7)
