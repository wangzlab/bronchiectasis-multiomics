library(data.table)
library(dplyr)

meta <- fread("../_data/combine_metadata(1).txt", data.table = F) 
head(meta)

species <- fread("../_data/combine_metagenome_rel(1).txt", data.table = F)
colnames(species)[!colnames(species) %in% meta$SampleID]

species.l <-
  species %>% 
  reshape2::melt(id.vars = "#NAME", variable.name = "Sample")
head(species.l)
colnames(species.l)[colnames(species.l) == "#NAME"] <- "Species"

# species to show 
species.l <- merge(species.l, meta %>% select(SampleID, Disease), by.x = 'Sample', by.y = "SampleID")

topSpecies <- species.l %>%
  group_by(Disease, Species) %>%
  summarise(avg = mean(value)) %>%
  ungroup %>%
  group_by(Disease) %>%
  slice_max(order_by = avg, n = 10) 
spcs2cluster <- unique(topSpecies$Species) %>% as.character() 


species.l.4cluster <- species.l %>%
  mutate(Species = ifelse(grepl("Neisseria", as.character(Species)), "Neisseria", as.character(Species)))  %>%
  mutate(Species.simp = ifelse(Species %in% c(spcs2cluster,"Neisseria"), Species, "others")) %>%
  group_by(Sample,Species.simp ) %>%
  summarise(value=sum(value)) %>%
  ungroup 


z <- species.l.4cluster %>%
  reshape2::dcast(Sample ~ Species.simp , value.var = "value") %>%
  tibble::column_to_rownames("Sample")%>%
  as.matrix()


# normalizaiton
means <- apply(z,2,mean)
sds <- apply(z,2,sd)
nor <- scale(z,center=means,scale=sds) 
apply(nor, 2, mean); apply(nor, 2, sd)

mtd="ward.D"  #   "ward.D"
mydata.hclust <- hclust(dist(nor),method=mtd)  
#plot(mydata.hclust, hang=-1, labels=rownames(z),main=mtd)
dd.row <- as.dendrogram(mydata.hclust)
od <- order.dendrogram(dd.row) 


# Cluster membership
if(T){
  member = cutree(mydata.hclust, 7)
  table(member)
  member #每个样品属于哪个group
  names(which(member == 1))
  tmp <- aggregate(nor,list(member), mean)  %>% data.frame() %>% t() # Characterizing clusters
  colnames(tmp) <-  tmp[1,]
  member.df <- tmp[-1,]
  
  afterClusterID <- as.data.frame(member) %>% tibble::rownames_to_column("SampleID")
}


# plotting --------------------
library(ggplot2)
library(ggdendro)

# ggplot to plot clust segment figure: p2
ddata_x <- dendro_data(dd.row) 
p2<- ggplot(segment(ddata_x)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_dendro() + theme(axis.title.x=element_blank())

xx <- scale(t(z))[, od] 
xx_names <- attr(xx, "dimnames")
xx_names[[2]]

# further simplify the species to show
topSpecies <- species.l %>%
  group_by(Disease, Species) %>%
  summarise(avg = mean(value)) %>%
  ungroup %>%
  group_by(Disease) %>%
  slice_max(order_by = avg, n = 8) 
spcs2show <- unique(topSpecies$Species) %>% as.character() 


species.l.simp <- species.l %>%
  mutate(Species.simp = ifelse(Species %in% spcs2show, as.character(Species), "others")) %>%
  group_by(Sample,Species.simp ) %>%
  summarise(value=sum(value)) %>%
  ungroup

plotDat <- species.l.simp %>% rename(Species = Species.simp)
plotDat$Sample <- factor(plotDat$Sample, levels = xx_names[[2]])

plotDat$Species <- as.character(plotDat$Species)
plotDat$Species <- sub(".*;s__","", plotDat$Species)
spcs2show <- sub(".*;s__","", spcs2show )
spcs2show <- spcs2show[order(spcs2show)]
spcs2show <- c( "others",spcs2show[spcs2show != "others"])  

plotDat$Species <- factor(plotDat$Species, levels = spcs2show)

# colors of species 
spcs2show <- spcs2show[order(spcs2show)]
spcs2show <- c(spcs2show[spcs2show != "others"], "others")

species.colors <- setNames(
  c("#e29192","#a4559d","#7abf98",
    "#25528D","#3ca6de","#a9dbed",
    "#1ea239","#C55102","#167968","#f8cedb",
    "#efbc1b","#dac475","#f7f39a","gray",
    "#B790D3","#664CFF","#ed7117","#b88d18","#f3a649", "#69869e"),
  nm = c("Haemophilus_influenzae","Haemophilus_parainfluenzae","Pseudomonas_aeruginosa",
         "Neisseria flavescens","Neisseria_mucosa","Neisseria_subflava",
         "Escherichia_coli", "Salmonella_enterica","Moraxella_catarrhalis","Rothia_mucilaginosa",
         "Streptococcus_mitis", "Streptococcus_pneumoniae", "Staphylococcus_aureus","others",
         "Veillonella_atypica", "Veillonella_dispar","Prevotella_melaninogenica","Prevotella_jejuni","Prevotella_intermedia","Schaalia_odontolytica" )
)

names(species.colors)[!names(species.colors) %in% spcs2show]
spcs2show[!spcs2show %in% names(species.colors)]


library(RColorBrewer)

p0<-ggplot(plotDat, aes(x = Sample, y = value, fill = Species)) + 
  geom_col() +  
  #scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(length(spcs2show))) +
  scale_fill_manual(values = species.colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  #facet_grid(cols=vars(sampleTypes),scales = "free_x", space = "free_x") + 
  xlab("") + ylab("Relative abundance") + 
  guides(fill = guide_legend(nrow = 4,
                             label.theme = element_text(size = 8))) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
p0


# p3: color of clusters 
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank()
)

afterClusterID$Sample <- factor(afterClusterID$Sample, levels = xx_names[[2]])
afterClusterID <- afterClusterID %>% arrange(Sample)
afterClusterID$member %>% table
write.csv(afterClusterID %>% select(-SampleID), file = "1b.sampleSeqs.csv", quote = F, row.names = F)

p3 <- ggplot(afterClusterID) +
  geom_tile(aes(x=Sample, y=0, fill = as.character(member))) +
  # geom_text(aes(x=Sample, y=0, label=label)) +
  #scale_fill_manual(values = Cluster.colors)+
  theme_none + theme(legend.position = "none")
p3


library(grid)
pdf(paste0("1b.cluster_v1.",mtd,".pdf"), width = 12, height = 7)
grid.newpage()
print(p3 + theme(legend.position = "none"), vp=viewport(width = 0.92, height = 0.4, x=0.5, y = 0.72 ))
print(p2, vp=viewport(width=1, height=0.2, x=0.5, y=0.87))
print( p0+theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()), 
       vp=viewport(width=0.92, height=0.7, x=0.5, y=0.36))
dev.off()


