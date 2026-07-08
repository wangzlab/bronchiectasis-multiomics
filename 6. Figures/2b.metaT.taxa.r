
library(readxl)
library(dplyr)

excel_sheets("metaGT.xlsx")

# metaG ===================================================
dat.metaG <- read_excel("metaGT.xlsx", sheet = "metaG")


# arrange bacteria by relative abundance
avg.dat <- 
  cbind(
    species = dat.metaG[1],
    avg = rowMeans(dat.metaG[-1])) %>%
  as.data.frame() %>%
  arrange(desc(avg))

speciesLvls <- c(avg.dat$MetaG[avg.dat$MetaG != "Others"], "Others")

#  microbGrouping  ===================================
microbGroup <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
colnames(microbGroup)[colnames(microbGroup) == "new_grouping_info_cutoff10"] <- "Grouping_new"


# arrange the samples within microbGroup -----------------
sps.Pa <- microbGroup$Sample[which(microbGroup$Grouping_new == "Pa")]
sps.Pa <- intersect(sps.Pa, colnames(dat.metaG)) 
dat.Pa <- dat.metaG %>% select(MetaG, all_of(sps.Pa))
dat.Pa.t <- dat.Pa %>% tibble::column_to_rownames("MetaG") %>% t() %>% as.data.frame()
dat.Pa.t <- dat.Pa.t %>% arrange(desc(`Pseudomonas aeruginosa`))
rownames(dat.Pa.t)


sps.Hi <- microbGroup$Sample[which(microbGroup$Grouping_new == "Hi")]
sps.Hi <- intersect(sps.Hi, colnames(dat.metaG)) 
dat.Hi <- dat.metaG %>% select(MetaG, all_of(sps.Hi))
dat.Hi.t <- dat.Hi %>% tibble::column_to_rownames("MetaG") %>% t() %>% as.data.frame()
dat.Hi.t <- dat.Hi.t %>% arrange(desc(`Haemophilus influenzae`))
rownames(dat.Hi.t)


sps.PPM <- microbGroup$Sample[which(microbGroup$Grouping_new == "PPM")]
sps.PPM <- intersect(sps.PPM, colnames(dat.metaG)) 
dat.PPM <- dat.metaG %>% select(MetaG, all_of(sps.PPM))
dat.PPM.t <- dat.PPM %>% tibble::column_to_rownames("MetaG") %>% t() %>% as.data.frame()
# the highest bacteria within PPM
colMeans(dat.PPM.t)[order(colMeans(dat.PPM.t), decreasing = T)]
dat.PPM.t <- dat.PPM.t %>% arrange(desc(`Neisseria subflava`))
rownames(dat.PPM.t)

sps.Commensal <- microbGroup$Sample[which(microbGroup$Grouping_new == "Commensal")]
sps.Commensal <- intersect(sps.Commensal, colnames(dat.metaG)) 
dat.Commensal <- dat.metaG %>% select(MetaG, all_of(sps.Commensal))
dat.Commensal.t <- dat.Commensal %>% tibble::column_to_rownames("MetaG") %>% t() %>% as.data.frame()
# the highest bacteria within Commensal 
colMeans(dat.Commensal.t)[order(colMeans(dat.Commensal.t), decreasing = T)]
dat.Commensal.t <- dat.Commensal.t %>% arrange(desc(`Neisseria subflava`))
rownames(dat.Commensal.t)


library(ggplot2)
plotDat <- dat.metaG %>% reshape2::melt(id.var="MetaG", variable.name="Sample")
plotDat$MetaG <- factor(plotDat$MetaG, levels = speciesLvls)
plotDat$Sample <- 
  factor(plotDat$Sample, 
         levels = c(rownames(dat.Pa.t),rownames(dat.Hi.t),rownames(dat.PPM.t),rownames(dat.Commensal.t)))

library(RColorBrewer)
colors <- setNames(c("#77b793","#d88c8e","#446a96","#3c9e9e","#5dbbbd",
                     "#93d2d5","#7b6a7f","#934e81","#c06445","#e37820",
                     "#ebac2c","#ebdf45","#d0b93b","#ae762d","#aa5b46",
                     "#c66b80","#d27ba6","#b6879f","#d5d5d5"),
                   c("Pseudomonas aeruginosa","Haemophilus influenzae","Neisseria subflava","Porphyromonas gingivalis","Haemophilus parainfluenzae",
                     "Prevotella melaninogenica","Prevotella intermedia","Prevotella jejuni","Escherichia coli","Neisseria flavescens",
                     "Klebsiella pneumoniae","Fusobacterium nucleatum","Aggregatibacter segnis","Veillonella atypica","Veillonella dispar",
                     "Moraxella catarrhalis","Prevotella oris","Veillonella parvula","Others"  ))



# the metaG taxonomy 
ggplot(plotDat) +
  geom_col(aes(x=Sample, y=value, fill=MetaG)) +
  scale_fill_manual(values = colors) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 90))

# MetaT ===============================================================
dat.metaT <- read_excel("metaGT.xlsx", sheet = "metaT")

library(ggplot2)

plotDat_metaT <- 
  dat.metaT %>%
  reshape2::melt(id.var="MetaT", variable.name="Sample")

# remove exacerbation in samples
fullsps.levels <- c(rownames(dat.Pa.t),rownames(dat.Hi.t),rownames(dat.PPM.t),rownames(dat.Commensal.t))
sps.exacerbation <- microbGroup$Sample[microbGroup$Group == "Exacerbations"]
sps.levels <- fullsps.levels[!fullsps.levels %in% sps.exacerbation]

plotDat_metaT <- plotDat_metaT %>% filter(!Sample %in% sps.exacerbation)
plotDat_metaT$Sample <- factor(plotDat_metaT$Sample, levels = sps.levels)

#  colors 
unique(plotDat_metaT$MetaT) %in% names(colors)

# same species arrangement as those in MetaG
plotDat_metaT$MetaT <- factor(plotDat_metaT$MetaT, levels = names(colors))


pdf("2b.metaT_taxonomy.pdf", width = 10, height = 4)
ggplot(plotDat_metaT) +
  geom_col(aes(x=Sample, y=value, fill=MetaT)) +
  scale_fill_manual(values = colors) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 90))
dev.off()
