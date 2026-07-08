library(data.table)
library(dplyr)
library(vegan)

# omic data =====================================
metaT.tx <- fread("../_dataGZ/metaT_taxa.txt", data.table=F)
metaT.ko <- fread("../_dataGZ/metaT_ko.txt", data.table=F)
metaG.tx <- fread("../_dataGZ/metag_taxa.txt",data.table = F)
metaG.ko <- fread("../_dataGZ/metag_ko.txt",data.table = F)

colnames(metaT.tx)
colnames(metaT.ko)
colnames(metaG.tx)
colnames(metaG.ko) 

colnames(metaT.ko) %in% colnames(metaT.tx)
colnames(metaG.tx) %in% colnames(metaT.tx)
colnames(metaG.ko) %in% colnames(metaT.tx)

matches <- fread("../_dataGZ/match.txt", data.table = F, header = F)
colnames(metaG.ko)[-1] <- sapply(colnames(metaG.ko)[-1], function(x) matches$V2[which(matches$V1 == x)])


# bray-curtis distance ------------------------
BC.metaT.tx <- vegdist(t(metaT.tx[,-1]), method="bray") %>% as.matrix()
BC.metaT.tx[1:3,1:3]
dim(BC.metaT.tx)

BC.metaT.tx <- BC.metaT.tx %>% reshape2::melt(value.name="BC.metaT.tx")
BC.metaT.ko <- vegdist(t(metaT.ko[,-1]), method="bray") %>% as.matrix() %>% reshape2::melt(value.name="BC.metaT.ko")
BC.metaG.tx <- vegdist(t(metaG.tx[,-1]), method="bray") %>% as.matrix() %>% reshape2::melt(value.name="BC.metaG.tx")
BC.metaG.ko <- vegdist(t(metaG.ko[,-1]), method="bray") %>% as.matrix() %>% reshape2::melt(value.name="BC.metaG.ko")

# paired sample IDs ----------------------------
meta_origin <- fread("../_dataGZ/meta_confounder.txt", data.table = F)
pairedIDs <- meta_origin %>%
  filter(Group %in% c("Exacerbations","Stable")) %>%
  reshape2::dcast(SubjectID~Group, value.var = "SampleID")
pairedIDs <- pairedIDs[complete.cases(pairedIDs),]

distances  <- merge(merge(merge(merge(pairedIDs, BC.metaT.tx, by.x = c("Exacerbations", "Stable"), by.y = c("Var1","Var2")),
                                BC.metaT.ko, by.x = c("Exacerbations", "Stable"), by.y = c("Var1","Var2")),
                          BC.metaG.tx, by.x = c("Exacerbations", "Stable"), by.y = c("Var1","Var2")),
                    BC.metaG.ko, by.x = c("Exacerbations", "Stable"), by.y = c("Var1","Var2"))


dist.l <- distances %>%
  reshape2::melt(id.vars=c("Exacerbations",  "Stable", "SubjectID"), value.var="BCdistances")
library(stringr)
dist.l$molecule = str_extract(dist.l$variable, "(?<=\\.)[^.]+(?=\\.)")
dist.l$type = sub("^[^.]*\\.[^.]*\\.", "", dist.l$variable)


library(ggpubr)
ggboxplot(dist.l, x="molecule", y="value", outlier.shape = NA) +
  geom_jitter(aes(x=molecule, y=value)) +
  facet_grid(.~type) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means()


ggpaired(dist.l, x = "molecule", y = "value",
         color = "molecule", line.color = "gray", line.size = 0.4,
        # palette = "jco",
         facet.by = "type")+
  stat_compare_means(paired = TRUE) +
  ylab("Bray-curtis distance")
ggsave(filename = "S3d.metaGT.BCdist.pdf", width = 4, height = 3)
