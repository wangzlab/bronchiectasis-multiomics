library(data.table)
library(dplyr)

dat.pval <- fread("PPMD_Health_MetaB_BSI.txt", data.table = F)

#feature match list
metabMatch <- fread("metab_matchlist.txt", data.table = F, header = F)
dat.pval$V1 %in% metabMatch$V2
dat.pval$Feature.full <- sapply(dat.pval$V1, function(x) metabMatch$V1[which(metabMatch$V2 == x)])
dat.pval$Feature.full <- paste0("ME", dat.pval$Feature.full)

dat <- fread("../_dataGZ/metaB_module.txt", data.table = F)
meta_origin <- fread("../_dataGZ/group_all.txt", data.table = F) %>% mutate(Group = paste(Type, MetaG, sep="."))
head(meta_origin)
# new microb grouping 
grouping <- fread("../_data/new_grouping_cutoff10_v2.txt",data.table = F)
meta_origin$`#NAME` %in% grouping$Sample
meta_origin$MicrobGrouping <- 
  sapply(meta_origin$`#NAME`,
         function(x) grouping$new_grouping_info_cutoff10[which(grouping$Sample == x)])

sps.PPM <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "PPM")]
sps.Pa <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "Pa")]
sps.Hi <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "Hi")]
sps.Commensal <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable" & meta_origin$MicrobGrouping == "Commensal")]

sps.Health <- meta_origin$`#NAME`[which(meta_origin$Type == "Health")]
sps.Stable <- meta_origin$`#NAME`[which(meta_origin$Type == "Stable")]
sps.Exacerbation <- meta_origin$`#NAME`[which(meta_origin$Type == "Exacerbation")]



Groups_df <- 
  bind_rows(
    cbind.data.frame(Sample=sps.Health, Group="Health", stringsAsFactors=F),
    cbind.data.frame(Sample=sps.Exacerbation, Group="Exacerbation", stringsAsFactors=F ),
    #cbind.data.frame(Sample=sps.Stable, Group="Stable", stringsAsFactors=F ),
    cbind.data.frame(Sample=sps.PPM, Group="PPM", stringsAsFactors=F ),
    cbind.data.frame(Sample=sps.Pa, Group="Pa", stringsAsFactors=F ),
    cbind.data.frame(Sample=sps.Hi, Group="Hi", stringsAsFactors=F ),
    cbind.data.frame(Sample=sps.Commensal, Group="Commensal", stringsAsFactors=F ))


dat <- dat %>%
  tibble::column_to_rownames("V1") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") 
dat$Sample %in% Groups_df$Sample


# zscore within the two groups
tmp <- dat %>% tibble::column_to_rownames("Sample") %>% as.matrix()

sps.keep <- c(sps.PPM, sps.Health)[c(sps.PPM, sps.Health) %in% rownames(tmp)] 
tmp  <- tmp[sps.keep,]

rm.modules <- which(apply(tmp,2,function(x)all(x==0)))
if(length(rm.modules) > 0) tmp <- tmp[,-rm.modules]

dat.st <- apply(tmp, 2, function(x) (x-mean(x))/sd(x)) 
dat.st <- as.data.frame(dat.st) %>% tibble::rownames_to_column("Sample")
sapply(dat.st, mean)
sapply(dat.st, sd)

dat.avg <- dat.st %>% 
  mutate(Group = sapply(Sample, function(x) Groups_df$Group[which(Groups_df$Sample == x)])) %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  tibble::column_to_rownames("Group") %>%
  t() %>% as.data.frame()

dat.need <- merge(dat.pval %>% select(-Health, -PPM), dat.avg , by.x="Feature.full", by.y = 0)
head(dat.need)
dat.need <- dat.need %>% arrange(PPM)
dat.need$dir <- ifelse(dat.need$Direction == "Down", "-", "+")


dat.hm <- dat.need %>% select(Feature.full,dir, PPM, Health) %>% reshape2::melt(id.var=c("Feature.full","dir"), value.name = "avg")
head(dat.hm)
dat.hm$dir[dat.hm$variable == "Health"] <- NA
dat.hm$Feature.full <- factor(dat.hm$Feature.full, levels = dat.need$Feature.full)
dat.hm$variable <- factor(dat.hm$variable, levels = c("Health","PPM"))

library(ggplot2)

p1 <- ggplot(dat.hm) +
  geom_tile(aes(x=variable, y=Feature.full, fill = avg)) +
  geom_text(aes(x=variable, y=Feature.full, label = dir)) +
  scale_fill_gradient2(low = "#E6949A", 
                       high = "#47B1B6", 
                       mid = "white",
                       midpoint =  0, 
                       limits = c(-1.28, 1.26), # unify scale for NPD and PPMD
                       na.value="gray") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "left") 
p1


dat.need$Feature.full <- factor(dat.need$Feature.full, levels = dat.need$Feature.full)

p2 <- ggplot(dat.need) +
  geom_tile(aes(x=1, y=Feature.full, fill=log10P)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(palette="Greys",
                       limits=c(3,6.5),  #  unify scale for NPD and PPMD
                       direction = 1 )
p2



# plot of estimate
p3 <- ggplot(dat.need) +
  geom_tile(aes(x=1, y=Feature.full, fill=Estimate))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient2(low = "#E6949A", 
                       high = "#47B1B6", 
                       mid = "white",
                       midpoint =  0, 
                       limits = c(-14.1,13.6), #  unify scale for NPD and PPMD
                       na.value="gray")



# module  composition =============================================

md.abb.pthy <- fread("../_dataGZ/metaB_module_description.txt", data.table = F, col.names = c("module","pathway"))
md.classes <- fread("../_dataGZ/metaB_module_class.txt",data.table = F)

md.compositions <- md.classes %>%
  filter(Class != "Unclassified") %>%
  group_by(Module, Class) %>%
  summarise(n=n()) %>%
  mutate(rel = n/sum(n)) %>%
  as.data.frame() %>%
  mutate(Module = paste0("ME", Module))

colors.mdClass <- c("#b5c3d5","#df9797","#9abca7","#cbe7ec","#f6eacc",
                    "#d5d1e5","#d1ebed","#eeced9",
                    "#f7e672", "#a1d6d4", "#e0d0bc","#d6d6d6")
names(colors.mdClass) <- c("Alkaloids","Amino acids, peptides, and analogues",   "Benzenoids","Carbohydrates","Indole and derivatives",
                           "Lipids and lipid-like molecules","Nucleosides, nucleotides, and analogues","Organic acids",  
                           "Organoheterocyclic compounds","Phenylpropanoids and polyketides","Terpenoids","Others"  )


pD.md <- md.compositions %>%  filter(Module %in%  dat.need$Feature.full) 
pD.md$Module <- factor(pD.md$Module, levels = dat.need$Feature.full)


class.abc <- unique(pD.md$Class)[order(unique(pD.md$Class))]
class.lvls <- c(class.abc[class.abc != "Others"], "Others")
pD.md$Class <- factor(pD.md$Class, levels = rev(class.lvls))
p4 <- ggplot(pD.md) +
  geom_col(aes(x=rel,y=Module, fill=Class)) +
  scale_fill_manual(values = colors.mdClass) +
  theme_bw() + theme(panel.grid = element_blank())
p4


# plot of adonis
p5 <- ggplot(dat.need) +
  geom_tile(aes(x=1,y=Feature.full, fill=adonis))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(
    limits=c(0,0.29),  # unify scale for NPD and PPMD
    palette="Blues",
    direction = 1 )

library(ggpubr)
ggarrange(p1,p2,p3,p5,p4,  ncol = 5, widths = c(0.25,0.1, 0.11,0.1,0.44), nrow = 1)
ggsave(filename = "3a_PPMD.health_FC.hm_moduleComposition.pdf", width = 12, height = 4.5)
