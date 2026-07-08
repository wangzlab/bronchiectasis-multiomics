library(data.table)
library(dplyr)

dat.pval <- fread("PAD_Health_HostT_BSI.txt", data.table = F); head(dat.pval)
colnames(dat.pval)[2] <- "Feature"


dat <- fread("../_dataGZ/hostT_pathway.txt", data.table = F)
dat <- dat %>%
  tibble::column_to_rownames("#") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") 

meta_origin <- fread("../_dataGZ/group_all.txt", data.table = F) %>% mutate(Group = paste(Type, MetaG, sep="."))
head(meta_origin)
# new microb grouping 
grouping <- fread("../_data/new_grouping_cutoff10_v2.txt",data.table = F)
meta_origin$`#NAME` %in% grouping$Sample
meta_origin$MicrobGrouping <- 
  sapply(meta_origin$`#NAME`,
         function(x) grouping$new_grouping_info_cutoff10[which(grouping$Sample == x)])

dat$Sample %in% meta_origin$`#NAME`

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

dat$Sample %in% Groups_df$Sample

dat.avg <- dat %>% 
  mutate(Group = sapply(Sample, function(x) Groups_df$Group[which(Groups_df$Sample == x)])) %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  tibble::column_to_rownames("Group") %>%
  t() %>% as.data.frame()

dat.need <- merge(dat.pval %>% select(-Pa, -Health), dat.avg %>% select(Pa, Health), by.x="Feature", by.y = 0)
head(dat.need)
dat.need <- dat.need %>% arrange(Pa)
dat.need$dir <- ifelse(dat.need$Direction  == "Down","-","+")


dat.hm <- dat.need %>% select(Feature, dir, Pa, Health) %>% reshape2::melt(id.var=c("Feature","dir"), value.name = "avg")
head(dat.hm)
dat.hm$dir[dat.hm$variable == "Health"] <- NA
dat.hm$Feature <- factor(dat.hm$Feature, levels = dat.need$Feature)
dat.hm$variable <- factor(dat.hm$variable, levels = c("Health","Pa"))

library(ggplot2)

p1 <- ggplot(dat.hm) +
  geom_tile(aes(x=variable, y=Feature, fill = avg)) +
  geom_text(aes(x=variable, y=Feature, label = dir)) +
  scale_fill_gradient2(low = "#E6949A", 
                       high = "#47B1B6", 
                       mid = "white",
                       midpoint =  0, 
                       limits = c(-0.58, 0.65), # NPD和PPMD统一scale
                       na.value="gray") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "left") 
p1

dat.need$Feature <- factor(dat.need$Feature, levels = dat.need$Feature)

p2 <- ggplot(dat.need) +
  geom_tile(aes(x=1, y=Feature, fill=log10P)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(palette="Greys",
                       limits=c(2.4, 14.6), # NPD和PPMD统一scale
                       direction = 1 )
p2


# plot of adonis 
p3 <- ggplot(dat.need) +
  geom_tile(aes(x=1,y=Feature, fill=Adonis))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(limits=c(0,0.22),  # NPD和PPMD统一scale
                       palette="Blues",
                       direction = 1 )

# plot of estimate
p4 <- ggplot(dat.need) +
  geom_tile(aes(x=1, y=Feature, fill=Estimate))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_gradient2(low = "#E6949A", 
                       high = "#47B1B6", 
                       mid = "white",
                       midpoint =  0, 
                       limits = c(-8.3,5.0), # NPD和PPMD统一scale
                       na.value="gray")


library(ggpubr)

ggarrange(p1,p2,p3,p4, widths = c(0.4,0.2,0.2,0.2), nrow = 1)
ggsave(filename = "S4c.Pa.health_hostTpathway.hm.pdf", width = 8, height = 4.5)

