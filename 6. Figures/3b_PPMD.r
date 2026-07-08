library(data.table)
library(dplyr)

dat.pval <- fread("PPMD_Health_HostT_BSI.txt", data.table = F); head(dat.pval)
colnames(dat.pval)[2] <- "Feature.full"

dat <- fread("../_dataGZ/hostT_pathway.txt", data.table = F)
meta_origin <- fread("../_dataGZ/group_all.txt", data.table = F) %>% mutate(Group = paste(Type, MetaG, sep="."))
head(meta_origin)
# new microb grouping 
grouping <- fread("../_data/new_grouping_cutoff10_v2.txt",data.table = F)
meta_origin$`#NAME` %in% grouping$Sample
meta_origin$MicrobGrouping <- 
  sapply(meta_origin$`#NAME`,
         function(x) grouping$new_grouping_info_cutoff10[which(grouping$Sample == x)])

colnames(dat) %in% meta_origin$`#NAME`

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

# dat转换成colnames是Module
dat <- dat %>%
  tibble::column_to_rownames("#") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") 
dat$Sample %in% Groups_df$Sample

dat.avg <- dat %>% 
  mutate(Group = sapply(Sample, function(x) Groups_df$Group[which(Groups_df$Sample == x)])) %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  tibble::column_to_rownames("Group") %>%
  t() %>% as.data.frame()


dat.need <- merge(dat.pval %>% select(-PPM, -Health), dat.avg %>% select(PPM, Health), by.x="Feature.full", by.y = 0)
head(dat.need)
dat.need <- dat.need %>% arrange(PPM)
dat.need$dir <- ifelse(dat.need$Direction == "Down","-","+")

dat.hm <- dat.need %>% select(Feature.full, dir, PPM, Health) %>% reshape2::melt(id.var=c("Feature.full","dir"), value.name = "avg")
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
                       limits = c(-0.58, 0.65), # unify scale for NPD and PPMD
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
                       limits=c(2.4, 14.6), # unify scale for NPD and PPMD
                       direction = 1 )
p2



# plot of adonis 
p3 <- ggplot(dat.need) +
  geom_tile(aes(x=1,y=Feature.full, fill=Adonis))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(
    limits=c(0,0.22),  # unify scale for NPD and PPMD
    palette="Blues",
    direction = 1 )

# plot of estimate
p4 <- ggplot(dat.need) +
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
                       limits = c(-8.3,5.0), # unify scale for NPD and PPMD
                       na.value="gray")




library(ggpubr)

ggarrange(p1,p2,p3,p4, widths = c(0.4,0.2,0.2,0.2), nrow = 1)
ggsave(filename = "3b_PPMD.health_hostTpathway.hm.pdf", width = 8, height = 4.5)
