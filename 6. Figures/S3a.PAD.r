library(data.table)
library(dplyr)

dat.pval <- fread("PAD_Healthy_Module_BSI.txt", data.table = F) 
colnames(dat.pval)[1] <- "Feature"

dat <- fread("../_dataGZ/metat_module_sum.txt", data.table = F)
IDmatch <- fread("../_dataGZ/metat_match.txt", data.table = F, header = F)
head(IDmatch)

meta_origin <- fread("../_dataGZ/group_all.txt", data.table = F) %>% mutate(Group = paste(Type, MetaG, sep="."))
head(meta_origin)

grouping <- fread("../_data/new_grouping_cutoff10_v2.txt",data.table = F)
meta_origin$`#NAME` %in% grouping$Sample
meta_origin$MicrobGrouping <- 
  sapply(meta_origin$`#NAME`, function(x) grouping$new_grouping_info_cutoff10[which(grouping$Sample == x)])

all(IDmatch$V2 %in% meta_origin$`#NAME`)
IDmatch$Group <- sapply(IDmatch$V2, function(x) meta_origin$MicrobGrouping[which(meta_origin$`#NAME` == x)])  #此处采用microbGrouping: Pa, Hi, PPM, Commensal
colnames(dat)[!colnames(dat) %in% IDmatch$V1] 

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
  tibble::column_to_rownames("#SampleID") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("V1") %>%
  mutate(Sample = sapply(V1, function(x)IDmatch$V2[which(IDmatch$V1 == x)] )) %>%
  dplyr::select(-V1)

tmp <- log1p(as.matrix(dat %>% tibble::column_to_rownames("Sample")) )
dim(tmp)

tmp  <- tmp[c(sps.Pa, sps.Health),]
if(any(is.na(tmp))) print("Samples not matched")


rm.modules <- which(apply(tmp,2,function(x)all(x==0)))
tmp <- tmp[,-rm.modules]

dat.st <- apply(log1p(tmp), 2, function(x) (x-mean(x))/sd(x)) 
dat.st <- as.data.frame(dat.st) %>% tibble::rownames_to_column("Sample")
sapply(dat.st, mean)
sapply(dat.st, sd)

dat.avg <- dat.st %>% 
  mutate(Group = sapply(Sample, function(x)Groups_df$Group[which(Groups_df$Sample == x)])) %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  tibble::column_to_rownames("Group") %>%
  t() %>% as.data.frame()

dat.need <- merge(dat.pval, dat.avg %>% dplyr::select(Pa, Health), by.x="Feature", by.y = 0)
dat.need <- dat.need %>% arrange(Pa)
dat.need$dir <- ifelse(dat.need$Direction == "Down","-","+")


dat.hm <- dat.need %>% dplyr::select(Feature, dir, Pa, Health) %>% reshape2::melt(id.var=c("Feature","dir"), value.name = "avg")
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
                       # limits = c(-1.02,1.09), 
                       na.value="gray") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "left") 


dat.need$Feature <- factor(dat.need$Feature, levels = dat.need$Feature)
p2 <- ggplot(dat.need) +
  geom_tile(aes(x=1, y=Feature, fill=log10P)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(limits=c(1.6,7), 
                       palette="Greys",
                       direction = 1 )


# plot of contribution
p3 <- ggplot(dat.need) +
  geom_tile(aes(x=1,y=Feature, fill=Contribution))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(limits=c(0,1),  
                       palette="Blues",
                       direction = 1 )

# plot of estimate
p4 <- ggplot(dat.need) +
  geom_tile(aes(x=1, y=Feature, fill=Estimate))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_gradientn(
    colours = c("#E6949A", "white", "#47B1B6"),
    values = scales::rescale(
      c(-0.37, 0, 0.29),  
      from = c(-0.37, 0.29)  
    ),
    limits = c(-0.37, 0.29),
    na.value = "gray",
    oob = scales::squish
  )

library(ggpubr)
ggarrange(p1,p2,p3,p4, widths = c(0.4,0.2,0.2,0.2), nrow = 1)
ggsave(filename = "S3a.PAD_Health_module.pdf", width = 8, height = 4.5)

