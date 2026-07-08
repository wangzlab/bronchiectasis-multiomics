library(data.table)
library(dplyr)

# meta --------------
vars <- c("AE_1year")
covars <- c("Age","Gender")


meta.617 <- fread("../_data/combinedData_617sps.csv", data.table = F) 
meta.617 <- meta.617 %>% 
  rename(AE_1year = Exacerbation_frequency, FEV1pred_perc = FEV1_percent_predicted, SampleID=Sequencing.ID) %>% 
  dplyr::select(SampleID, Study, all_of(vars), all_of(covars))
sapply(meta.617, class)


meta.118 <- readxl::read_excel("../_data/metadata_118sps_more.xlsx")  
meta.118$Group %>% table()
meta.118.stable <- meta.118 %>% filter(Group == "Stable")  %>% dplyr::select(SampleID, all_of(vars), all_of(covars))
sapply(meta.118.stable, class)
sapply(meta.118.stable, class)
meta.118.stable$AE_1year <- as.integer(meta.118.stable$AE_1year)


meta <- bind_rows(
  meta.118.stable %>% mutate(Study = "China"),
  meta.617
)

# groups ------------------------------
ClusterID <- data.table::fread("../_data/new_grouping_cutoff10_v2.txt", data.table=F)
colnames(ClusterID) <- c("Sample", "Clusters","Cohort","Disease")
ClusterID$Clusters %>% unique

# 合并source data供计算画图  ------------------------------
sourceDat <- merge(meta, ClusterID, by.x = "SampleID",by.y = "Sample")
sourceDat$AE_1year.cat <- cut(sourceDat$AE_1year,
                              breaks = c(0,1,2,3,Inf),
                              labels = c("0","1","2",">2"),
                              include.lowest = T, right = F)


#  categorical variables: -----------
pD.cat <- sourceDat %>%
  dplyr::select(SampleID, Study, Clusters, AE_1year.cat) %>%
  reshape2::melt(id.vars=c("SampleID","Study","Clusters"), variable.name="Variable") %>% 
  filter(!is.na(value)) %>%
  group_by(Study, Clusters, Variable, value) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))

head(pD.cat)
table(pD.cat$value)
pD.cat$value <- factor(pD.cat$value, levels =  c("0","1","2",">2"))
pD.cat$Study <- factor(pD.cat$Study, levels = c("China","CAMEB2","EMBARC"))
pD.cat$Clusters <- factor(pD.cat$Clusters, levels = c("Pa", "Hi", "PPM", "Commensal"))


library(ggplot2)
cat.colors <- setNames(c("#CFE3A0","#99CE9C","#006b6d","#09404f"),
                       nm=c("0","1","2",">2"))



ggplot(pD.cat) +
  geom_col(aes(x=Clusters, y=perc, fill = value)) +
  facet_grid(Variable~Study) +
  scale_fill_manual(values = cat.colors) +
  theme_bw()
ggsave(filename = "1e.AE_seperateCohorts.pdf", width = 9, height = 3)



#  不分cohort 
pD2.cat <- sourceDat %>%
  dplyr::select(SampleID, Study, Clusters, AE_1year.cat) %>%
  reshape2::melt(id.vars=c("SampleID","Study","Clusters"), variable.name="Variable") %>% 
  filter(!is.na(value)) %>%
  group_by( Clusters, Variable, value) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))
head(pD2.cat)
pD2.cat$value <- factor(pD2.cat$value, levels = c("0","1","2",">2"))
pD2.cat$Clusters <- factor(pD2.cat$Clusters, levels = c("Pa", "Hi", "PPM", "Commensal"))
ggplot(pD2.cat) +
  geom_col(aes(x=Clusters, y=perc, fill = value)) +
#  facet_grid(Variable~.) +
  scale_fill_manual(values = cat.colors) +
  theme_bw()
ggsave(filename = "1e.AE_full.pdf", width = 4.5, height = 3)
