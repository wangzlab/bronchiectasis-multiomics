library(data.table)
library(dplyr)


# meta--------------
vars <- c("BSI","AE_1year","FEV1pred_perc","MMRC")
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
meta.118.stable$BSI <- as.integer(meta.118.stable$BSI)
sapply(meta.118.stable, class)
meta.118.stable$AE_1year <- as.integer(meta.118.stable$AE_1year)
meta.118.stable$MMRC <- as.integer(meta.118.stable$MMRC)


meta <- bind_rows(
  meta.118.stable %>% mutate(Study = "China"),
  meta.617
)

# groups ------------------------------
ClusterID <- data.table::fread("../_data/new_grouping_cutoff10_v2.txt", data.table=F)
head(ClusterID)
colnames(ClusterID)[1:2] <- c("Sample", "Clusters")
ClusterID$Clusters %>% unique



# merge source data to plot  ------------------------------
sourceDat <- merge(meta, ClusterID, by.x = "SampleID",by.y = "Sample")
sourceDat$BSI.cat <-  cut(sourceDat$BSI,
                          breaks = c(-Inf,5,9,Inf),
                          labels = c("mild","moderate","severe"),
                          include.lowest = T, right = F)


sourceDat$AE_1year.cat <- cut(sourceDat$AE_1year,
                              breaks = c(0,1,2,3,Inf),
                              labels = c("0","1","2",">2"),
                              include.lowest = T, right = F)




#  BSI : -----------
pD.cat <- sourceDat %>%
  dplyr::select(SampleID, Study, Clusters, BSI.cat) %>%
  reshape2::melt(id.vars=c("SampleID","Study","Clusters"), variable.name="Variable") %>% 
  filter(!is.na(value)) %>%
  group_by(Study, Clusters, Variable, value) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))
head(pD.cat)
pD.cat$value <- factor(pD.cat$value, levels = c("mild","moderate","severe"))
pD.cat$Study <- factor(pD.cat$Study, levels = c("China","CAMEB2","EMBARC"))
pD.cat$Clusters <- factor(pD.cat$Clusters, levels = c("Pa", "Hi", "PPM", "Commensal"))


library(ggplot2)
cat.colors <- setNames(c("#f7c1c3","#c25f77","#64213d"),
                       nm=c("mild","moderate","severe"))

ggplot(pD.cat) +
  geom_col(aes(x=Clusters, y=perc, fill = value)) +
  facet_grid(Variable~Study) +
  scale_fill_manual(values = cat.colors) +
  theme_bw()
ggsave(filename = "1e.BSI_separateCohorts.pdf", width = 9, height = 3)

#  full data 
pD2.cat <- sourceDat %>%
  dplyr::select(SampleID, Study, Clusters, BSI.cat) %>%
  reshape2::melt(id.vars=c("SampleID","Study","Clusters"), variable.name="Variable") %>% 
  filter(!is.na(value)) %>%
  group_by( Clusters, Variable, value) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))
head(pD2.cat)
pD2.cat$value <- factor(pD2.cat$value, levels = c("mild","moderate","severe"))
pD2.cat$Clusters <- factor(pD2.cat$Clusters, levels = c("Pa", "Hi", "PPM", "Commensal"))
ggplot(pD2.cat) +
  geom_col(aes(x=Clusters, y=perc, fill = value)) +
  scale_fill_manual(values = cat.colors) +
  theme_bw()
ggsave(filename = "1e.BSI_fullData.pdf", width = 4.5, height = 3)


