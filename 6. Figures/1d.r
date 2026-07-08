library(data.table)
library(dplyr)
library(ggplot2)


# microbe Clusters
microbGrouping <- fread("../_data/new_grouping_cutoff10_v2.txt", data.table = F)
colnames(microbGrouping)[2] <- "cluster"

cluster.colors <- setNames(c("#7abf98","#8cc2d8","#ecbf71","#e29192"),
                           nm=c("Pa","PPM","Commensal","Hi"))


cluster.cohort1 <- 
  microbGrouping %>%
  group_by(Cohort, cluster) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))
cluster.cohort1$Cohort <- factor(cluster.cohort1$Cohort, levels = c("China","CAMEB2","EMBARC"))
cluster.cohort1$cluster <- factor(cluster.cohort1$cluster, levels = c("Pa","Hi","PPM","Commensal"))
ggplot(cluster.cohort1) +
  geom_col(aes(x=Cohort,y=perc,fill = cluster)) +
  scale_fill_manual(values = cluster.colors) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = "1d.Cohort_MicrobCluster_bar.pdf", width = 4, height = 4)


# regions
meta.617 <- fread("../_data/combinedData_617sps.csv", data.table = F) 
microbGrouping$Region <- 
  sapply(1:nrow(microbGrouping),
         function(i){
           if(microbGrouping$Cohort[i] == "China"){
             ""
           }else{
             meta.617$European_Region[which(meta.617$Sequencing.ID == microbGrouping$Sample[i])]
           }
         })
table(microbGrouping$Region )


# EMBARC
cluster.EMBARC <- microbGrouping %>%
  filter(Cohort == "EMBARC") %>%
  group_by(Region, cluster) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))

cluster.EMBARC$cluster <- factor(cluster.EMBARC$cluster, levels = c("Pa","Hi","PPM","Commensal"))

ggplot(cluster.EMBARC) +
  geom_col(aes(x=Region,y=perc,fill = cluster)) +
  scale_fill_manual(values = cluster.colors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))
ggsave(filename = "1d.EMBARCRegion_MicrobCluster_bar.pdf", width = 4, height = 5)


# CAMEB2
cluster.CAMEB2 <- microbGrouping %>%
  filter(Cohort == "CAMEB2") %>%
  group_by(Region, cluster) %>%
  summarise(n=n()) %>%
  mutate(perc = n/sum(n))

cluster.CAMEB2$cluster <- factor(cluster.CAMEB2$cluster, levels = c("Pa","Hi","PPM","Commensal"))

ggplot(cluster.CAMEB2) +
  geom_col(aes(x=Region,y=perc,fill = cluster)) +
  scale_fill_manual(values = cluster.colors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))
ggsave(filename = "1d.CAMEB2Region_MicrobCluster_bar.pdf", width = 4, height = 4)
