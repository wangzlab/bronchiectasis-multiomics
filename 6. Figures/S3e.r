list.files()

library(data.table)

dat <- fread("S3e.txt", data.table = F); head(dat)
dat.l <- dat %>%
   reshape2::melt(id.var="V1", variable.name="Cluster")
dat.l$panel = sapply(strsplit(dat.l$V1, "_", fixed = T),"[[", 1)
dat.l$V1 <- factor(dat.l$V1, 
                   levels = c("HID_up_stable","HID_up_exac","NPD_up_stable" ,"NPD_up_exac","PAD_up_stable","PAD_up_exac","PPMD_up_stable","PPMD_up_exac"  ))


library(ggplot2)

library(ggalluvial)


ggplot(dat.l,
         aes(x = V1,       
             y = value,            
             stratum = Cluster,      
             alluvium = Cluster,    
             fill = Cluster,         
             label = Cluster)) +   
  geom_alluvium(                
    aes(fill = Cluster),
    alpha = 0.7,                  
    width = 0.4,                 
    color = NA
  ) +
  geom_stratum(width = 0.4, color=NA)+
  facet_wrap(vars(panel),ncol=2,nrow=2, scale="free")  +   
  geom_text(stat = "stratum",    
            size = 3.5) +
  scale_fill_manual(values = c("Pa" = "#91c4a7", 
                               "Hi" = "#e1a3a3", 
                               "PPM" = "#a0c8d8",
                               "Commensal" = "#eac78b")) +
  theme_minimal() 
ggsave(filename = "S3e.barplot.with.connection.pdf", width = 5, height = 4)
