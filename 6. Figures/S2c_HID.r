library(data.table)
library(dplyr)

dat <- fread("HID_Health_taxonomy.txt", data.table = F)
head(dat)

dat <- dat %>% arrange(metaG_FC)

dat$species <- factor(dat$species, levels = c(dat$species))

dat.bb <- dat %>%
  select(species, metaG_FC, metaT_FC) %>%
  reshape2::melt(id.var="species", value.name = "FC")

dat.bb$species <- factor(dat.bb$species, levels = c(dat$species))


cp1.col <- "#E6949A"; cp2.col <- "#47B1B6"

values.full <- c(dat$metaG_FC, dat$metaT_FC)

p1 <- ggplot(dat) +
  geom_point(aes(x=metaG_FC, y=species, size=abs(metaG_FC), fill=metaG_FC), shape=21) +
  geom_point(aes(x=metaT_FC, y=species, size=abs(metaT_FC), fill=metaT_FC), shape=23) +
  theme_bw() +
  theme(legend.position = "left") +
  scale_size(range = c(2,6)) +
  # ggtitle("NPD_Health") +
  xlab("metaG/metaT FC" )+
  #scale_fill_gradient2(low = cp1.col, # low for cp[1]
  #                     high = cp2.col, # high for cp[2]
  #                     mid = "white",
  #                     midpoint = 0,
  #                     #limits = c(-1.53,2.45),
  #                     na.value="white")  #数值0作为中点
  scale_fill_gradientn(
    colors = c("#E6949A", "#E9E9E9","#47B1B6" ),  # 自定义颜色
    values = scales::rescale(c(min(min(values.full), 0)-0.1,
                               0,  
                               max(0,max(values.full))+0.1 )),  
    limits = c(min(min(values.full), 0)-0.1,   max(0,max(values.full))+0.1 ),
    na.value = "white"  # 设置 NA 值的颜色
  )

dat.hm <-
  dat %>%
  select(species, metaG_Pval, metaT_Pval) %>%
  reshape2::melt(id.var="species", value.name = "pval")
dat.hm$species <- factor(dat.hm$species, levels = dat$species)

p2 <- ggplot(dat.hm) +
  geom_tile(aes(x=variable, y=species, fill=pval)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_distiller(palette="Greys", 
                       limits=c(0,7.48),
                       oob = scales::squish, 
                       direction = 1)


library(ggpubr)
ggarrange(p1, p2, widths = c(0.8,0.2))
ggsave(filename = "S2c.HID_Health.pdf", width = 7, height = 4)

