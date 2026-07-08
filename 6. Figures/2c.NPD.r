
library(data.table)
library(dplyr)

dat <- fread("NPD_Health_taxonomy.txt",data.table = F) 
dat$metaG_Pval <- -log10(dat$metaG_Pval)
dat$metaT_Pval <- -log10(dat$metaT_Pval)

dat <- dat %>% arrange(metaG_FC)
dat$species <- factor(dat$species, levels = c(dat$species))

dat.bb <- dat %>%
  select(species, metaG_FC, metaT_FC) %>%
  reshape2::melt(id.var="species", value.name = "FC")

dat.bb$species <- factor(dat.bb$species, levels = c(dat$species))

cp1.col <- "#E6949A"; cp2.col <- "#47B1B6"
min_val <- min(dat.bb$FC)
max_val <- max(dat.bb$FC)
break1 <- -4.5   
break2 <- 0    
break3 <- 4.5    

rescale01 <- function(x) (x - min_val) / (max_val - min_val)

pos1 <- rescale01(break1)  # -4 -> 4/13 ≈ 0.3077
pos2 <- rescale01(break2)  # 0  -> 8/13 ≈ 0.6154
pos3 <- rescale01(break3)  # 4  -> 12/13 ≈ 0.9231

colors <- c("#E6949A", "#E6949A", "white", "#47B1B6", "#47B1B6")
values <- c(0, pos1, pos2, pos3, 1)

p1 <- ggplot(dat) +
  geom_point(aes(x=metaG_FC, y=species, size=abs(metaG_FC), fill=metaG_FC), shape=21) +
  geom_point(aes(x=metaT_FC, y=species, size=abs(metaT_FC), fill=metaT_FC), shape=23) +
  theme_bw() +
  theme(legend.position = "left") +
  scale_size(range = c(2,6)) +
  # ggtitle("NPD_Health") +
  xlab("metaG/metaT FC" )+
  scale_fill_gradientn(
    colors = colors,
    values = values,
    limits = c(min_val, max_val),
    oob = scales::squish,          # 将超出 limits 的值钳制到边界颜色
    name = "Value"
  ) 
p1

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
  scale_fill_distiller(palette="Greys", direction = 1, limits=c(0,7.48))
p2



library(ggpubr)
ggarrange(p1, p2, widths = c(0.8,0.2))
ggsave(filename = "2c.NPD_health.pdf", width = 7, height = 4)
