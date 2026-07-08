
library(data.table)
library(dplyr)


dat <- fread("metaGT_FC.txt", data.table = F)
dat.l <- dat %>%
  reshape2::melt(id.var="V1")
head(dat.l)

dat.l$V1 <- factor(dat.l$V1, levels = rev(dat$V1))

library(ggplot2)

ggplot(dat.l) +
  geom_tile(aes(x=variable, y=V1, fill = value), color="#E3E3E3") +
  scale_fill_gradient2(high = "#E6949A", 
                       low = "#47B1B6", 
                       mid = "white",
                       midpoint =  0, 
                       # limits = c(-1.28, 1.28), # NPD和PPMD统一scale
                       na.value="gray") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(
  filename = "S7a.species_AM_heatmap.pdf", width = 4, height = 5.5
)
