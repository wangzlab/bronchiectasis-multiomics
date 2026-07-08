library(data.table)
library(dplyr)
dat<- fread("CI_Shannon_Stable_only.txt",data.table = F)

head(dat)

colors <- setNames(c("#ECBF71","#E29192","#8DC0A3","#8CC2D8"),
                   nm = c("Commensal","Hi","Pa","PPM" ))

library(ggplot2)
ggplot(dat) +
  geom_point(aes(y=BDI, x=Shannon, color=Group)) +
  geom_smooth(aes(y=BDI, x=Shannon), method = "lm") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = "S7d.pdf", width = 5, height = 3.5)
cor.test(x=dat$Shannon, y=dat$BDI, method = "spearman")
