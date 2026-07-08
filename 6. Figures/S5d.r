list.files()
library(data.table)
library(dplyr)

dat <- fread("Pm.txt", data.table = F)
head(dat)

library(ggplot2)

ggplot() +
  geom_point(data = dat, aes(x=rsq, y=rmse), color="gray")+
  geom_point(data = dat %>% filter(!is.na(Highlight)),
             aes(x=rsq, y = rmse), color="#C30D23", size=2) +
  theme_bw()+
  theme(panel.grid = element_blank())
