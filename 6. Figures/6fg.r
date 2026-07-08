library(data.table)
library(ggplot2)
library(dplyr)

# 6f =====================================================
dat <- fread("LT_master_data.txt", data.table = F)


ggplot(dat) +
  geom_line(aes(x=TypeSamples, y=CI, group=SputumSampleNo, color=TmToNxtEx)) +
  geom_point(aes(x=TypeSamples, y=CI, color=TmToNxtEx), color="gray", size=2) +
  theme_bw()+
  ylab("BDI")+
  scale_color_manual(values = c("indianred","lightgray"))
ggsave(filename = '6f.CI_BSI.EX1.P1.pdf', width = 5.5, height = 4)


# 6f --scatter_WksToNxtEx.vs.CIdiff  --------------------------------
dat <- fread("CI_diff_exac_assoc(1).txt", data.table = F)
head(dat)

covars2Select <- c("FEV1","BSI","Antibiotic","Antibiotic_class")

Covars <- c(as.list(covars2Select),
            combn(covars2Select, m = 2, simplify = F),
            combn(covars2Select, m = 3, simplify = F))
Covars[c(10,13,14)] <- NULL


dat$TmToNxtEx %>% unique
dat$TmToNxtEx <- factor(dat$TmToNxtEx , levels = c("LessThan12w" , "MoreThan12w"))
#dat$Antibiotic <- as.factor(dat$Antibiotic)
#dat$Antibiotic_class <- as.factor(dat$Antibiotic_class)

ggplot(dat) +
  geom_point(aes(x=WksToNxtEx, y=`CI diff`, color=TmToNxtEx), size=2) +
  geom_smooth(aes(x=WksToNxtEx, y=`CI diff`), method = "lm", linetype="dashed", color="black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("indianred","gray"))
ggsave("6f.scatter_WksToNxtEx.vs.CIdiff.pdf",width = 4.5, height = 3)



# 6g ========================================================

dat <- fread("Greek_data.txt",data.table = F)
head(dat)


dat$Timepoint <- factor(dat$Timepoint, levels = c("Pre","Post"))

dat <- dat %>% arrange(Timepoint)

ggpaired(dat, x = "Timepoint", y = "CI",
         color = "Timepoint", 
         line.color = "gray", line.size = 0.4,
         point.size = 2)+
  stat_compare_means(paired = TRUE) +
  ylab("BDI")+
  scale_color_manual(values = c("indianred","lightgray"))
ggsave(filename = "6g.BDI_pre.vs.post.pdf", width = 3.5, height = 4)

ggpaired(dat, x = "Timepoint", y = "Pa",
         color = "Timepoint",
         line.color = "gray", line.size = 0.4, 
         point.size = 2)+
  stat_compare_means(paired = TRUE)  +
  ylab("Pa") +
  scale_color_manual(values = c("indianred","lightgray"))
ggsave(filename = "6g.Pa_pre.vs.post.pdf", width = 3.5, height = 4)
