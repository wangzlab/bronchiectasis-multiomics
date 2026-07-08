list.files()

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

dat<- fread("ABO_BCO_COPD_Asthma.txt", data.table = F)
head(dat)

dat$disease <- factor(dat$disease, levels = c("Healthy", "BCO", "ABO", "COPD", "Asthma", "Bronchiectasis"))

ggplot(dat, aes(x=disease, y=CI))  +
  #geom_hline(yintercept = 0, linetype="dashed", color="gray", lwd=1) +
  geom_violin(aes( fill=disease), alpha=0.6, trim = F) +
  geom_boxplot( width=0.1, outlier.shape = NA) +
  #scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Set3"))(22)) +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_text(angle = 90))  +      
 stat_compare_means(label = "p.signif", method =  "wilcox.test",
                    ref.group = "Healthy",  
                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                       symbols = c("***", "**", "*", "+", "ns")))               

ggsave(filename = "6e.violin.pdf", width = 4.5, height = 3)
