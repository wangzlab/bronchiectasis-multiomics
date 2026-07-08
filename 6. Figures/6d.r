library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


dat<- fread("external_CI_plot_violin6.txt", data.table = F)
head(dat)


dat$disease <- factor(dat$disease, levels = c("Health", "COPD", "Asthma", "Bronchiectasis", "CF", "Lung cancer","Lung transplantation","Pulmonary TB"))

ggplot(dat, aes(x=disease, y=BDI))  +
  geom_violin(aes( fill=disease), alpha=0.6, trim = F) +
  geom_boxplot( width=0.1, outlier.shape = NA) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Set3"))(22)) +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_text(angle = 90))  +      
  facet_wrap(vars(BioProject), scales="free", nrow=1 ) +
  stat_compare_means(label = "p.signif", method =  "wilcox.test",
                     ref.group = "Health",  # manually change to "GeneralHealth" and "NoSymptom" and "largeCAT" and revise stars accordingly
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                        symbols = c("***", "**", "*", "+", "ns")))               
ggsave(filename = "6d.externalData.pdf", width = 7, height = 3)
