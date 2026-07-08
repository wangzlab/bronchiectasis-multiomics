
library(data.table)
library(ggplot2)
library(dplyr)

dat <- fread("glm_clinicalVars_vs_CI.Shannon.csv", data.table = F)

head(dat)
dat$Estimate[dat$X == "Shannon"] <- -dat$Estimate[dat$X == "Shannon"]  #只把Shannon的Association符号改变
# dat$Estimate <- -dat$Estimate

table(dat$Group)

dat$`-logP` <-  -log10(dat$p_value)
dat$`-logP*Direction` <-  -log10(dat$p_value) * sign(dat$Estimate)

dat$dataSet <- factor(dat$dataSet, levels = c("All","EMBARC","CAMEB2","China"))

dat$var %>% unique

dat <- dat %>%
  filter(var %in% c("BSI","FEV1pred_perc","AE_1year","MMRC"))
dat$var <- factor(dat$var, levels = rev(c("BSI","FEV1pred_perc","AE_1year","MMRC")))

dat$sigLabels0.1 <- ifelse(dat$p_value < 0.1, ifelse(dat$Estimate>0 , "+","-"), "")
dat$sigLabels0.05 <- ifelse(dat$p_value < 0.05, ifelse(dat$Estimate>0 , "+","-"), "")

# 系列颜色 ----------------
library(scales)
color_gradient <- colorRampPalette(c("white",  "#79bdbf")) # 使用 colorRampPalette() 函数生成渐变色函数
greens <- color_gradient(10) # # 生成渐变色向量, 这里的10表示生成10个颜色
show_col(greens);greens

color_gradient <- colorRampPalette(c("white",   "#e4a0a4")) # 使用 colorRampPalette() 函数生成渐变色函数
reds <- color_gradient(10) # # 生成渐变色向量, 这里的10表示生成10个颜色
show_col(reds);reds


# 分All和Commensal画 --------------------------------------
# all 
dat.all <- dat %>% filter(Group == "All")
CIvalues <- dat.all$`-logP*Direction`[dat.all$X == "CI"]
Shannonvalues <- dat.all$`-logP*Direction`[dat.all$X == "Shannon"]

# 大于0 的数值
CIvalues[CIvalues > 0] %>% summary
Shannonvalues[Shannonvalues > 0] %>% summary

CIvalues[CIvalues > 0][order(CIvalues[CIvalues > 0])]
Shannonvalues[Shannonvalues > 0][order(Shannonvalues[Shannonvalues > 0])]

# 小于0 的数值
CIvalues[CIvalues < 0] %>% summary
Shannonvalues[Shannonvalues < 0] %>% summary

CIvalues[CIvalues < 0][order(CIvalues[CIvalues < 0])]
Shannonvalues[Shannonvalues < 0][order(Shannonvalues[Shannonvalues < 0])]


Values.all <- c(-14.0773528, -10.2768319, -7.552528, 0,  3.66676 , 5.48775  , 8.62052  )
Colors.all <- c(greens[c(10,8,4)],"white", reds[c(5,7,10)])

ggplot(dat.all) +
  geom_tile(aes(x=dataSet, y=var, fill = `-logP*Direction`), color="#E3E3E3") +
  facet_grid(Group~X) +
  scale_fill_gradientn(values = scales::rescale(Values.all),
                       colours = Colors.all) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(aes(x=dataSet, y=var, label = sigLabels0.05))
ggsave(filename = "6c.GroupAll.pdf", width = 6, height = 2)

# Commensal 
dat.commensal <- dat %>% filter(Group == "Commensal")
CIvalues <- dat.commensal$`-logP*Direction`[dat.commensal$X == "CI"]
Shannonvalues <- dat.commensal$`-logP*Direction`[dat.commensal$X == "Shannon"]


# 大于0 的数值
CIvalues[CIvalues > 0] %>% summary
Shannonvalues[Shannonvalues > 0] %>% summary

CIvalues[CIvalues > 0][order(CIvalues[CIvalues > 0])]
Shannonvalues[Shannonvalues > 0][order(Shannonvalues[Shannonvalues > 0])]


# 小于0 的数值
CIvalues[CIvalues < 0] %>% summary
Shannonvalues[Shannonvalues < 0] %>% summary

CIvalues[CIvalues < 0][order(CIvalues[CIvalues < 0])]
Shannonvalues[Shannonvalues < 0][order(Shannonvalues[Shannonvalues < 0])]


Values.commensal <- c(-3.7922934, -1.7456, -0.6493,  0,   2.79824, 3.5,  5.5336 )
Colors.commensal <- c(greens[c(10,8,5)], "white", reds[c(6, 8, 10)])


ggplot(dat.commensal) +
  geom_tile(aes(x=dataSet, y=var, fill = `-logP*Direction`), color="#E3E3E3") +
  facet_grid(Group~X) +
  scale_fill_gradientn(values = scales::rescale(Values.commensal),
                       colours = Colors.commensal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text(aes(x=dataSet, y=var, label = sigLabels0.05))
ggsave(filename = "6c.GroupNPD.pdf", width = 6, height = 2)
