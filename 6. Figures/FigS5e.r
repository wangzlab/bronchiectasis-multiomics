dat <- fread("FigS5e.txt", data.table = F)
ggplot(dat) +
geom_point(aes(x=dat$Pa, y=dat$Vp), color="#E6949A") +
geom_smooth(aes(x=x=dat$Pa, y=dat$Vp), linetype="dashed", method = "lm", color="black", se = F) +
theme_classic()
ggplot(dat) +
geom_point(aes(x=dat$Pa, y=dat$Vd), color="#E6949A") +
geom_smooth(aes(x=x=dat$Pa, y=dat$Vd), linetype="dashed", method = "lm", color="black", se = F) +
theme_classic()
ggplot(dat) +
geom_point(aes(x=dat$Pa, y=dat$Pm), color="#E6949A") +
geom_smooth(aes(x=x=dat$Pa, y=dat$Pm), linetype="dashed", method = "lm", color="black", se = F) +
theme_classic()
ggplot(dat) +
geom_point(aes(x=dat$Hi, y=dat$Rm), color="#E6949A", size=2.5, shape=15) +
geom_smooth(aes(x=dat$Hi, y=dat$Rm), linetype="dashed", method = "lm", color="black", se = F) +
theme_classic()