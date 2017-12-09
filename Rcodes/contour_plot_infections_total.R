library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(lattice)
library(reshape)
library(RColorBrewer)   # for brewer.pal(...)
library(cowplot)
library(Hmisc)
library(plotly)
library(plot3D)
################################################################
dt1 <- read.csv("../data/Incidence_8Dec2017.csv", header=T)
dt1$infection_per_million <- dt1$total_infections/1000000
dt1 <- dt1[,c("relative_coverage", "vaccine_efficacy", "infection_per_million")]

dt2 <- read.csv("../data/Hospitalizations_8Dec2017.csv", header=T)
str(dt2)

dh <- read.csv("../data/historical_seasons_data.csv", header=T)
dh$infection_per_million <- dh$infections
str(dh)
dh<- na.omit(dh)
#####################################################################
## some pretty colors
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))



#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
p1 <- ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = infection_per_million))+
  geom_tile() + xlab("Vaccine efficacy") + ylab("Overall coverage")+ 
  geom_vline(xintercept = 0.1, color ="black")+
  scale_fill_gradientn(colours=my.cols, breaks=c(10,50, 100,  150))

p1
ggsave("infections_2017.eps", p1, height = 5, width =8)
#p1 <- p1 +  geom_point(data=dh, aes(x=vaccine_efficacy, y=relative_coverage, fill=(infection_per_million)), pch=21)

library(reshape2)
df <-acast(dt1, vaccine_efficacy~ relative_coverage, value.var="infection_per_million")
vaccine_efficacy <- seq(0.001, 0.601, by=0.01)
relative_coverage <- seq(0,0.99, by=0.01)

persp3D(x=vaccine_efficacy, y=relative_coverage, z=df, 
        xlab = "Vaccine efficacy", ylab="Vaccination coverage", zlab="Infections")







dt2$hosp_permillion <- dt2$total_hospitalizations/1000000
p2 <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = hosp_permillion))+
  geom_tile() + xlab("Vaccine efficacy") + ylab("Relative coverage (as compared to typical data)")+
  scale_fill_gradientn(colours=my.cols)

p2
ggsave("hospitalization_2017.eps", p1, height = 5, width =8)

#p<- grid.arrange(p1, p2,  p3, p4, nrow=2, ncol=2)
#p
#cairo_ps("Qthreshold_vs_pathogen_robustness.eps",  height =3, width = 3.42)
#p<- grid.arrange(p1, p2,  p3, p4, nrow=2, ncol=2)
#dev.off()
#ggsave("3_Contour_disease_burden.eps", p, height = 7.5, width = 8.7, units = "cm", dpi = 1200)
