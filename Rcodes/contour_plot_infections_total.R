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
dt1a <- read.csv("../data/Incidence_8Dec2017.csv", header=T)
dt1b <- read.csv("../data/Incidence_8Dec2017_v2.csv", header=T)
dt1 <- rbind(dt1a,dt1b)
dt1$infection_per_million <- dt1$total_infections/1000000
dt1 <- dt1[,c("relative_coverage", "vaccine_efficacy", "infection_per_million")]

dt2a <- read.csv("../data/Hospitalizations_8Dec2017.csv", header=T)
dt2b <- read.csv("../data/Hospitalizations_8Dec2017_v2.csv", header=T)
dt2 <- rbind(dt2a, dt2b)
dt2$hospitalization_per_million <- dt2$total_hospitalizations/1000000
dt2 <- dt2[,c("relative_coverage", "vaccine_efficacy", "hospitalization_per_million")]
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
  scale_fill_gradientn(colours=my.cols, breaks=c(0, 25,50, 100, 150,200,250))+
  theme(legend.title=element_blank())+ggtitle("Influenza cases")

p1
ggsave("infections_2017_50iter.eps", p1, height = 5, width =8)


p2 <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = hospitalization_per_million))+
  geom_tile() + xlab("Vaccine efficacy") + ylab("Overall coverage")+ 
  geom_vline(xintercept = 0.1, color ="black")+
  scale_fill_gradientn(colours=my.cols, breaks=c(0, 1,2, 3,4))+
  theme(legend.title=element_blank())+ggtitle("Case hospitalizations")

p2
ggsave("hospitalization_2017_50iter.eps", p2, height = 5, width =8)

#p1 <- p1 +  geom_point(data=dh, aes(x=vaccine_efficacy, y=relative_coverage, fill=(infection_per_million)), pch=21)





