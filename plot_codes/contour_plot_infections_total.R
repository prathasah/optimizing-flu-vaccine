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
dt11 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v1.csv", header=T)
dt12 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v2.csv", header=T)
dt13 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v3.csv", header=T)
dt14 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v4.csv", header=T)
dt15 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v5.csv", header=T)
dt16 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v6.csv", header=T)
dt17 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v7.csv", header=T)
dt18 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v8.csv", header=T)
dt19 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v9.csv", header=T)
dt110 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v10.csv", header=T)
dt111 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v11.csv", header=T)
dt112 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v12.csv", header=T)
dt113 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v13.csv", header=T)
dt114 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v14.csv", header=T)
dt115 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v15.csv", header=T)
dt116 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v16.csv", header=T)
dt117 <- read.csv("../data/full_model_dec18_500iter/Incidence_11Dec2017_v17.csv", header=T)


dt1 <- rbind(dt11, dt12,dt13,dt14,dt15,dt16,dt17,dt18,dt19,dt110, dt111,dt112,dt113,dt114,dt115,dt116,dt117)
dt1$infection_per_million <- dt1$total_infections/1000000
dt1 <- subset(dt1, dt1$vaccine_efficacy>=0.1)
dt1 <- subset(dt1, dt1$vaccine_efficacy <= 0.6)
dt1 <- subset(dt1, dt1$relative_coverage < 1)
dt1$vaccine_efficacy <- dt1$vaccine_efficacy*100
dt1$relative_coverage <- dt1$relative_coverage*100
dt1hat <- dt1[,c("relative_coverage", "vaccine_efficacy", "infection_per_million")]
dt1_dups <- dt1hat[,c("relative_coverage", "vaccine_efficacy")]
dt1 <- dt1[ order(dt1["vaccine_efficacy"], dt1["relative_coverage"]), ]
write.csv(dt1,"check_infection_fromR.csv")
dt1 <- dt1hat[!duplicated(dt1_dups),]
table(dt1$vaccine_efficacy)
table(dt1$relative_coverage)
dt1$relative_coverage <- as.numeric(dt1$relative_coverage)
dt1$vaccine_efficacy <- as.numeric(dt1$vaccine_efficacy)


table(dt1$vaccine_efficacy)
summary(factor(dt1$relative_coverage))
dt21 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v1.csv", header=T)
dt22 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v2.csv", header=T)
dt23 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v3.csv", header=T)
dt24 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v4.csv", header=T)
dt25 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v5.csv", header=T)
dt26 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v6.csv", header=T)
dt27 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v7.csv", header=T)
dt28 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v8.csv", header=T)
dt29 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v9.csv", header=T)
dt210 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v10.csv", header=T)
dt211 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v11.csv", header=T)
dt212 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v12.csv", header=T)
dt213 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v13.csv", header=T)
dt214 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v14.csv", header=T)
dt215 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v15.csv", header=T)
dt216 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v16.csv", header=T)
dt217 <- read.csv("../data/full_model_dec18_500iter/Mortality_11Dec2017_v17.csv", header=T)

dt2 <- rbind(dt21, dt22,dt23,dt24,dt25,dt26,dt27,dt28,dt29,dt210 , dt211,dt212,dt213,dt214,dt215,dt216,dt217)
dt2 <- subset(dt2, dt2$vaccine_efficacy>=0.1)
dt2 <- subset(dt2, dt2$vaccine_efficacy <= 0.6)
dt2$mortality_per_thousand <- dt2$total_mortality/1000
dt2$vaccine_efficacy <- dt2$vaccine_efficacy*100
dt2$relative_coverage <- dt2$relative_coverage*100
dt2hat <- dt2[,c("relative_coverage", "vaccine_efficacy", "mortality_per_thousand")]
dt2_dups <- dt2hat[,c("relative_coverage", "vaccine_efficacy")]
dt2 <- dt2hat[!duplicated(dt2_dups),]
dt2 <- dt2[ order(dt2["vaccine_efficacy"], dt2["relative_coverage"]), ]

######
## pull some numbers for the paper

## infections, high efficacy vaccine
base.infection <- dt1[(dt1$relative_coverage=="54") & (dt1$vaccine_efficacy=="50"),]$infection_per_million
high.coverage.infection<- dt1[(dt1$relative_coverage=="60") & (dt1$vaccine_efficacy=="50"),]$infection_per_million
perc.decrease <- ((base.infection - high.coverage.infection)*100)/base.infection
perc.decrease


### for low efficacy
base.infection <- dt1[(dt1$relative_coverage=="54") & (dt1$vaccine_efficacy=="15"),]$infection_per_million
high.coverage.infection<- dt1[(dt1$relative_coverage=="60") & (dt1$vaccine_efficacy=="15"),]$infection_per_million
((base.infection - high.coverage.infection)*100)/base.infection

### mortality, high efficacy vaccine
base.mortality <- dt2[(dt2$relative_coverage=="54") & (dt2$vaccine_efficacy=="50"),]$mortality_per_thousand
high.coverage.mortality<- dt2[(dt2$relative_coverage=="60") & (dt2$vaccine_efficacy=="50"),]$mortality_per_thousand
((base.mortality - high.coverage.mortality)*100)/base.mortality


### for low efficacy
base.mortality <- dt2[(dt2$relative_coverage=="54") & (dt2$vaccine_efficacy=="15"),]$mortality_per_thousand
high.coverage.mortality<- dt2[(dt2$relative_coverage=="60") & (dt2$vaccine_efficacy=="15"),]$mortality_per_thousand
((base.mortality - high.coverage.mortality)*100)/base.mortality
#####################################################################
## some pretty colors
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

str(dt1)
summary(factor(dt1$vaccine_efficacy))
summary(factor(dt1$relative_coverage))

ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = infection_per_million))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("Overall coverage (%)")+ 
  scale_fill_gradientn(colours=my.cols, breaks=c(0, 25,50, 100, 150,200,250))+
#  theme(panel.background=element_rect(fill="blue", colour="blue"))+
  theme(legend.title=element_blank())+ggtitle("Influenza infections (millions)")+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))


#ggsave("infections_2017_10iter.png", p1, height = 5, width =7)


p2 <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = mortality_per_thousand))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("Overall coverage (%)")+ 
  scale_fill_gradientn(colours=my.cols, breaks=c(0, 100, 200, 300, 400))+
  theme(legend.title=element_blank())+ggtitle("Influenza mortality (thousands)")

p2
#ggsave("hospitalization_2017_10iter.png", p2, height = 5, width =7)

library(gridExtra)
p <- grid.arrange(p1,p2, nrow=1)

ggsave("effect_of_efficacy_coverage.png", p,height=4.2, width=11)


