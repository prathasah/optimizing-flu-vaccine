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
## demography
pop.school <- 57934551
pop.young.adults <- 83852575
pop.middle.adults <- 36585866
pop.old.adults <- 59695905
pop.elderly <- 49244195

###############3
dt1<- read.csv("../data/vaccination_results/full_model_dec19_100iter/combined_incidence_11Dec2017.csv", header=T)
dt1$infection_per_million <- dt1$total_infections/1000000
dt1 <- subset(dt1, dt1$vaccine_efficacy>=0.1)
dt1 <- subset(dt1, dt1$vaccine_efficacy <=0.6)
dt1 <- subset(dt1, dt1$relative_coverage < 1)
dt1$vaccine_efficacy <- dt1$vaccine_efficacy*100
dt1$relative_coverage <- dt1$relative_coverage*100
dt1$school_children <- (100*(dt1$age5 + dt1$age10 + dt1$age15))/pop.school 
dt1$young_adults <- (100*(dt1$age20 + dt1$age25+ dt1$age30 + dt1$age35))/pop.young.adults
dt1$middle_adults <- (100*(dt1$age40 + dt1$age45))/pop.middle.adults
dt1$old_adults <- (100*(dt1$age50 + dt1$age55))/pop.old.adults
dt1$elderly <- (100*(dt1$age60 + dt1$age65 + dt1$age70 + dt1$age75))/pop.elderly
#dplyr::filter(dt1,is.na(dt1$elderly)) 

dt2 <- read.csv("../data/vaccination_results/full_model_dec19_100iter/combined_mortality_11Dec2017.csv", header=T)
dt2 <- subset(dt2, dt2$vaccine_efficacy>=0.1)
dt2$mortality_per_thousand <- dt2$total_mortality/1000
dt2$vaccine_efficacy <- dt2$vaccine_efficacy*100
dt2$relative_coverage <- dt2$relative_coverage*100
#dt2 <- dt2[,c("relative_coverage", "vaccine_efficacy", "mortality_per_thousand")]
dt2$school_children <- (100*(dt2$age5 + dt2$age10 + dt2$age15))/pop.school 
dt2$young_adults <- (100*(dt2$age20 + dt2$age25+ dt2$age30 + dt2$age35))/pop.young.adults
dt2$middle_adults <- (100*(dt2$age40 + dt2$age45))/pop.middle.adults
dt2$old_adults <- (100*(dt2$age50 + dt2$age55))/pop.old.adults
dt2$elderly <- (100*(dt2$age60 + dt2$age65 + dt2$age70 + dt2$age75))/pop.elderly

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
library(viridis)
k <- 9
my.cols <- rev(brewer.pal(k, "BrBG"))
jet.colors <- rich.colors(100)

summary(factor(dt1$relative_coverage))
summary(factor(dt1$vaccine_efficacy))
p1a <- ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = school_children))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("Overall coverage (%)")+ 
  scale_fill_viridis(option = "magma", limits= c(0,80))+
  #scale_fill_gradientn(colours=jet.colors, breaks=c(0, 20, 40,60,80), limits= c(0,80))+
  theme(legend.title=element_blank())+ggtitle("School-aged children (5-19 years)")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(axis.text.y = element_text(size=20))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))


p1b <- ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = young_adults))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  scale_fill_viridis(option = "magma", limits = c(0,80))+
  #scale_fill_gradientn(colours=jet.colors, breaks=c(0,20, 40,60,80), limits= c(0,80))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme(axis.title.y = element_text(size=20))+
  theme(legend.title=element_blank())+ggtitle("Young adults (20-39 years)")+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))
p1b

p1c <- ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = middle_adults))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  theme(axis.title.y=element_blank())+
  scale_fill_viridis(option = "magma", limits = c(0,80))+
  #scale_fill_gradientn(colours=jet.colors, breaks=c(0,20,  40,60,80), limits= c(0,80))+
  theme(legend.title=element_blank())+ggtitle("Middle-aged adults (40-59 years)")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))

p1d <- ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = old_adults))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  theme(axis.title.y=element_blank())+
  scale_fill_viridis(option = "magma", limits = c(0,80))+
  #scale_fill_gradientn(colours=jet.colors, breaks=c(0,20,  40,60,80), limits= c(0,80))+
  theme(legend.title=element_blank())+ggtitle("Older-aged adults (50-64 years)")+
  theme(legend.position = "none")+
  theme(axis.text.y=element_blank())+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank())+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))


p1e <- ggplot(dt1, aes(x=vaccine_efficacy, y=relative_coverage, fill = elderly))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  theme(axis.title.y=element_blank())+
  scale_fill_viridis(option = "magma", limits = c(0,80), name ="Infections (%)")+
  #scale_fill_gradientn(colours=jet.colors, breaks=c(0,20,  40,60,80), limits= c(0,80),name ="Infections (%)")+
  ggtitle("Elderly (65+ years)")+
  theme(axis.title.x=element_blank())+
  theme(legend.text = element_text(size=17))+
  theme(plot.title = element_text(size=17))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme(legend.title = element_text(size=18))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 13))


###################################################################

p2a <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = school_children))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("Overall coverage (%)")+ 
  scale_fill_viridis(option = "viridis",  breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits = c(0,0.135))+
  #scale_fill_gradientn(colours=my.cols, breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits= c(0,0.33))+
 # theme(legend.title=element_blank())+ggtitle("School-aged children (5-19 years)")+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_text(size=20),
        axis.text.y = element_text(size=20))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))

p2b <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = young_adults))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  #scale_fill_gradientn(colours=my.cols, breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits= c(0,0.33))+
  scale_fill_viridis(option = "viridis",  breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125),limits = c(0,0.135))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_text(size=20))+
#  theme(legend.title=element_blank())+ggtitle("Young adults (20-29 years)")+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))

p2c <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = middle_adults))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  theme(axis.title.y=element_blank())+
  #scale_fill_gradientn(colours=my.cols, breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits= c(0,0.33))+
  scale_fill_viridis(option = "viridis", breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits = c(0,0.135))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  # theme(legend.title=element_blank())+ggtitle("Middle-aged adults (30-59 years)")+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_text(size=20))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))

p2d <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = old_adults))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  theme(axis.title.y=element_blank())+
  #scale_fill_gradientn(colours=my.cols, breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits= c(0,0.33))+
  scale_fill_viridis(option = "viridis",  breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125),limits = c(0,0.135))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
 # theme(legend.title=element_blank())+ggtitle("Middle-aged adults (30-59 years)")+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_text(size=20))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))

p2e <- ggplot(dt2, aes(x=vaccine_efficacy, y=relative_coverage, fill = elderly))+
  geom_tile() + xlab("Vaccine efficacy (%)") + ylab("")+ 
  theme(axis.title.y=element_blank())+
  scale_fill_viridis(option = "viridis", name ="Mortality (%)", breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits = c(0,0.135))+
  #scale_fill_gradientn(colours=my.cols, breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125), limits= c(0,0.33), name ="Mortality (%)")+
  #theme(legend.title=element_blank())+ggtitle("Elderly (60+ years)")+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y=element_blank())+
  theme(legend.text = element_text(size=17))+
  theme(axis.text.x=element_text(size=20),
        legend.title = element_text(size=18))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))+
 guides(fill = guide_colorbar(barwidth = 1, barheight = 13))

library(gridExtra)
p2 <- grid.arrange(p1a,p1b,p1c,p1d,p1e, p2a,p2b,p2c,p2d, p2e,nrow = 2, ncol=5, widths= c(1.1,1,1,1,1.45), heights = c(1,1.02))
g <- arrangeGrob(p2,  bottom=grid::textGrob("Vaccine efficacy (%)",gp = gpar(fontsize = 20)),
                 left=grid::textGrob("Baseline coverage (%)", rot=90, gp = gpar(fontsize = 20)))
ggsave(filename="../plots/age_specific_efficacy_coverage_100iter.png", plot=g, width = 21, height = 7)
ggsave(filename="../plots/age_specific_efficacy_coverage_100iter.pdf", plot=g, width = 21, height = 7)

#ggsave("age_specific_mortality.png", p2, width = 17, height = 7)
#grid::grid.newpage()
#grid::grid.draw(g)