library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(lattice)
library(reshape)
library(RColorBrewer)   # for brewer.pal(...)
library(cowplot)
library(Hmisc)
################################################################
dt11 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v1.csv", header=T)
dt12 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v2.csv", header=T)
dt13 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v3.csv", header=T)
dt14 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v4.csv", header=T)
dt15 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v5.csv", header=T)
dt16 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v6.csv", header=T)
dt17 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v7.csv", header=T)
dt18 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v8.csv", header=T)
dt19 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v9.csv", header=T)
dt110 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v10.csv", header=T)
dt111 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v11.csv", header=T)
dt112 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v12.csv", header=T)
dt113 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v13.csv", header=T)
dt114 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v14.csv", header=T)
dt115 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v15.csv", header=T)
dt116 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v16.csv", header=T)
dt117 <- read.csv("../data/full_model_dec18_500iter/vaccinated_comparisons_11Dec2017_v17.csv", header=T)
dt <- rbind(dt11, dt12,dt13,dt14,dt15,dt16,dt17,dt18,dt19,dt110, dt111,dt112,dt113,dt114,dt115,dt116,dt117)
dt <- subset(dt, dt$vaccine_efficacy>=0.1)
dt <- subset(dt, dt$vaccine_efficacy <=0.6)
dt <- subset(dt, dt$relative_coverage < 1)
dt$relative_coverage <- dt$relative_coverage*100
dt$vaccine_efficacy <- dt$vaccine_efficacy *100
dt <- subset(dt, dt$vaccine_efficacy>9)
dt_dups <- dt[c("relative_coverage", "vaccine_efficacy")]
dt <- dt[!duplicated(dt_dups),]
dt <- dt[ order(dt["vaccine_efficacy"], dt["relative_coverage"]), ]


dt$vac_prob_infection <- dt$vaccinated_total_infections/dt$total_vaccinated
dt$unvac_prob_infection <- dt$unvaccinated_total_infections/dt$total_unvaccinated
dt$OR.vax <- dt$vac_prob_infection/dt$unvac_prob_infection
dt$OR.unvax <- dt$unvac_prob_infection/dt$vac_prob_infection
### plot odds ration

dt1 <- subset(dt, dt$relative_coverage == 10)
dt2 <- subset(dt, dt$relative_coverage == 20)
dt3 <- subset(dt, dt$relative_coverage == 30)
dt4 <- subset(dt, dt$relative_coverage == 40)
dt5 <- subset(dt, dt$relative_coverage == 50)
dt6 <- subset(dt, dt$relative_coverage == 60)
dt7 <- subset(dt, dt$relative_coverage == 70)
dt8 <- subset(dt, dt$relative_coverage == 80)

dx <- rbind(dt2, dt4,dt6, dt8)
summary(factor(dx$vaccine_efficacy))
dx$relative_coverage <- factor(dx$relative_coverage)

p1 <- ggplot(dx, aes(x= vaccine_efficacy, y = OR.vax , color=relative_coverage))+
  geom_point()+ 
  #geom_smooth(method = "lm") + 
  ylab("Odds ratio for infection of vaccinated")+xlab("Vaccine efficacy (%)")+
  theme(panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
         text = element_text(size=12),legend.key = element_blank(), legend.title = element_text(size=10))+
 labs(color='Baseline coverage (%)')+
  theme(legend.key = element_blank())+
  theme(legend.position =c(0.08, 0.2))+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))
p1

dx2 <- rbind(dt4,dt5)
summary(lm(unvac_prob_infection ~ vaccine_efficacy, data=dt2))
summary(lm(unvac_prob_infection ~ vaccine_efficacy, data=dt4))
summary(lm(unvac_prob_infection ~ vaccine_efficacy, data=dt6))
summary(lm(unvac_prob_infection ~ vaccine_efficacy, data=dt8))

p2 <- ggplot(dx, aes(x= vaccine_efficacy, y = unvac_prob_infection,  color=relative_coverage))+
  geom_point()+ 
  #geom_smooth(method = "lm") + 
  theme(panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=12),legend.key = element_blank(), legend.title = element_text(size=10))+
  ylab("Probability of infection for unvaccinated")+
  theme(legend.position = "none")+
  xlab("Vaccine efficacy (%)")+
  scale_x_continuous(breaks=c(10, 20, 30, 40, 50, 60))

#ggsave("herd_immunity_10iter.png", p2, height=4, width = 5)

library(gridExtra)
p <- grid.arrange(p1, p2, nrow=1)
ggsave("../plots/individual_level_risks_500iter.png",p, height= 4, width=8)

