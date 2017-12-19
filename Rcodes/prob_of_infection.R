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
dta <- read.csv("../data/vaccinated_comparisons_8Dec2017.csv", header=T)
dtb <- read.csv("../data/vaccinated_comparisons_8Dec2017_v2.csv", header=T)
dt <- rbind(dta,dtb)
str(dt)

dt$vac_prob_infection <- dt$vaccinated_total_infections/dt$total_vaccinated
dt$unvac_prob_infection <- dt$unvaccinated_total_infections/dt$total_unvaccinated
dt$OR.vax <- dt$vac_prob_infection/dt$unvac_prob_infection
dt$OR.unvax <- dt$unvac_prob_infection/dt$vac_prob_infection
### plot odds ration

dt1 <- subset(dt, dt$relative_coverage == 0.1)
dt2 <- subset(dt, dt$relative_coverage == 0.2)
dt3 <- subset(dt, dt$relative_coverage == 0.3)
dt4 <- subset(dt, dt$relative_coverage == 0.4)
dt5 <- subset(dt, dt$relative_coverage == 0.5)
dt6 <- subset(dt, dt$relative_coverage == 0.6)
dt7 <- subset(dt, dt$relative_coverage == 0.7)
dt8 <- subset(dt, dt$relative_coverage == 0.8)

dx <- rbind(dt2,dt4,dt6,dt8)
dx$relative_coverage <- factor(dx$relative_coverage)
p1 <- ggplot(dx, aes(x= vaccine_efficacy, y = OR.vax,  color = relative_coverage))+
  geom_point()+ geom_smooth(method = "loess") +ylim(0.5,1)+
  ylab("Odds ratio of being infected for vaccinated")+xlab("Vaccine efficacy")
p1
ggsave("odds_ratio_vaccinated_being_infected_50iter.pdf", p1, height=6, width = 8)

dx2 <- rbind(dt3,dt4,dt5)
m1 <- lm(unvac_prob_infection ~ vaccine_efficacy, data=dx2)
summary(m1)
p2 <- ggplot(dx2, aes(x= vaccine_efficacy, y = unvac_prob_infection))+
  geom_point()+ geom_smooth(method = "lm") + ylim(0,1)+
  ylab("P(infection) for unvaccinated individuals")+
  xlab("Vaccine efficacy")
p2
ggsave("herd_immunity_50iter.pdf", p2, height=4, width = 5)

dt1 <- dt[,c("relative_coverage", "vaccine_efficacy", "vac_prob_infection")]
dt2 <- dt[,c("relative_coverage", "vaccine_efficacy", "unvac_prob_infection")]
dt1$type <- "vaccinated"
dt2$type <- "unvaccinated"

colnames(dt1)[3] <- "prob_infection"
colnames(dt2)[3] <- "prob_infection"

df <- rbind(dt1, dt2)

df1 <- subset(df, df$relative_coverage == 0.1)
df2 <- subset(df, df$relative_coverage == 0.2)
df3 <- subset(df, df$relative_coverage == 0.3)
df4 <- subset(df, df$relative_coverage == 0.4)
df5 <- subset(df, df$relative_coverage == 0.5)
df6 <- subset(df, df$relative_coverage == 0.6)
df7 <- subset(df, df$relative_coverage == 0.7)
df8 <- subset(df, df$relative_coverage == 0.8)

de <- rbind(df1,df2,df3,df4,df5,df6,df7,df8)

p1 <- ggplot(de, aes(x=vaccine_efficacy, y=prob_infection, color = type))+
    theme_bw()+  geom_point() + geom_line()+
  xlab("vaccine efficacy")+ ylab("Probability of infection")+ facet_grid(relative_coverage~.)

