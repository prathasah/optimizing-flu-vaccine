library(plotly)


dt1 <- read.csv("../data/Incidence_8Dec2017.csv", header=T)
dt1$infection_per_million <- dt1$total_infections/1000000
dt1 <- dt1[,c("relative_coverage", "vaccine_efficacy", "infection_per_million")]


df <-acast(dt1,  relative_coverage ~vaccine_efficacy, value.var="infection_per_million")
vaccine_efficacy <- seq(0.001, 0.601, by=0.01)
relative_coverage <- seq(0,0.99, by=0.01)
##check dimension
dim(df)
length(vaccine_efficacy)
length(relative_coverage)


plot_ly(type = 'surface' , x=relative_coverage, y=vaccine_efficacy, z = (df))%>%
  layout(
  scene = list(
    xaxis=list(title="Vaccination coverage"),
    yaxis=list(title="Vaccine efficacy"),
    zaxis=list(title ="Infections (millions)")))
