list.files()

files <- list.files("metat_metab_mediation_table/", full.names = T)

library(data.table)
library(mediation)

for(f in files){
  # f =  files[1]
  
  data<-read.table(f, sep="\t",header=T)
  model.m=lm(MetaB~MetaT+Age+Gender+BMI,data)
  model.y=lm(Pa_Hi~MetaT+MetaB+Age+Gender+BMI,data)
  summary=summary(mediate(model.m,model.y,treat="MetaT",mediator="MetaB",boot=F,sims=1000))
  capture.output(summary,file=paste0("forward/", basename(f)),append=FALSE)
  
  
  model.m=lm(MetaB~Pa_Hi+Age+Gender+BMI,data)
  model.y=lm(MetaT~MetaB+Pa_Hi+Age+Gender+BMI,data)
  summary=summary(mediate(model.m,model.y,treat="Pa_Hi",mediator="MetaB",boot=F,sims=1000))
  capture.output(summary,file=paste0("reverse/", basename(f)),append=FALSE)
  
}
