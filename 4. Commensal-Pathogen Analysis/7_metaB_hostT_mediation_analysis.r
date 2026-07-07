list.files()

files <- list.files("metab_hostt_mediation_table/", full.names = T)

library(data.table)
library(mediation)

for(f in files){
  # f =  files[1]
  
  data<-read.table(f, sep="\t",header=T)
  model.m=lm(HostT~MetaB+Age+Gender+BMI,data)
  model.y=lm(Pa_Hi~HostT+MetaB+Age+Gender+BMI,data)
  summary=summary(mediate(model.m,model.y,treat="MetaB",mediator="HostT",boot=F,sims=1000))
  capture.output(summary,file=paste0("forward/", basename(f)),append=FALSE)
  
  
  model.m=lm(HostT~Pa_Hi+Age+Gender+BMI,data)
  model.y=lm(MetaB~HostT+Pa_Hi+Age+Gender+BMI,data)
  summary=summary(mediate(model.m,model.y,treat="Pa_Hi",mediator="HostT",boot=F,sims=1000))
  capture.output(summary,file=paste0("reverse/", basename(f)),append=FALSE)
  
}
