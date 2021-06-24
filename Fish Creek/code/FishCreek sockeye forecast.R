setwd("H:/Upper Cook Inlet/Forecasts_RB/Preseason Forecasts/Data")
mydata=read.table("fishcrkBYTB2020.csv", sep = ",", skip = 1, header=T)
names(mydata)[1] <- "YEAR"

fishcrk12<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Fish Creek sockeye forecast 
  #input data
  #d=read.table("fishcrkBYTB2020.csv", sep = ",", skip = 1, header=T)
  #names(d)[1] <- "YEAR"

  #using log(R1_2)~log(R1_1) (silbing model) 
  myvars <- c("YEAR","R1_1", "R1_2")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  #fish=fish[which(fish$R1_2 != "NA"),]
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm(log(R1_2) ~ log(R1_1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #using (R1_2)~(R1_1) (silbing model without log transformation)   
  myvars <- c("YEAR","R1_1", "R1_2")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm((R1_2) ~ (R1_1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  sibling2=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R1_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R1_2
    #exponential smoothing from library("smooth")
    mod.es=es(run, holdout = F)
    f.es=forecast(mod.es, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.es$forecast)
    lower90=as.vector(f.es$lower)
    upper90=as.vector(f.es$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.es=cbind(f.yr, forecast.point)
    d.es=as.data.frame(d.es)
    m=rbind(m,d.es)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  exsmooth=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #simple moving avearage from library("smooth")
  myvars <- c("YEAR", "R1_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R1_2
    mod.sma= sma(run)
    f.sma=forecast(mod.sma, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.sma$forecast)
    lower90=as.vector(f.sma$lower)
    upper90=as.vector(f.sma$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.sma=cbind(f.yr, forecast.point)
    d.sma=as.data.frame(d.sma)
    m=rbind(m,d.sma)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  ma=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  accuracy(sibling2$pred, sibling2$R1_2)
  accuracy(exsmooth$forecast.point, exsmooth$R1_2)
  accuracy(ma$forecast.point, ma$R1_2)
  accuracy(sibling$pred, sibling$R1_2)
  
  comb<-cbind(obs$YEAR,obs$R1_2,sibling$pred, sibling2$pred,exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 4 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_2", "sibling.log", "sibling", "exsmooth", "ma")
  comb
}

fishcrk13<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Fish Creek sockeye forecast 
  #input data
  #d=read.table("fishcrkBYTB2020.csv", sep = ",", skip = 1, header=T)
  #names(d)[1] <- "YEAR"
  
  #using log(R1_3)~log(R1_2) (silbing model) 
  myvars <- c("YEAR","R1_2", "R1_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm(log(R1_3) ~ log(R1_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #using (R1_3)~(R1_2) (silbing model without log transformation)   
  myvars <- c("YEAR","R1_2", "R1_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm((R1_3) ~ (R1_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  sibling2=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R1_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R1_3
    #exponential smoothing from library("smooth")
    mod.es=es(run, holdout = F)
    f.es=forecast(mod.es, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.es$forecast)
    lower90=as.vector(f.es$lower)
    upper90=as.vector(f.es$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.es=cbind(f.yr, forecast.point)
    d.es=as.data.frame(d.es)
    m=rbind(m,d.es)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  exsmooth=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  #exsmooth=na.omit(exsmooth)
  
  #simple moving avearage from library("smooth")
  myvars <- c("YEAR", "R1_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R1_3
    mod.sma= sma(run)
    f.sma=forecast(mod.sma, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.sma$forecast)
    lower90=as.vector(f.sma$lower)
    upper90=as.vector(f.sma$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.sma=cbind(f.yr, forecast.point)
    d.sma=as.data.frame(d.sma)
    m=rbind(m,d.sma)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  ma=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  accuracy(sibling$pred, sibling$R1_3)
  accuracy(sibling2$pred, sibling2$R1_3)
  accuracy(exsmooth$forecast.point, exsmooth$R1_3)
  accuracy(ma$forecast.point, ma$R1_3)
  comb<-cbind(obs$YEAR,obs$R1_3,sibling$pred, sibling2$pred,exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 5 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_3", "sibling.log", "sibling", "exsmooth", "ma")
  comb
}  

fishcrk22<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Fish Creek sockeye forecast 
  #input data
  #d=read.table("fishcrkBYTB2020.csv", sep = ",", skip = 1, header=T)
  #names(d)[1] <- "YEAR"
  
  #using log(R2_2)~log(R2_1) (silbing model) 
  myvars <- c("YEAR","R2_1", "R2_2")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm(log(R2_2) ~ log(R2_1+1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #using (R2_2)~(R2_1) (silbing model without log transformation)   
  myvars <- c("YEAR","R2_1", "R2_2")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm((R2_2) ~ (R2_1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  sibling2=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R2_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R2_2
    #exponential smoothing from library("smooth")
    mod.es=es(run, holdout = F)
    f.es=forecast(mod.es, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.es$forecast)
    lower90=as.vector(f.es$lower)
    upper90=as.vector(f.es$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.es=cbind(f.yr, forecast.point)
    d.es=as.data.frame(d.es)
    m=rbind(m,d.es)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  exsmooth=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #simple moving avearage from library("smooth")
  myvars <- c("YEAR", "R2_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R2_2
    mod.sma= sma(run)
    f.sma=forecast(mod.sma, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.sma$forecast)
    lower90=as.vector(f.sma$lower)
    upper90=as.vector(f.sma$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.sma=cbind(f.yr, forecast.point)
    d.sma=as.data.frame(d.sma)
    m=rbind(m,d.sma)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  ma=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  accuracy(sibling$pred, sibling$R2_2)
  accuracy(sibling2$pred, sibling2$R2_2)
  accuracy(exsmooth$forecast.point, exsmooth$R2_2)
  accuracy(ma$forecast.point, ma$R2_2)
  
  comb<-cbind(obs$YEAR,obs$R2_2,sibling$pred, sibling2$pred,exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 5 #change blood year to return year
  colnames(comb)<-c("run.year", "R2_2", "sibling.log", "sibling", "exsmooth", "ma")
  comb
}

fishcrk23<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Fish Creek sockeye forecast 
  #input data
  #d=read.table("fishcrkBYTB2020.csv", sep = ",", skip = 1, header=T)
  #names(d)[1] <- "YEAR"
  
  #using log(R2_3)~log(R2_2) (silbing model) 
  myvars <- c("YEAR","R2_2", "R2_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=14 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm(log(R2_3+1) ~ log(R2_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #using (R2_3)~(R2_2) (silbing model without log transformation)   
  myvars <- c("YEAR","R2_2", "R2_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly"))
  year=fish$YEAR
  n=14 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm((R2_3) ~ (R2_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit, new)
    pred<-(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  sibling2=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R2_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=14 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R2_3
    #exponential smoothing from library("smooth")
    mod.es=es(run, holdout = F)
    f.es=forecast(mod.es, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.es$forecast)
    lower90=as.vector(f.es$lower)
    upper90=as.vector(f.es$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.es=cbind(f.yr, forecast.point)
    d.es=as.data.frame(d.es)
    m=rbind(m,d.es)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  exsmooth=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #simple moving avearage from library("smooth")
  myvars <- c("YEAR", "R2_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=14 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R2_3
    mod.sma= sma(run)
    f.sma=forecast(mod.sma, h=1, level =0.9) #one-year ahead forecast
    forecast.point=as.vector(f.sma$forecast)
    lower90=as.vector(f.sma$lower)
    upper90=as.vector(f.sma$upper)
    f.yr=year[i]  #forecast one-year ahead
    d.sma=cbind(f.yr, forecast.point)
    d.sma=as.data.frame(d.sma)
    m=rbind(m,d.sma)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  ma=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  accuracy(sibling$pred, sibling$R2_3)
  accuracy(sibling2$pred, sibling2$R2_3)
  accuracy(exsmooth$forecast.point, exsmooth$R2_3)
  accuracy(ma$forecast.point, ma$R2_3)
  
  comb<-cbind(obs$YEAR,obs$R2_3,sibling$pred, sibling2$pred,exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 6 #change blood year to return year
  colnames(comb)<-c("run.year", "R2_3", "sibling.log", "sibling", "exsmooth", "ma")
  comb
}

#function to calculate forecast errors
mod.err<-function(d){ #d is the model predictions and actual runs
  #d<-fishcrk12
  d10<-d[1:10,] #10 years compared to actual runs.
  n<-ncol(d10)
  mod.names<-colnames(d10)[-1:-2] #list of forecast methods only
  
  merr=matrix(NA, nrow=1, ncol=5)
  colnames(merr)=c("ME", "RMSE","MAE","MPE","MAPE")
  merr=as.data.frame(merr)
  d10[,2] <-d10[,2] +1 # in case actual run is 0
  #calculate forecast errors
  for(i in 3:n ) { 
    a<-accuracy(d10[,i], d10[,2])
    merr=rbind(merr,a)
  }
  merr <- na.omit(merr) 
  rownames(merr)<-mod.names #rename rownames as forecast method
  merr
}

#function to find best model by measures, like RMSE, MAPE...
mod.select<-function (merr,d){
  n<-ncol(merr)
  best<-rep(NA, n) #find best model by minumum error
  for (i in 1:n){
    min.err<-min(merr[,i])
    b<-merr[which(merr[,i]==min.err),]
    best[i]<-rownames(b) #name of best model
  }
  s<-cbind(colnames(merr),best) #error type and its best model
  s<-as.data.frame(s)
  colnames(s)<-c("type", "best") #best model for each error type
  s
}  

#function to forecast a run with selected model
#d is mod.out file; s is best model file
myforecast<-function(s, d, type.err ='RMSE', run.y=2021 ) {  
  mod.name<-s[which(s$type==type.err), 2] #if choose RMSE, the model
  d<-as.data.frame(d) #must be dataframe to do next coding line
  f<-d[which(d$run.year== run.y), mod.name] #forecast for 2021 using selected model
  f
}

#mod.out<-fishcrk12(mydata)
#mod.out<-fishcrk13(mydata)
#mod.out<-fishcrk22(mydata)
#mod.out<-fishcrk23(mydata)
yr<-2021
run.predict<-function(yr){
  #predict age12 group
  mod.out<-fishcrk12(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r12<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  #predict age13 group
  mod.out<-fishcrk13(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r13<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  #predict age22 group
  mod.out<-fishcrk22(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r22<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  #predict age23 group
  mod.out<-fishcrk23(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r23<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  
  totalr<-r12 +r13 +r22 +r23 #sum of major age classes
  totalr <-1.111 * totalr #expand to include other age classes
  totalr
}

run.predict(2021)


#Write results to a file
path <- "H:/Upper Cook Inlet/Forecasts_RB/Preseason Forecasts/Results/Fish Creek"
write.csv(mod.out, file.path(path, "R1.2 model predictions.csv"), row.names=FALSE)
write.csv(merr, file.path(path, "R1.2 model errors.csv"))
write.csv(s, file.path(path, "R1.2 best model by errors.csv"))


