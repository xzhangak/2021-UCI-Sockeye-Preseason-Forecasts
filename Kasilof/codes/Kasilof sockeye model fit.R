mydata=read.table("KA_BRTB_ADJ.csv", sep = ",", header=T)

kasilof12<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kasilof sockeye R1_2 forecast 
  #input data
  #d=read.table("KA_BRTB_ADJ.csv", sep = ",", header=T)
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  myvars <- c("YEAR","ESCP","SM1","R1_1", "R1_2")
  d=d[myvars] #select the variables
  
  #R1_2 using log(R1_1) (silbing model with log transformation)   
  myvars <- c("YEAR","R1_1", "R1_2")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  #fish=fish[which(fish$R1_2 != "NA"),]
  year=fish$YEAR
  n=11 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit12.log <- lm(log(R1_2) ~ log(R1_1 + 1), data=fish.sub) # R1_1 + 1 to avoid log(0)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit12.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #R1_2 forecast using escapment: LNR1_2=LNESCP   
  myvars <- c("YEAR","ESCP","R1_2")
  fish<-d[myvars]
  year=fish$YEAR
  n=11#ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit12.esc <- lm(log(R1_2) ~ log(ESCP), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit12.esc, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  lnesc=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #Forecast using regression on escapment with ARIMA error
  myvars <- c("YEAR","ESCP","R1_2")
  fish<-d[myvars]
  year=fish$YEAR
  n=11
  year<-tail(year,n)
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit22ar<-auto.arima(log(fish.sub$R1_2), xreg=log(fish.sub$ESCP))
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-forecast(fit22ar, xreg=log(new$ESCP))
    pred<-exp(lnpred$mean)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  escp.ar=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #R1_2 forecast using smolt model  
  myvars <- c("YEAR","SM1", "R1_2")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=11 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit12 <- lm(log(R1_2) ~ log(SM1), data=fish.sub) 
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit12, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  smolt=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R1_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=11 #ten-year comparison with run data and one-year ahead forecast
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
  n=11 #ten-year comparison with run data and one-year ahead forecast
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
  
  comb<-cbind(obs$YEAR,obs$R1_2,sibling$pred, lnesc$pred, escp.ar$pred,
              smolt$pred, exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 4 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_2", "sibling.log", "lnesc", "lnesc.ar",
                    "smolt", "exsmooth", "ma")
  comb
}

kasilof13<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kasilof sockeye R1_3 forecast 
  #input data
  #d=read.table("KA_BRTB_ADJ.csv", sep = ",", header=T)
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  myvars <- c("YEAR","SM1","R1_2", "R1_3")
  d=d[myvars] #select the variables
  
  #R1_3 using log(R1_3)~log(R1_2) (silbing model)   
  myvars <- c("YEAR","R1_2", "R1_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit.log <- lm(log(R1_3) ~ log(R1_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #R1_3 forecast using smolt model:log(R1_2)~log(SM1)
  myvars <- c("YEAR","SM1", "R1_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm(log(R1_3) ~ log(SM1), data=fish.sub) 
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
  smolt=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R1_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
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
  
  #simple moving avearage from library("smooth")
  myvars <- c("YEAR", "R1_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
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
  
  comb<-cbind(obs$YEAR,obs$R1_3,sibling$pred, 
              smolt$pred, exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 5 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_3", "sibling.log", 
                    "smolt", "exsmooth", "ma")
  comb
}

kasilof22<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kasilof sockeye R2_2 forecast 
  #input data
  #d=read.table("KA_BRTB_ADJ.csv", sep = ",", header=T)
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  myvars <- c("YEAR","ESCP","R2_1", "R2_2")
  d=d[myvars] #select the variables
  
  #using log(R2_2)~log(R2_1) (silbing model)   
  myvars <- c("YEAR","R2_1", "R2_2")
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
    fit.log <- lm(log(R2_2) ~ log(R2_1+1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #using log(R2_2)~log(ESCP)   
  myvars <- c("YEAR","ESCP", "R2_2")
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
    fit.log <- lm(log(R2_2) ~ log(ESCP), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  esc=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R2_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
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
  n=12 #ten-year comparison with run data and one-year ahead forecast
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
  
  comb<-cbind(obs$YEAR,obs$R2_2, sibling$pred, esc$pred, 
              exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 5 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_3", "sibling.log", "lnesc",  "exsmooth", "ma")
  comb
}

kasilof23<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kasilof sockeye R2_3 forecast 
  #input data
  #d=read.table("KA_BRTB_ADJ.csv", sep = ",", header=T)
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  myvars <- c("YEAR","SM2","R2_2", "R2_3")
  d=d[myvars] #select the variables
  
  #R2_3 using log(R2_3)~log(R2_2) (silbing model)   
  myvars <- c("YEAR","R2_2", "R2_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit.log <- lm(log(R2_3) ~ log(R2_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #R2_3 using R2_3~R2_2 (silbing model without log transformation)   
  myvars <- c("YEAR","R2_2", "R2_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit<- lm(R2_3 ~ R2_2, data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    pred<-predict(fit, new)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  sibling2=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #R2_3 forecast using smolt model  
  myvars <- c("YEAR","SM2", "R2_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit <- lm(log(R2_3) ~ log(SM2), data=fish.sub) 
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
  smolt=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R2_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=13 #ten-year comparison with run data and one-year ahead forecast
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
  n=13 #ten-year comparison with run data and one-year ahead forecast
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
  
  comb<-cbind(obs$YEAR,obs$R2_3, sibling$pred,sibling2$pred, smolt$pred, 
              exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 6 #change blood year to return year
  colnames(comb)<-c("run.year", "R2_3", "sibling.log", "sibling", "smolt", 
                    "exsmooth", "ma")
  comb
}

#function to calculate forecast errors
mod.err<-function(d){ #d is the model predictions and actual runs
  #d<-kasilof12
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
myforecast<-function(s, d, type.err ='RMSE') {  
  mod.name<-s[which(s$type==type.err), 2] #if RMSE, best model is selected
  d<-as.data.frame(d) #must be dataframe to do next coding line
  comb<-cbind(d[,1:2],d[,mod.name])
  names(comb)[3] <- mod.name #rename the column 3
  comb
}

#choose a type of error for model selection
#RMSE, MAPE, MAE, MPE, and ME
m=matrix(NA, nrow=1, ncol=4)
colnames(m)=c("run.year", "obs","fitted.all", "err.type")
m=as.data.frame(m)
err.group<-c("RMSE", "MAPE", "MAE", "MPE", "ME")
for(er in err.group){
  myerror<-er
  mod.out<-kasilof12(mydata)
  mod.out<-mod.out[1:10,] #the most recent 10-year fitted vaules with actual runs
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  fit.r12<-myforecast(s=s, d=mod.out, type.err = myerror)
  #predict age13 group
  mod.out<-kasilof13(mydata)
  mod.out<-mod.out[1:10,]#the most recent 10-year fitted vaules with actual runs
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  fit.r13<-myforecast(s=s, d=mod.out, type.err =myerror)
  #predict age22 group
  mod.out<-kasilof22(mydata)
  mod.out<-mod.out[1:10,]#the most recent 10-year fitted vaules with actual runs
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  fit.r22<-myforecast(s=s, d=mod.out, type.err =myerror)
  #predict age23 group
  mod.out<-kasilof23(mydata)
  mod.out<-mod.out[1:10,]#the most recent 10-year fitted vaules with actual runs
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  fit.r23<-myforecast(s=s, d=mod.out, type.err =myerror)
  
  obs<-fit.r12[,2]+fit.r13[,2]+fit.r22[,2]+fit.r23[,2]#actual runs
  fitted.all<-fit.r12[,3]+fit.r13[,3]+fit.r22[,3]+fit.r23[,3]#fitted vaule
  run.year<-fit.r12[,1]#retrieve the firt column (run.year)
  
  d<-cbind.data.frame(run.year, obs, fitted.all)
  d$err.type<-myerror
  m=rbind(m,d)
} 
staking <- na.omit(m) 
#staking all fitted values including 10-yr actual run and fitted with different Min.errors 

#Write results to a file
path <- "H:/Upper Cook Inlet/Forecasts_RB/Preseason Forecasts/Results/Kasilof"
write.csv(staking, file.path(path, "model fitted.csv"))


#ggplot: model fitted vs. actual run
library(tidyverse)
library(scales)
mytitle <- "Model Fitted value v.s. Actual Run(black dots)"
myplot<-ggplot(data = staking) + 
  geom_line(mapping = aes(x = run.year, y = fitted.all, colour=err.type ))+
  geom_point(mapping = aes(x = run.year, y = obs ))+
  labs(
    title = mytitle,
    x="Year",
    y="Run Size",
    colour ="Error Type")+   #don't show legend title
  theme(legend.position="right",
        legend.justification = c("right", "top"))+
  scale_x_continuous(breaks= pretty_breaks())#library(scales) needed

print(myplot)

