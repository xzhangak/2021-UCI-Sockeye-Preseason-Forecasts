mydata=read.table("KE_BRTB_ADJ.csv", sep = ",", header=T)

#forecast return R1_2 
kenai12<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kenai sockeye R1_2 forecast 
  #input data
  #d=read.table("KE_BRTB_ADJ.csv", sep = ",", header=T)
  myvars <- c("YEAR","ESCP","R1_2")
  d=d[myvars] #select the variables
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  d$R12S<-d$R1_2/d$ESCP
  
  #Create data for brood interaction Model: MODEL LNRS=INTAC;
  mylag <- function(x,k) c(rep(NA,k),head(x,-k)) #lag function
  nrs <- d$ESCP # equivalent to datalines
  one <- data.frame(
    x = nrs,
    lag1 = mylag(nrs,1),
    lag2 = mylag(nrs,2),
    year = 999  # R automatically loops, so no extra command needed
  )
  d$lag1 <- one$lag1
  d$INTAC=d$ESCP*d$lag1
  
  #Kenai sockeye R1_2 forecast using escapment: LNR1_2=LNESCP   
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
  
  #Kenai sockeye R1_2 forecast using Ricker model: LNRS=ESCP  
  fish=d
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=11 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit12.rick <- lm(log(R12S) ~ ESCP, data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit12.rick, new)
    pred<-exp(lnpred)*new$ESCP #calculate predicted R1_2 by R/S * S
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  rick=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #Kenai sockeye R1_2 forecast using Brood interaction Model: LNRS=ESCP *ESCP(-1) (LNRS=INTAC ) 
  fish=d
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=11 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit12.intac <- lm(log(R12S) ~ INTAC, data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit12.intac, new)
    #pred<-exp(lnpred)*$ESCP[i] #calculate predicted R1_2 by R/S * S
    pred<-exp(lnpred)*new$ESCP #calculate predicted R1_2 by R/S * S
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_2")
  obs=obs[myvars]
  intac=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
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
  #obs=fish[which(fish$YEAR > 2002),]
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
  
  accuracy(lnesc$pred, lnesc$R1_2)  #Model:LNR1_2=LNESCP 
  accuracy( rick$pred, rick$R1_2) #Ricker model:LNRS=ESCP 
  accuracy( intac$pred, intac$R1_2) #Brood interaction model:LNRS=ESCP *ESCP(-1)
  accuracy( exsmooth$forecast.point, exsmooth$R1_2) #Exponential smoothing
  accuracy( ma$forecast.point, ma$R1_2) #Moving average
  
  comb<-cbind(obs$YEAR, obs$R1_2, lnesc$pred,rick$pred, intac$pred, 
              exsmooth$forecast.point, ma$forecast.point)
  comb[,1] = comb[,1]+ 4 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_2", "lnesc", "ricker", "brd-inter", 
                    "exsmooth", "ma")
  comb
}

kenai13<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kenai sockeye R1_3 forecast 
  #input data
  #d=read.table("KE_BRTB_ADJ.csv", sep = ",", header=T) #they are same as SAS data
  #names(d)[1] <- "YEAR" #rename first colume's name
  myvars <- c("YEAR","FABND0","FSWT0","R1_2", "R1_3")
  d=d[myvars] #select the variables
  is.na(d) <- d=="." #change missing value . to NA
  #d <- na.omit(d) 
  #remove rows with missing values(uncompleted returns)
  #for model comparisons below, must use the same dataset, so removing all missing values
  #if check with Mark's SAS code's fry and sibling models, comment it out
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  d$CFSWT[d$FSWT0<0.9]<-1 #create variable CFSWT (=1 if FSWT0<0.9)
  d$CFSWT[d$FSWT0>=0.9] <- 2 #CFSWT=2 if FSWT0<0.9
  d$CFSWT <- factor(d$CFSWT)#take CFSWT as a factor or "class" in SAS
  
  #Kenai sockeye R1_3 forecast using fry(age0) abundance and weight
  myvars <- c("YEAR","FABND0","CFSWT", "R1_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and two-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    #fish.sub=fish[which(fish$YEAR < year[i]),]
    fit13.log <- lm(log(R1_3) ~ log(FABND0)+CFSWT, data=fish.sub)
    #f.yr=i+1  #forecast one-year ahead
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit13.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m, d.f)
  }
  m <- na.omit(m) 
  #obs=fish[which(fish$YEAR > 2002),]
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  fry=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  #fry <- na.omit(fry)
  
  #Kenai sockeye R1_3 forecast using R1_2 (silbing model)   
  #the best model for 2020 is sibling model:log(R1_3)~log(R1_2)
  myvars <- c("YEAR","R1_2", "R1_3")
  fish=d[myvars] #select the variables
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  #fish=fish[which(fish$R1_2 != "NA"),]
  year=fish$YEAR
  n=12 #ten-year comparison with run data and two-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit13.log <- lm(log(R1_3) ~ log(R1_2), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit13.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R1_3")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R1_3")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #ten-year comparison with run data and two-year ahead forecast
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
    f.yr=year[i] #forecast one-year ahead
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
  n=12 #ten-year comparison with run data and two-year ahead forecast
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
  
  accuracy(sibling$pred, sibling$R1_3) #calculate forecast errors
  accuracy(fry$pred, fry$R1_3)
  accuracy(exsmooth$forecast.point, exsmooth$R1_3)
  accuracy(ma$forecast.point, ma$R1_3)
  
  comb<-cbind(obs$YEAR,obs$R1_3, sibling$pred, fry$pred, 
              exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 5 #change blood year to return year
  colnames(comb)<-c("run.year", "R1_3", "sibling.log", "fry", 
                    "exsmooth", "ma")
  comb
}

kenai22<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kenai sockeye R2_2 forecast 
  #Input data
  d=read.table("KE_BRTB_ADJ.csv", sep = ",", header=T)
  #names(d)[1] <- "YEAR" #rename first colume's name
  myvars <- c("YEAR","ESCP","R2_1","R2_2")
  d=d[myvars] #select the variables
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  
  #Kenai sockeye R2_2 forecast using escapment: LNR2_2=LNESCP  
  myvars <- c("YEAR","ESCP","R2_2")
  fish<-d[myvars]
  #fish <- na.omit(d)
  year=fish$YEAR
  n=12 #ten-year comparison with run data and two-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit22.esc <- lm(log(R2_2) ~ log(ESCP), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit22.esc, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  lnesc=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #Kenai sockeye R2_2 forecast using sibling model:log(R2_2)=log(R2_1)
  myvars <- c("YEAR","R2_1","R2_2")
  fish<-d[myvars]
  #fish <- subset(fish, fish$R2_1 != 0) #remove R2_1=0
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=12 #ten-year comparison with run data and one-year ahead forecast
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit22.sibling <- lm(log(R2_2) ~ log(R2_1 + 1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit22.sibling, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_2")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R2_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12
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
  #exsmooth=na.omit(exsmooth)
  
  #simple moving avearage from library(smooth)
  myvars <- c("YEAR", "R2_2")
  fish=d[myvars] #select the variables
  year=fish$YEAR
  n=12 #how many years to forecast
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
  
  accuracy(lnesc$pred, lnesc$R2_2)
  accuracy(sibling$pred, sibling$R2_2)
  accuracy(exsmooth$forecast.point, exsmooth$R2_2)
  accuracy(ma$forecast.point, ma$R2_2)
  
  comb<-cbind(obs$YEAR,obs$R2_2, sibling$pred, lnesc$pred, 
              exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 5 #change blood year to return year
  colnames(comb)<-c("run.year", "R2_2", "sibling.log", "lnesc",  "exsmooth", "ma")
  comb
}

kenai23<-function(d=mydata){
  library(forecast)
  library(smooth)
  #Kenai sockeye R2_3 forecast 
  #input data
  #d=read.table("KE_BRTB_ADJ.csv", sep = ",", header=T)
  #names(d)[1] <- "YEAR" #rename first colume's name
  myvars <- c("YEAR","ESCP", "FABND1","R2_2", "R2_3")
  d=d[myvars] #select the variables
  is.na(d) <- d=="." #change missing value . to NA
  d<-as.data.frame(lapply(d, as.numeric)) #change all varialbes as numeric
  d.omit4 <- subset(d, YEAR!="1983"&YEAR!="1987"&YEAR!="2000"&YEAR!="2005")
  #d.omit4 <- d  #if I want to know what if without omit4, run it.
  rownames(d.omit4) <- NULL #reset rownames integer
  
  #Forecast using regression on fry(age1) abundance with ARIMA error
  #best model for 2020 forecast
  fish <- d.omit4
  #fish <- subset(fish, fish$FABND1 != 0)
  #fish <- subset(fish, fish$R2_2 != "NA")
  year=fish$YEAR
  n=13
  year<-tail(year,n)
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit23ar<-auto.arima(log(fish.sub$R2_3), xreg=log(fish.sub$FABND1+1))
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-forecast(fit23ar, xreg=log(new$FABND1 +1))
    pred<-exp(lnpred$mean)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  fry=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #Kenai sockeye R2_3 forecast using R2_2 (silbing model log transformation)  
  #log(R2_3) ~ log(R2_2)
  myvars <- c("YEAR","R2_2", "R2_3")
  fish=d.omit4[myvars] #select the variables
  #fish <- subset(fish, fish$R2_2 != "NA")
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=13
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit13.log <- lm(log(R2_3) ~ log(R2_2+1), data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    lnpred<-predict(fit13.log, new)
    pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  sibling.log=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  #sibling <- na.omit(sibling) 
  #accuracy(sibling.log$R2_3, sibling.log$pred)
  
  #Kenai sockeye R2_3 forecast using R2_2 (no log transformation) 
  #R2_3 ~ R2_2
  myvars <- c("YEAR","R2_2", "R2_3")
  fish=d.omit4[myvars] #select the variables
  #fish <- subset(fish, fish$R2_2 != "NA")
  options(contrasts=c("contr.SAS","contr.poly")) #make parameter vlaues equal to one from SAS. Compare to Mark's SAS
  year=fish$YEAR
  n=13
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "pred") 
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    fit13 <- lm(R2_3 ~ R2_2, data=fish.sub)
    f.yr=year[i]  #forecast one-year ahead
    new<-fish[which(fish$YEAR==f.yr),]
    pred<-predict(fit13, new)
    #pred<-exp(lnpred)
    d.f=cbind(f.yr, pred)
    m=rbind(m,d.f)
  }
  m <- na.omit(m) 
  obs=fish[which(fish$YEAR >= year[1]),]
  myvars <- c("YEAR", "R2_3")
  obs=obs[myvars]
  sibling=merge(obs, m, by.x="YEAR", by.y="f.yr", all=T)
  
  #forecast using exponential smoothing
  myvars <- c("YEAR", "R2_3")
  fish=d.omit4[myvars] #select the variables
  #fish <- na.omit(fish) #remove rows with missing values(uncompleted returns)
  year=fish$YEAR
  n=13
  year<-tail(year,n) #get the last 10 years for comparison
  m=matrix(NA, nrow=1, ncol=2)
  colnames(m)=c("f.yr", "forecast.point")
  m=as.data.frame(m)
  for(i in 1:n ) {
    fish.sub=fish[which(fish$YEAR < year[i]),]
    run=fish.sub$R2_3
    #exponential smoothing from library(smooth)
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
  fish=d.omit4[myvars] #select the variables
  #fish <- na.omit(fish) #remove rows with missing values(uncompleted returns)
  year=fish$YEAR
  n=13
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
  #ma=na.omit(ma)
  
  
  accuracy(sibling$pred,sibling$R2_3) #calculate forecast errors
  accuracy( sibling.log$pred, sibling.log$R2_3) #calculate forecast errors
  accuracy(fry$pred, fry$R2_3)
  accuracy(exsmooth$forecast.point, exsmooth$R2_3)
  accuracy(ma$forecast.point, ma$R2_3)
  
  comb<-cbind(obs$YEAR,obs$R2_3, sibling.log$pred, sibling$pred, fry$pred,
              exsmooth$forecast.point,ma$forecast.point)
  comb[,1] = comb[,1]+ 6 #change blood year to return year
  colnames(comb)<-c("run.year", "R2_3", "sibling.log", "sibling", "fry",  "exsmooth", "ma")
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
myforecast<-function(s, d, type.err ='RMSE', run.y=2021 ) {  
  mod.name<-s[which(s$type==type.err), 2] #if choose RMSE, the model
  d<-as.data.frame(d) #must be dataframe to do next coding line
  f<-d[which(d$run.year== run.y), mod.name] #forecast for 2021 using selected model
  f
}

# yr=2021
run.predict<-function(yr){
  #predict age12 group
  mod.out<-kenai12(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r12<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  #predict age13 group
  mod.out<-kenai13(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r13<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  #predict age22 group
  mod.out<-kenai22(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r22<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  #predict age23 group
  mod.out<-kenai23(mydata)
  merr<-mod.err(d=mod.out)
  s<-mod.select(merr=merr,d=mod.out)
  r23<-myforecast(s=s, d=mod.out, type.err ='RMSE', run.y=yr)
  
  totalr<-r12 +r13 +r22 +r23 #sum of major age classes
  totalr <- totalr #major age classes only
  totalr
}

run.predict(yr=2021) #forecast next year's run size

