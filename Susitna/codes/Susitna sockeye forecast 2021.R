setwd("H:/Upper Cook Inlet/Forecasts_RB/Preseason Forecasts/Data")
#Susitna sockeye forecast 
#input data
d=read.table("SusitnaRun2020.csv", sep = ",", skip = 1, header=T)
attach(d)
#run by age class (proportions from genetic analaysis of CF harvest)
d$R02=P02*Run; d$R11=P11*Run; d$R03=P03*Run 
d$R12=P12*Run; d$R21=P21*Run; d$R04=P04*Run
d$R13=P13*Run; d$R22=P22*Run; d$R31=P31*Run
d$R14=P14*Run; d$R23=P23*Run; d$R32=P32*Run
d$R24=P24*Run; d$R33=P33*Run
d$subrun<-d$R03 + d$R12 + d$R13 + d$R22 + d$R23
d$ratio<-d$Run/d$subrun #need it to calculate the total run later
detach(d)
myvars <- c("Year","Run","Spawners","R03","R12","R13","R22","R23", "ratio")
rt<-d[myvars] #Run by Age Class table
spawners<-d[c("Year","Spawners")]#retrive spawners data and year
R0312<-d[c("Year","R03","R12")] #retrieve age03 and 12 in the run
R0312$Year = R0312$Year - 4 #change the run year to blood year related to age03 and 12
comb1=merge(spawners, R0312, by="Year", all=T) #build the blood-year table

R1322<-d[c("Year","R13","R22")] #retrieve the age classes
R1322$Year = R1322$Year - 5 #change the run year to blood year related to ages
comb2=merge(comb1, R1322, by="Year", all=T) #build the blood-year table

R23<-d[c("Year","R23")] #retrieve the age class
R23$Year = R23$Year - 6 #change the run year to blood year related to ages
comb3=merge(comb2, R23, by="Year", all=T) #build the blood-year table

bt<-comb3 #the final blood-year table
bt$R03S<-bt$R03/bt$Spawners #R03 per spawner
bt$R12S<-bt$R12/bt$Spawners #R12 per spawner
bt$R13S<-bt$R13/bt$Spawners #R13 per spawner
bt$R22S<-bt$R22/bt$Spawners #R22 per spawner
bt$R23S<-bt$R23/bt$Spawners #R23 per spawner

#forecast of 2021 return by Age Class:r03~r23 by
#multiplying spawners by average of retrun-per-spawner at that age class
f.yr = 2021
r03=bt$Spawners[which(bt$Year==f.yr-4)]*mean(bt$R03S, na.rm=T)#na.rm=T to remove NA's
r12=bt$Spawners[which(bt$Year==f.yr-4)]*mean(bt$R12S, na.rm=T)
r13=bt$Spawners[which(bt$Year==f.yr-5)]*mean(bt$R13S, na.rm=T)
r22=bt$Spawners[which(bt$Year==f.yr-5)]*mean(bt$R22S, na.rm=T)
r23=bt$Spawners[which(bt$Year==f.yr-6)]*mean(bt$R23S, na.rm=T)
subtotal<-r03+r12+r13+r22+r23
ave.ratio = mean(d$ratio, na.rm=T) #average of ratios(total_run/subtotal)
total.run<-subtotal*ave.ratio
f2021<-cbind(f.yr,r03,r12,r13,r22,r23,subtotal, total.run)
print(f2021, digits=5)

#forecast of 2022 return by Age Class:r03~r23 by
#multiplying spawners by average of retrun-per-spawner at that age class
f.yr = 2022
r03=bt$Spawners[which(bt$Year==f.yr-4)]*mean(bt$R03S, na.rm=T)#na.rm=T to remove NA's
r12=bt$Spawners[which(bt$Year==f.yr-4)]*mean(bt$R12S, na.rm=T)
r13=bt$Spawners[which(bt$Year==f.yr-5)]*mean(bt$R13S, na.rm=T)
r22=bt$Spawners[which(bt$Year==f.yr-5)]*mean(bt$R22S, na.rm=T)
r23=bt$Spawners[which(bt$Year==f.yr-6)]*mean(bt$R23S, na.rm=T)
subtotal<-r03+r12+r13+r22+r23
ave.ratio = mean(d$ratio, na.rm=T) #average of ratios(total_run/subtotal)
total.run<-subtotal*ave.ratio
f2022<-cbind(f.yr,r03,r12,r13,r22,r23, subtotal, total.run)
print(f2022, digits=5)

#Write results to a file
path <- "H:/Upper Cook Inlet/Forecasts_RB/Preseason Forecasts/Results/Susitna"
write.csv(bt, file.path(path, "brood-yr table.csv"))
write.csv(f2021, file.path(path, "run forecast.csv"))


