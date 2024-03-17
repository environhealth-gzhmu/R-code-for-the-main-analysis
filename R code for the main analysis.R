##############################################################################
#  "Drivers of associations between daytime-nighttime compound temperature 
#   extremes and mortality in China"
#   Jun Yang, et al.
#   Communications Medicine, 2024
#   Date: 16th March, 2024
#   Email: yangjun_eci@jnu.edu.cn; yangjun@gzhmu.edu.cn
##############################################################################
###1. Install and load the packages 
install.packages(c("splines","mgcv","nlme","dlnm","tsModel","metafor")
library(splines);library(nlme);library(mgcv);library(dlnm);library(tsModel);
library(MASS);library(metafor)

##############################################################################
##2. Load the data
datas<-as.data.frame(get(load("daily_temperature_mortality.Rdata")))
code<-unique(datas$city_code) #"city_code" is the ID code for each city
datas$date<-as.Date(datas$date,origin="1960-1-1")
datas$dow<-format(datas$date,"%w") 
datas$year<-format(datas$date,"%Y")
datas$month<-format(datas$date,"%m")
datas$tmean<-datas$mean_temperature #"mean_temperature" is daily mean temperature.
datas$tmax<-datas$maximum_temperature #"maximum_temperature" is daily maximum temperature.
datas$tmin<-datas$minimum_temperature #"minimum_temperature" is daily minimum temperature.
dlist<-lapply(code,function(x) datas[datas$CODE==x,])

##############################################################################
###3.1 Creat functions for identify heat wave and cold spell
#3.1.1 Function fun.hw(x,dur,percent) to define heat wave 
#"x" is temperature series; "dur" and "percent" are the duration and threshold of heat wave;
fun.hw <- function(x,dur,percent){
s.hw<-apply(Lag(as.numeric(x>=quantile(x,percent,na.rm=T)),0:(dur-1)),1,sum,na.rm=T)
s.hw[s.hw<dur]<-0
for (i in (dur:length(s.hw))){
  if (s.hw[i]==dur){
      s.hw[c((i-dur+1):(i-1))]<-1
}
}
hw<-ifelse(s.hw>0,1,0)
return(hw)
}

#3.1.2 Function fun.cs() to define cold spell 
#"x" is temperature series; "dur" and "percent" are the duration and threshold of cold spell;
fun.cs <- function(x,dur,percent){
s.cs<-apply(Lag(as.numeric(x<=quantile(x,percent,na.rm=T)),0:(dur-1)),1,sum,na.rm=T)
s.cs[s.cs<dur]<-0
for (i in (dur:length(s.cs))){
  if (s.cs[i]==dur){
      s.cs[c((i-dur+1):(i-1))]<-1
}
}
cs<-ifelse(s.cs>0,1,0)
return(cs)
}

###3.2 Creat functions for calculation of the Q-AIC
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

##############################################################################
##4.1 Estimate the effect of heat wave at lag 0-1
#Use definition of at least 3 days with daily maximum temperature ¡Ý 95th for heat wave as example
#Data during hot season (May to September) are used for heat wave
HW.est01<-HW.est.se01<-HW.est02<-HW.est.se02<-HW.est03<-HW.est.se03<-rep(NA,length(dlist))
for (i in seq(dlist)){
  sub1<-dlist[[i]]
  sub<-subset(sub1,sub1$month%in%c(5,6,7,8,9))
  duration<-3;percent<-0.95;
  sub$HW01<-fun.hw(sub$tmax,duration,percent)
  sub$HW02<-fun.hw(sub$tmin,duration,percent)
  sub$HW03<-sub$HW01-sub$HW02;sub$HW03[sub$HW03<0]<-0 #"HW03" is used for defining  heat wave only during daytime
  sub$HW04<-sub$HW02-sub$HW01;sub$HW04[sub$HW04<0]<-0 #"HW04" is used for defining heat wave only during nighttime
  sub$HW05<-sub$HW02+sub$HW01;sub$HW05[sub$HW05<2]<-0
  sub$HW05[sub$HW05==2]<-1  #"HW05" is used for defining heat wave occurred both at daytime and nighttime                        
  sub$HW03<-runMean(sub$HW03,0:1) #Moving average of heat wave at current day and previous day
  sub$HW04<-runMean(sub$HW04,0:1)
  sub$HW05<-runMean(sub$HW05,0:1)
  knots<-equalknots(sub$hum,fun="ns",df=5)
  logknots<-logknots(10,df=4)  
  argvar=list(knots=knots);arglag=list(knots=logknots)
  hum.basis<-crossbasis(sub$hum,argvar=argvar,arglag=arglag,lag=10) #Defining crossbasis for relative humidity using maximum lag of 10 days
  press.basis<-crossbasis(sub$press,argvar=list(knots=equalknots(sub$press,fun="ns",df=5)),arglag=arglag,lag=10) #Defining crossbasis for air pressure using maximum lag of 10 days
  #Calculate effect of heat wave occurred only during daytime for each city and store it in "est01", and its standard error in "est.se01"
  model01<-glm(death ~HW03+ns(date,df=4)+hum.basis+as.factor(year)+press.basis+as.factor(dow)+as.factor(holiday),family=quasipoisson,data=sub);
  HW.est01[k]<-summary(model01)$coefficients[2,1];HW.est.se01[k]<-summary(model01)$coefficients[2,2]
  #Calculate effect of heat wave occurred only during nighttime for each city and store it in "est02", and its standard error in "est.se02"
  model02<-glm(death ~HW04+ns(date,df=4)+hum.basis+as.factor(year)+press.basis+as.factor(dow)+as.factor(holiday),family=quasipoisson,data=sub);
  HW.est02[k]<-summary(model02)$coefficients[2,1];HW.est.se02[k]<-summary(model02)$coefficients[2,2]
  #Calculate effect of heat wave occurred both at daytime and nighttime for each city and store it in "est03", and its standard error in "est.se03"
  model03<-glm(death ~HW05+ns(date,df=4)+hum.basis+as.factor(year)+press.basis+as.factor(dow)+as.factor(holiday),family=quasipoisson,data=sub);
  HW.est03[k]<-summary(model03)$coefficients[2,1];HW.est.se03[k]<-summary(model03)$coefficients[2,2]
}
#rma.uni function is used to pool city-specific estimates
model01<-rma.uni(yi=HW.est01,sei=HW.est.se01,method="REML");model02<-rma.uni(yi=HW.est02,sei=HW.est.se02,method="REML");model03<-rma.uni(yi=HW.est03,sei=HW.est.se03,method="REML");
#Calculate the percentage change and its 95%CI for heat wave at daytime-only, nigttime-only and those both at daytime and nighttime
paste0(round(exp(as.numeric(model01$b))-1)*100,2),"(",round(exp(as.numeric(model01$ci.lb))-1)*100,2),",",round(exp(as.numeric(model01$ci.ub))-1)*100,2),")");
paste0(round(exp(as.numeric(model02$b))-1)*100,2),"(",round(exp(as.numeric(model02$ci.lb))-1)*100,2),",",round(exp(as.numeric(model02$ci.ub))-1)*100,2),")");
paste0(round(exp(as.numeric(model03$b))-1)*100,2),"(",round(exp(as.numeric(model03$ci.lb))-1)*100,2),",",round(exp(as.numeric(model03$ci.ub))-1)*100,2),")");

###
##4.2 Estimate the effect of cold spell at lag 0-14
#Data during cold season (November to March) are used for cold spell
#Use definition of at least 3 days with daily minimum temperature ¡Ü 5th for cold spell as example
CS.est01<-CS.est.se01<-CS.est02<-CS.est.se02<-CS.est03<-CS.est.se03<-rep(NA,length(dlist))
for (i in seq(dlist)){
  sub1<-dlist[[i]]
  sub<-subset(sub1,sub1$month%in%c(11,12,1,2,3))
  duration<-3;percent<-0.05;
  sub$CS01<-fun.cs(sub$tmax,duration,percent)
  sub$CS02<-fun.cs(sub$tmin,duration,percent)
  sub$CS03<-sub$CS01-sub$CS02;sub$CS03[sub$CS03<0]<-0 #"CS03" is used for defining  cold spell only during daytime
  sub$CS04<-sub$CS02-sub$CS01;sub$CS04[sub$CS04<0]<-0 #"CS04" is used for defining cold spell only during nighttime
  sub$CS05<-sub$CS02+sub$CS01;sub$CS05[sub$CS05<2]<-0
  sub$CS05[sub$CS05==2]<-1  #"CS05" is used for defining cold spell occurred both at daytime and nighttime                        
  sub$CS03<-runMean(sub$CS03,0:14) #Moving average of cold spell at lag 0-14 days
  sub$CS04<-runMean(sub$CS04,0:14)
  sub$CS05<-runMean(sub$CS05,0:14)
  knots<-equalknots(sub$hum,fun="ns",df=5)
  logknots<-logknots(10,df=4)  
  argvar=list(knots=knots);arglag=list(knots=logknots)
  hum.basis<-crossbasis(sub$hum,argvar=argvar,arglag=arglag,lag=10) #Defining crossbasis for relative humidity using maximum lag of 10 days
  press.basis<-crossbasis(sub$press,argvar=list(knots=equalknots(sub$press,fun="ns",df=5)),arglag=arglag,lag=10) #Defining crossbasis for air pressure using maximum lag of 10 days
  #Calculate effect of cold spell occurred only during daytime for each city and store it in "est01", and its standard error in "est.se01"
  model01<-glm(death ~CS03+ns(date,df=4)+hum.basis+as.factor(year)+press.basis+as.factor(dow)+as.factor(holiday),family=quasipoisson,data=sub);
  CS.est01[k]<-summary(model01)$coefficients[2,1];CS.est.se01[k]<-summary(model01)$coefficients[2,2]
  #Calculate effect of cold spell occurred only during nighttime for each city and store it in "est02", and its standard error in "est.se02"
  model02<-glm(death ~CS04+ns(date,df=4)+hum.basis+as.factor(year)+press.basis+as.factor(dow)+as.factor(holiday),family=quasipoisson,data=sub);
  CS.est02[k]<-summary(model02)$coefficients[2,1];CS.est.se02[k]<-summary(model02)$coefficients[2,2]
  #Calculate effect of cold spell occurred both at daytime and nighttime for each city and store it in "est03", and its standard error in "est.se03"
  model03<-glm(death ~CS05+ns(date,df=4)+hum.basis+as.factor(year)+press.basis+as.factor(dow)+as.factor(holiday),family=quasipoisson,data=sub);
  CS.est03[k]<-summary(model03)$coefficients[2,1];CS.est.se03[k]<-summary(model03)$coefficients[2,2]
}
#rma.uni function is used to pool city-specific estimates
model01<-rma.uni(yi=CS.est01,sei=CS.est.se01,method="REML");model02<-rma.uni(yi=CS.est02,sei=CS.est.se02,method="REML");model03<-rma.uni(yi=CS.est03,sei=CS.est.se03,method="REML");
#Calculate the percentage change and its 95%CI for cold spell at daytime-only, nigttime-only and those both at daytime and nighttime
paste0(round(exp(as.numeric(model01$b))-1)*100,2),"(",round(exp(as.numeric(model01$ci.lb))-1)*100,2),",",round(exp(as.numeric(model01$ci.ub))-1)*100,2),")");
paste0(round(exp(as.numeric(model02$b))-1)*100,2),"(",round(exp(as.numeric(model02$ci.lb))-1)*100,2),",",round(exp(as.numeric(model02$ci.ub))-1)*100,2),")");
paste0(round(exp(as.numeric(model03$b))-1)*100,2),"(",round(exp(as.numeric(model03$ci.lb))-1)*100,2),",",round(exp(as.numeric(model03$ci.ub))-1)*100,2),")");

##############################################################################
#5. Meta-regression analysis on city-level characteristics
#5.1 Explore the influence of city-level characteristics on mortality risk of heat wave
#Take number of population as moderator on the heat wave occurred both at daytime and nighttime
socioeconomic_data<-read.csv("socioeconomic data for each city.csv")
y<-socioeconomic_data$Population #Number of population for each city
HW.moderator.pop<-rep(NA,6) #"moderator.pop" is used to store result for meta-regression
model<-rma.uni(yi=HW.est03,sei=HW.est.se03,mods=cbind(y),method="REML")
HW.moderator.popr[1]<-IQR(y1,na.rm=T);
HW.moderator.pop[2]<-paste(round((exp(summary(model)$b[2,1])-1)*100,2),"(",round((exp(summary(model)$ci.lb[2])-1)*100,2),",",
round((exp(summary(model)$ci.ub[2])-1)*100,2),")",sep="")
HW.moderator.pop[3]<-round(summary(model)$pval[2],3)
HW.moderator.pop[4:6]<-c(round(summary(model)$tau2,6),round(summary(model)$I2,6),round(summary(model)$R2,6))

#5.1 Explore the influence of city-level characteristics on mortality risk of cold spell
#Take number of population as moderator on the cold spell occurred both at daytime and nighttime
CS.moderator.pop<-rep(NA,6) #"moderator.pop" is used to store result for meta-regression
model<-rma.uni(yi=CS.est03,sei=CS.est.se03,mods=cbind(y),method="REML")
CS.moderator.pop[1]<-IQR(y1,na.rm=T);
CS.moderator.pop[2]<-paste(round((exp(summary(model)$b[2,1])-1)*100,2),"(",round((exp(summary(model)$ci.lb[2])-1)*100,2),",",
round((exp(summary(model)$ci.ub[2])-1)*100,2),")",sep="")
CS.moderator.pop[3]<-round(summary(model)$pval[2],3)
CS.moderator.pop[4:6]<-c(round(summary(model)$tau2,6),round(summary(model)$I2,6),round(summary(model)$R2,6))


