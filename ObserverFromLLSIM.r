library(ggplot2)
library(tidyverse)
library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)
setwd("~/Box Sync/bycatch project (ebabcock@miami.edu)/LLSim")
source("3.BycatchFunctions.r")
theme_set(theme_bw())

llsets<-read.csv("catchbyset.csv",as.is=TRUE)
dim(llsets)
summary(llsets)
table(llsets$year,llsets$month)
names(llsets)

## Raw CPUE for each population 
llsets$cpue.SWO=llsets$c.SWO/llsets$hooks*1000
llsets$cpue.BUM=llsets$c.BUM/llsets$hooks*1000

#Make all the factor variables factor
llsets<-llsets %>% mutate_at(c("gear","light","hook","bait","fleet"), factor)
summary(llsets)

#Plot raw catch rates
g1<-ggplot(llsets,aes(x=year,y=cpue.SWO,color=fleet,pch=fleet))+stat_summary()
g1
g2<-ggplot(llsets,aes(x=year,y=cpue.BUM,color=fleet,pch=fleet))+stat_summary()
g2

# Allocate sets to trips 
llsets$lat5<-trunc(llsets$lat/5)*5
llsets$lon5<-trunc(llsets$lon/5)*5
table(llsets$lon5,llsets$lat5)

llsets$trip<-paste(llsets$gear,llsets$year,llsets$month,
  llsets$fleet,llsets$lat5,llsets$lon5,sep=".")
trips<-unique(llsets$trip)
length(trips)  
x<-table(llsets$trip)
length(x)
summary(as.vector(x))
hist(x)
llsets$trip1<-llsets$trip
# 20 sets per trip should be median, 100 maximum, based on US fleet
# See Gruss et al. (2019) for details. 
# Divide up "trips" more than 8 sets
a<-which(x>100)
length(a)
NumCores<-detectCores()
cl2<-makeCluster(NumCores-2)
registerDoParallel(cl2)
foreach(i = a)  %do% {
  b=which(llsets$trip %in% names(x)[i])
  d=trunc((-1+1:length(b))/20)
  llsets$trip[b]=paste(llsets$trip[b],d,sep=".")    
}
stopCluster(cl2)
trips<-unique(llsets$trip)
x<-table(llsets$trip)
summary(as.vector(x))
length(x)
ggplot(data.frame(x),aes(x=Freq))+geom_histogram()
#Plot checks that trips are reasonable 
length(unique(llsets$trip))

## Make subsets for observer coverage
#Observer coverage, 5%, 10%, 20%, random by trip and by set

trip<-sort(unique(llsets$trip))
n<-length(trip)
trip.05=sample(trip,size=n*0.05,replace=FALSE)
#trip.10=sample(trip,size=n*0.10,replace=FALSE)
#trip.20=sample(trip,size=n*0.20,replace=FALSE)
#trip.75=sample(trip,size=n*0.75,replace=FALSE)

llsets$trip.05=ifelse(llsets$trip %in% trip.05,1,0)
#llsets$trip.10=ifelse(llsets$trip %in% trip.10,1,0)
#llsets$trip.20=ifelse(llsets$trip %in% trip.20,1,0)
#llsets$trip.75=ifelse(llsets$trip %in% trip.75,1,0)

#Generate trip by trip observer dataset, at 5% coverage
obstrip<-filter(llsets,trip.05==1) %>% 
  group_by(trip) %>%
  summarize(Year=median(year),
    month=mostfreqfunc(month),
    gear=mostfreqfunc(gear),
    light=mostfreqfunc(light),
    fleet=mostfreqfunc(fleet),
    bait=mostfreqfunc(bait),
    hook=mostfreqfunc(hook),
    hooks=sum(hooks)/1000,
    sets=n(),
    SWO=sum(c.SWO),
    BUM=sum(c.BUM),
    lat5=median(lat5),
    lon5=median(lon5),
    lat=median(lat),
    lon=median(lon),
    hbf=median(hbf),
    habSWO=median(w.SWO),
    habBUM=median(w.BUM))

#Generate logbook trip by trip data, for bycatch expansion. 
logtrip<-llsets %>% 
  group_by(trip) %>%
  summarize(Year=median(year),
    month=mostfreqfunc(month),
    gear=mostfreqfunc(gear),
    light=mostfreqfunc(light),
    fleet=mostfreqfunc(fleet),
    bait=mostfreqfunc(bait),
    hook=mostfreqfunc(hook),
    hooks=sum(hooks)/1000,
    sets=n(),
    SWO=sum(c.SWO),
    BUM=sum(c.BUM),
    lat5=median(lat5),
    lon5=median(lon5),
    lat=median(lat),
    lon=median(lon),
    hbf=median(hbf),
    habSWO=median(w.SWO),
    habBUM=median(w.BUM),
    trip.05=median(trip.05))
#Check tables
head(obstrip)
head(logtrip)
dim(obstrip)
dim(logtrip)

#Make season and North vs. South Atlantic variables
logtrip<-logtrip %>%mutate(trips=1,season=factor(seasonfunc(month,4)),area=ifelse(lat>0,"N","S"))
obstrip<-obstrip %>%mutate(season=factor(seasonfunc(month,4)),area=ifelse(lat>0,"N","S"))
head(obstrip)
summary(obstrip)
summary(logtrip)
length(unique(logtrip$area))
length(unique(obstrip$area))

#Fix trips with effort of zero hooks (30 trips, not sure why)
table(logtrip$hooks==0)
obstrip$hooks[obstrip$hooks==0]<-0.01
logtrip$hooks[logtrip$hooks==0]<-0.01

#Summary plot 
x<-obstrip %>% group_by(Year,fleet) %>% summarize(n=n())
g3<-ggplot(x,aes(x=Year,y=n,color=fleet,lty=fleet))+geom_line(lwd=2)+ylab("Number of sets")
grid.arrange(g3, g1,g2, ncol=1)

#Make set by set observer data
obsset<-filter(llsets,trip.05==1)
dim(obsset)
dim(llsets)
logset<-llsets %>%mutate(hooks=hooks/1000,season=factor(seasonfunc(month,4)),area=ifelse(lat>0,"N","S"))
obsset<-obsset %>%mutate(hooks=hooks/1000,season=factor(seasonfunc(month,4)),area=ifelse(lat>0,"N","S"))
table(obsset$hooks==0)
obsset$hooks[obsset$hooks==0]=1/1000

#For including observer bycatch in results, specify sampled and unsampled effort
logtrip$unsampledEffort<-logtrip$hooks
x<-match(obstrip$trip,logtrip$trip)
summary(x)
logtrip$unsampledEffort[x]<-logtrip$unsampledEffort[x]-obstrip$hooks
sum(logtrip$unsampledEffort)
sum(logtrip$unsampledEffort)+sum(obstrip$hooks)
sum(logtrip$hooks)
sum(obstrip$hooks)/sum(logtrip$hooks)
sum(logset$hooks)
logset$unsampledEffort<-logset$hooks*(1-logset$trip.05)

#Make set by set logbook, aggregated by values of X variables, 
logsetAgg<-logset %>% group_by(area,season,fleet,year,w.BUM,hbf) %>%
  summarize(c.BUM=sum(c.BUM),hooks=sum(hooks),sets=n(),
    unsampledEffort=sum(unsampledEffort))
nrow(logsetAgg)
nrow(logsetAgg)/nrow(logset)
sum(logsetAgg$unsampledEffort)
sum(logsetAgg$unsampledEffort)+sum(obsset$hooks)
sum(logset$hooks)

#For validation, make total catch/bycatch tables by year
logyearsum<-llsets%>%group_by(year) %>%
  summarize(Total.BUM=sum(c.BUM),
    Total.SWO=sum(c.SWO)) %>%
  rename(Year=year)
logyearsum

#Print out files 
write.csv(obstrip,"obstrip05.csv")  #Trip by trip observer with 5% effort
write.csv(logtrip,"logtrip05.csv")  #By trip logbook
write.csv(obsset,"obsset05.csv")  #By set observer
write.csv(logset,"logset05.csv")  #By set logbook
write.csv(logsetAgg,"logsetAgg05.csv")  #By set logbook, aggregated by x variables to reduce size
write.csv(logyearsum,"TotalAnnualCatches.csv")

##################OBserver effect bias

## Set up sampling probability proportional to catch of BUM (observer effect bias)
tripSummary<- llsets %>% group_by(trip) %>%
  summarize(BUM=sum(c.BUM),SWO=sum(c.SWO,hooks=sum(hooks),sets=n()))
head(tripSummary)

#Generate probability of being sampled from BUM catch
#p1 is the desired mean coverage and p2 is the minimum
#probability of being sampled, x is the catch
func1<-function(p1,p2,x) {
  m=mean(x)
  M=max(x)
  b=(log(p1/(1-p1))-log(p2/(1-p2)))/(m-M)
  a=log(p1/(1-p1))-b*m  
  plot(seq(0,M,.1),1/(1+exp(-(a+b*seq(0,M,.1)))),xlab="Catch",ylab="Probability of observation",type="l")
  1/(1+exp(-(a+b*x)))
}

# Randomly sample trips based on BUM catch, observer effect bias
tripSummary$sample1<-rbinom(nrow(tripSummary),1,func1(0.05,.01,tripSummary$BUM))
llsets$sample1<-ifelse(llsets$trip%in% tripSummary$sample1[tripSummary$sample1==1],1,0)
#Fraction of trips sampled
sum(tripSummary$sample1)/nrow(tripSummary)
#Fraction of sets sampled
sum(llsets$sample1)/nrow(llsets)

# Check that relationship between catch and P(sampling) exists
glm1<-glm(sample1~BUM,data=tripSummary,family="binomial")
summary(glm1)

## Probability of detection increasing with SWO catch
x<-func1(0.05,.1,tripSummary$SWO)

