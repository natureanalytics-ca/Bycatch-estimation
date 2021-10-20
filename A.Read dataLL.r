# This file reads in the observer data and summarizes by trip for further analysis
# For reef longline fishery. 
# !diagnostics off  
library(tidyverse)
library(haven)
library(ggplot2)

#Set the working director to the correct location
setwd("C:/Users/ebabcock/Box Sync/bycatch project (ebabcock@miami.edu)/Current R code")
#Get Beth's functions
source("3.BycatchFunctions.r")

#Uncomment to read in  data
 logbook <- read_sas("data/coastal_logbook_15oct19.sas7bdat")
 trip<-read_sas("data/ll_trip_10jul20.sas7bdat")  #This is a by trip by trip summary of observed effort
 llobsmatch<-read.csv("data/longline_obs_log_match_04oct19.csv",as.is=TRUE)
 rfopll<-read.csv("data/rfop_ll_09oct19.CSV",as.is=TRUE)
##END reading data

###Use all observer trips that have key "GL" in trip number. 
# This is regular longline trips in Gulf of Mexico
obsll<-rfopll
dim(rfopll)
obsll<-obsll %>% mutate(ObsType=substring(TRIPNUMBER,1,2))
table(obsll$ObsType)
obsll<-obsll %>% filter(ObsType=="GL")  %>% droplevels()
length(unique(obsll$TRIPNUMBER))  #423 trips

#Extract just reef Longline data from logbook, and exclude areas outside GOM  
#After converting areas to the old areas for consistency across time series (1-21)
logll<- logbook %>% mutate(area=areaGOM(AREA))%>% 
  filter(gear_type=="ll_reeffish" & !is.na(area))
length(unique(logll$LOGBOOK_KEY)) #11299 trips
### End definition of full dataset ####

#Make table with whole to gutted conversion factors from logbook
conversion<-logbook %>% filter(!is.na(CONVERSION)&!duplicated(COMMON)) %>% 
  dplyr::select( SPECIES,COMMON,SCIENTIFIC,CONVERSION)
dim(conversion)
summary(conversion)
write.csv(conversion,"Conversion.csv")

# Add totwt_kg to observer dataset
x<-conversion$CONVERSION[match(obsll$SCIENTIFIC,conversion$SCIENTIFIC)]
x[obsll$WEIGHT_TYPE!="GUTTED"]=1
obsll<-obsll %>% mutate(whole_kg=WEIGHT_KG * x)
summary(obsll)

# Print species summaries, and select a subset of species to analyze 
spSummary<-obsll %>% group_by(COMMON,SCIENTIFIC) %>% 
  summarize(Kept.kg=sum(whole_kg[FATE=="KEPT" | FATE=="KEPT AS BAIT"],na.rm=TRUE),
            Kept.number=length(whole_kg[(FATE=="KEPT"| FATE=="KEPT AS BAIT" )& !is.na(SCIENTIFIC)]),
            Dalive.kg=sum(whole_kg[FATE=="RELEASED ALIVE"],na.rm=TRUE),
            Dalive.number=length(whole_kg[FATE=="RELEASED ALIVE" & !is.na(SCIENTIFIC)]),
            Ddead.kg=sum(whole_kg[FATE=="RELEASED DEAD" | FATE=="RELEASED UNKNOWN"], na.rm=TRUE),
            Ddead.number=length(whole_kg[(FATE=="RELEASED DEAD" | FATE=="RELEASED UNKNOWN" )& !is.na(SCIENTIFIC)]),
            discard.num.all=length(whole_kg[is.na(whole_kg) & DISPOSITION=="DISCARD"]),
            kept.num.all=length(whole_kg[is.na(whole_kg) & DISPOSITION=="KEPT"])) 

# Pull out species groups for testing
commonKeptFinfish<- spSummary %>%
  filter(!grepl("Genus",COMMON),!grepl("Shark",COMMON),!grepl("Dogfish",COMMON),
         !grepl("Skate",COMMON),Kept.kg>10000,
         SCIENTIFIC %in% logll$SCIENTIFIC) 
commonKeptFinfish$COMMON
commonBycatchFinfish<- spSummary %>%
  filter(!grepl("Genus",COMMON),!grepl("Shark",COMMON),!grepl("Dogfish",COMMON),
         !grepl("Skate",COMMON),Ddead.number>100,
         SCIENTIFIC %in% logll$SCIENTIFIC)   
commonBycatchFinfish$COMMON
Sharks<- spSummary %>%
  filter(grepl("Shark",COMMON),!grepl("Grouped",COMMON), !grepl("Order",COMMON),
    !grepl("Genus",COMMON), !grepl("sucker",COMMON), !(COMMON=="Shark, Dogfish"),
    discard.num.all>100) 
Sharks$COMMON

#Pick species category to run
spSummary<-Sharks
#spSummary<-commonKeptFinfish
#spSummary<-commonBycatchFinfish
logll$SpeciesSelect<-ifelse(logll$SCIENTIFIC %in% spSummary$SCIENTIFIC,logll$SCIENTIFIC,"Other")
obsll$SpeciesSelect<-ifelse(obsll$SCIENTIFIC %in% spSummary$SCIENTIFIC,obsll$SCIENTIFIC,"Other")

#Summarize logbook by trip because trip is sample unit for bycatch estimation. 
logtripll <- logll %>% 
  group_by(LOGBOOK_KEY) %>%  #Check variable names
  summarize(sets=median(numgear,na.rm=TRUE),
            UNLOAD=median(UNLOAD,na.rm=TRUE),
            depth=median(DEPTH,na.rm=TRUE)*0.3048,  #convert to meters
            area=as.numeric(mostfreqfunc(area))) %>%
   mutate(Year=as.numeric(format(UNLOAD,"%Y")),
           month=as.numeric(format(UNLOAD,"%m")),
           depth=ifelse(!is.na(depth),depth,median(depth,na.rm=TRUE))) %>%
   mutate(season=ifelse(month<=6,1,2),EW=ifelse(area>=11,"W","E"),trips=1)
#Note: This version defines 2 seasons, and EW as the area strata
#consistent with the turtle bycatch estimates. Edit if necessary
summary(logtripll)


# Summarize observer data by trip as sample unit for bycatch estimation
obslltrip<-obsll %>% group_by(TRIPNUMBER) %>% 
  summarize(Year=max(YEAR_LANDED,na.rm=TRUE),
            area=as.numeric(mostfreqfunc(STATZONE)), 
            month=max(MONTH_LANDED,na.rm=TRUE),
            sampled.sets=median(SAMPLED_SETS,na.rm=TRUE),
            unsampled.sets=median(UNSAMPLED_SETS,na.rm=TRUE)) %>%
          mutate(season=ifelse(month<=6,"1","2"), #Defined season as in turtle refs
          EW=ifelse(area>=11,"W","E"))   #Defined East vs. West
dim(obslltrip)
summary(obslltrip)

#Add selected catch type as column to trip summary
#catchtype<-"kept.kg"
#catchtype<-"discard.dead.num"
catchtype<-"all.num"

if(catchtype=="kept.kg") {
  catch<-obsll %>% group_by(TRIPNUMBER,SpeciesSelect) %>% 
  summarize(catch=sum(whole_kg[FATE %in% c("KEPT","KEPT AS BAIT")],na.rm=TRUE)) 
}
if(catchtype=="discard.dead.num") {
  catch<-obsll %>% group_by(TRIPNUMBER,SpeciesSelect) %>% 
  summarize(catch=length(whole_kg[FATE %in% c("RELEASED DEAD","RELEASED UNKNOWN")])) 
}
if(catchtype=="all.num") {
  catch<-obsll %>% group_by(TRIPNUMBER,SpeciesSelect) %>% 
    summarize(catch=length(whole_kg)) 
}

# Add species catch
catch<-catch %>% pivot_wider(names_from = SpeciesSelect,values_from = catch)
catch[is.na(catch)]<-0
names(catch)[-1] <-paste0("catch.",names(catch)[-1])
obslltrip<-merge(obslltrip,catch)
names(obslltrip)


#Summarize logbook data by year and stratum
logyearsum<-logtripll %>% group_by(Year) %>% 
  summarize(effort=sum(sets,na.rm=TRUE))
logyearsum
spyearsum<-logll %>% 
  mutate(Year=as.numeric(format(UNLOAD,"%Y"))) %>% 
  group_by(Year,SpeciesSelect) %>%
  summarize(Kept.kg=0.453592*sum(tot_wholelbs,na.rm=TRUE)) 
spyearsum<-spyearsum%>% pivot_wider(names_from = SpeciesSelect,values_from = Kept.kg)
spyearsum[is.na(spyearsum)]<-0
names(spyearsum)[-1]<-paste0("catch.",names(spyearsum)[-1])  
logyearsum<-merge(logyearsum,spyearsum)
logyearsum
write.csv(logyearsum,"ReefLLkept.csv")

#Summary of log data by stratum
logyearstratsum<-logtripll %>% 
  group_by(Year,EW,season) %>% 
  summarize(trips=length(Year),sets=sum(sets,na.rm=TRUE))
logyearstratsum


# Make data by set to include depth as a variable

obsllset<-obsll %>% group_by(TRIPNUMBER,SETNUMBER) %>% 
  summarize(Year=max(YEAR_LANDED,na.rm=TRUE),
            area=as.numeric(mostfreqfunc(STATZONE)),
            month=max(MONTH_LANDED,na.rm=TRUE),
            depth=mean(FISHINGDEPTH_M)) %>%
  mutate(season=ifelse(month<=6,"1","2"),
         EW=ifelse(area>=11,"W","E"),
         depth=ifelse(is.na(depth),median(depth,na.rm=TRUE),depth),
         sets=rep(1,length(Year))) 
dim(obsllset)
summary(obsllset)
#Add selected catch type as column to set summary (should be same as above)
catchtype<-"kept.kg"
#catchtype<-"discard.dead.num"
#catchtype<-"all.encountered.num"
if(catchtype=="kept.kg") {
  catch<-obsll %>% group_by(TRIPNUMBER,SETNUMBER,SpeciesSelect) %>% 
    summarize(catch=sum(whole_kg[FATE %in% c("KEPT","KEPT AS BAIT")],na.rm=TRUE)) 
}
if(catchtype=="discard.dead.num") {
  catch<-obsll %>% group_by(TRIPNUMBER,SETNUMBER,SpeciesSelect) %>% 
    summarize(catch=length(whole_kg[FATE %in% c("RELEASED DEAD","RELEASED UNKNOWN")])) 
}
if(catchtype=="all.encountered.num") {
  catch<-obsll %>% group_by(TRIPNUMBER,SETNUMBER,SpeciesSelect) %>% 
    summarize(catch=length(whole_kg)) 
}
catch<-catch %>% pivot_wider(names_from = SpeciesSelect,values_from = catch)
catch[is.na(catch)]<-0
names(catch)[c(-1,-2)] <-paste0("catch.",names(catch)[c(-1,-2)])
obsllset<-merge(obsllset,catch)
##
dim(obslltrip)
summary(obsllset)
dim(obsllset)
sum(obslltrip$sampled.sets)
names(obsllset)


# Add column to match logbook to observer data
a<-match(obslltrip$TRIPNUMBER,llobsmatch$TripNumber)
table(obslltrip$Year,is.na(a))
obslltrip$LOGBOOK_KEY<-llobsmatch$LOGBOOK_KEY[a]
b<-match(logtripll$LOGBOOK_KEY,llobsmatch$LOGBOOK_KEY)
table(logtripll$Year,is.na(b))
logtripll$TRIPNUMBER<-llobsmatch$TripNumber[b]
summary(obslltrip$LOGBOOK_KEY)
obslltripmatch=obslltrip[!is.na(obslltrip$LOGBOOK_KEY),]
dim(obslltrip)
dim(obslltripmatch)

logtripll$SampledEffort<-rep(NA,nrow(logtripll))
for(i in 1:nrow(obslltripmatch)) {
  logtripll$SampledEffort[logtripll$LOGBOOK_KEY==obslltripmatch$LOGBOOK_KEY[i]]<-obslltripmatch$sampled.sets[i]
}
ggplot(filter(logtripll,!is.na(SampledEffort)),aes(y=SampledEffort,x=sets))+
  geom_point()+geom_abline(intercept=0,slope=1) + xlab("Logbook effort")+ ylab("Observed effort")+
  ggtitle("Observed and logbook number of hooks, same trips")
table(logtripll$SampledEffort>logtripll$sets)
#87 matched trips have more observer effort than logbook effort. 



