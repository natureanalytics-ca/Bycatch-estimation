# This file reads in the observer data and summarizes by trip for further analysis
# For reef longline fishery. 
# !diagnostics off  

library(tidyverse)
library(haven)
library(ggplot2)

#Get Beth's functions
setwd("C:/Users/ebabcock/dropbox/Bycatch Project/Current R code")
source("3.BycatchFunctions.r")

#Uncomment to read in  data
logbook <- read_sas("data/coastal_logbook_15oct19.sas7bdat") #For whole data set
logbooktrip<-read_sas("data/vl_trip_10jul20.sas7bdat")  #This is a by trip summary observed trips only
vlobsmatch<-read.csv("data/VerticalLine_obs_log_match_07oct19.csv",as.is=TRUE)
rfopvl<-read.csv("data/rfop_vl_09oct19.CSV",as.is=TRUE) #FUll observer data
##END reading data

###Use all observer trips from GB and GH categories (check)
obsvl<-rfopvl
dim(rfopvl)
obsvl<-obsvl %>% mutate(ObsType=substring(TRIPNUMBER,1,2))
table(obsvl$ObsType)
obsvl<-obsvl %>% filter(ObsType %in% c("GB","GH"))  %>% droplevels()
length(unique(obsvl$TRIPNUMBER))  #1350 trips

#Extract just reef vertical line data from logbook, and exclude areas outside GOM  (CHeck)
#After converting areas to the old areas for consistency across time series (1-24)
logvl<- logbook %>% mutate(area=areaGOM(AREA)) %>% filter(gear_type %in% c("bandit","hl") & !is.na(area))
length(unique(logvl$LOGBOOK_KEY)) #86623 trips
### End definition of full dataset ####

#Make table with whole to gutted conversion factors from logbook
#only one conversion per species so don't need to match other variables. 
conversion<-logbook %>% filter(!is.na(CONVERSION)&!duplicated(COMMON)) %>% 
  dplyr::select( SPECIES,COMMON,SCIENTIFIC,CONVERSION)
dim(conversion)
summary(conversion)
write.csv(conversion,"Conversion.csv")

# Add totwt_kg to observer dataset
table(is.na(match(obsvl$SCIENTIFIC,conversion$SCIENTIFIC)))
x<-conversion$CONVERSION[match(obsvl$SCIENTIFIC,conversion$SCIENTIFIC)]
x[obsvl$WEIGHT_TYPE!="GUTTED"]=1

#Add whole weight and fix date formats
obsvl<-obsvl %>% mutate(whole_kg=WEIGHT_KG * x,
                        STARTDATE=getdatefunc2(STARTDATE),
                        ENDDATE=getdatefunc2(ENDDATE),
                        PROP_DROPS_SAMPLED=as.numeric(PROP_DROPS_SAMPLED)) 
summary(obsvl)

# Print species summaries, and select a subset of species to analyze (finfish with >10000 kg kept for model testing)
spSummary<-obsvl %>% group_by(COMMON,SCIENTIFIC) %>% 
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
         SCIENTIFIC %in% logvl$SCIENTIFIC) 
commonKeptFinfish$SCIENTIFIC
commonKeptFinfish$COMMON

commonBycatchFinfish<- spSummary %>%
  filter(!grepl("Genus",COMMON),!grepl("Shark",COMMON),!grepl("Dogfish",COMMON),
         !grepl("Skate",COMMON),Ddead.number>100,
         SCIENTIFIC %in% logvl$SCIENTIFIC)   
commonBycatchFinfish$COMMON

Sharks<- spSummary %>%
  filter(grepl("Shark",COMMON),!grepl("Grouped",COMMON), !grepl("Order",COMMON),
    !grepl("Genus",COMMON), !grepl("sucker",COMMON), !(COMMON=="Shark, Dogfish"),
    discard.num.all>100) 
Sharks$COMMON

#Pick species category to run and make sure scientific names match
spSummary<-commonKeptFinfish
#spSummary<-commonBycatchFinfish
x<-match(spSummary$SCIENTIFIC,logvl$SCIENTIFIC)
logvl$SpeciesSelect<-ifelse(logvl$SCIENTIFIC %in% spSummary$SCIENTIFIC,logvl$SCIENTIFIC,"Other")
obsvl$SpeciesSelect<-ifelse(obsvl$SCIENTIFIC %in% spSummary$SCIENTIFIC,obsvl$SCIENTIFIC,"Other")

#Summarize logbook by trip because trip is sample unit for bycatch estimation. 
logtripvl <- logvl %>%
  group_by(LOGBOOK_KEY) %>%  #Check variable names
  summarize(away=median(AWAY,na.rm=TRUE),  
            depth=median(DEPTH,na.rm=TRUE),
            area=as.numeric(mostfreqfunc(area)),
            LANDED=median(LANDED,na.rm=TRUE)) %>%
   mutate(Year=as.numeric(format(LANDED,"%Y")),
           month=as.numeric(format(LANDED,"%m"))) %>%
   mutate(season=ifelse(month<=6,1,2),
           EW=ifelse(area>=11,"W","E"),
           trips=1)
#Note: This version defines 2 seasons, and EW as the area strata
#consistent with the turtle bycatch estimates. Edit if necessary
summary(logtripvl)

# Summarize observer data by trip as sample unit for bycatch estimation
#Effort is days at sea, per Smith. May be divided up based on fraction sampled
x<-obsvl %>% group_by(TRIPNUMBER,SETNUMBER) %>%
  mutate(PROP_DROPS_SAMPLED=ifelse(is.na(PROP_DROPS_SAMPLED),1,PROP_DROPS_SAMPLED),
    REELNO=ifelse(is.na(REELNO),1, REELNO)) %>%
  summarize(proportionreels=mean(PROP_DROPS_SAMPLED[unique(REELNO)],na.rm=TRUE)) %>% 
  group_by(TRIPNUMBER) %>%
   summarize(estsampledsets=sum(proportionreels,na.rm=TRUE))
summary(x)

obsvltrip<-obsvl %>% group_by(TRIPNUMBER) %>% 
  summarize(Year=max(YEAR_LANDED,na.rm=TRUE),
            area=as.numeric(mostfreqfunc(STATZONE)),
            month=max(MONTH_LANDED,na.rm=TRUE),
            enddate=median(ENDDATE,na.rm=TRUE),
            startdate=median(STARTDATE,na.rm=TRUE),
            seadays=median(SEADAYS,na.rm=TRUE),
            totalsets=median(SAMPLED_SETS+UNSAMPLED_SETS,na.rm=TRUE)) %>%
  mutate(season=ifelse(month<=6,"1","2"),
          EW=factor(ifelse(area>=11,"W","E")),
          seadays=ifelse(!is.na(seadays)&seadays!=0,seadays,enddate-startdate+1))
obsvltrip<-merge(obsvltrip,x)
summary(obsvltrip)
#There are 13 trips with no data on area and no sets sampled. Removed them.
obsvltrip<-dplyr::filter(obsvltrip,totalsets!=0) %>%  
   mutate(decDaysSampled=(estsampledsets/totalsets)*seadays,
        decDaysNotSampled=(1-estsampledsets/totalsets)*seadays,
        SeaDaysNotSampled=0)
dim(obsvltrip)
summary(obsvltrip)

#Add selected catch type as column to trip summary
catchtype<-"kept.kg"
#catchtype<-"discard.dead.num"
if(catchtype=="kept.kg") {
  catch<-obsvl %>% group_by(TRIPNUMBER,SpeciesSelect) %>% 
  summarize(catch=sum(whole_kg[FATE %in% c("KEPT","KEPT AS BAIT")],na.rm=TRUE)) 
}
if(catchtype=="discard.dead.num") {
  catch<-obsvl %>% group_by(TRIPNUMBER,SpeciesSelect) %>% 
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
obsvltrip<-merge(obsvltrip,catch)
#summary(obsvltrip$propdropssampled) #Check that none have zero catch

#Summarize logbook data by year 
logyearsum<-logtripvl %>% group_by(Year) %>% 
  summarize(effort=sum(away,na.rm=TRUE))
logyearsum
spyearsum<-logvl %>% 
  mutate(Year=as.numeric(format(LANDED,"%Y"))) %>% 
  group_by(Year,SpeciesSelect) %>%
  summarize(Kept.kg=0.453592*sum(tot_wholelbs,na.rm=TRUE)) 
spyearsum<-spyearsum%>% pivot_wider(names_from = SpeciesSelect,values_from = Kept.kg)
spyearsum[is.na(spyearsum)]<-0
names(spyearsum)[-1]<-paste0("catch.",names(spyearsum)[-1])  
logyearsum<-merge(logyearsum,spyearsum)
logyearsum

summary(obsvltrip)
