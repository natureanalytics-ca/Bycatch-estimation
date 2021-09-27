#This sets global conditions and makes a preliminary data summary
#Load required libraries
library(tidyverse)
library(ggplot2)
library(MASS)
library(lme4)
library(cplm)
library(DHARMa)
library(quantreg)
library(tidyselect)
library(MuMIn)
library(gridExtra)
library(foreach)
library(doParallel)
library(pdftools)
library(tweedie)
library(reshape2)
library(glmmTMB)

source("3.BycatchFunctions.r")

#Set global conditions
theme_set(theme_bw()) #ggplot theme
defaultW <- 0
options(warn=defaultW)
nSims<-1000  #For simulating CIs where needed
NumCores<-detectCores()  #Check if machine has multiple cores for parallel processing
#Make sure binomial is included if either of the delta models is
if(("Delta-Lognormal" %in% modelTry |"Delta-Gamma" %in% modelTry) & !"Binomial" %in% modelTry)
  modelTry<-c("Binomial",modelTry)

# Set up variables
obsdat<-obsdat %>%   ungroup() %>%
  rename(Year=!!yearVar)
logdat<-logdat %>%  ungroup() %>%
  rename(Year=!!yearVar)
requiredVarNames<-as.vector(getAllTerms(simpleModel))
allVarNames<-as.vector(getAllTerms(complexModel))  
allVarNames<-allVarNames[grep(":",allVarNames,invert=TRUE)]
allVarNames<-allVarNames[grep("I(*)",allVarNames,invert=TRUE)]
if(!all(allVarNames %in% names(obsdat))) 
  print(paste0("Variable ", allVarNames[!allVarNames%in% names(obsdat) ], " not found in observer data"))
if(!all(allVarNames %in% names(logdat))) 
  print(paste0("Variable ", allVarNames[!allVarNames%in% names(logdat) ], " not found in logbook data"))
#It's all right not to see the variable name if it is a function of another variable that is present
indexVarNames<-as.vector(getAllTerms(indexModel))
if(!"Year" %in% indexVarNames) indexVarNames<-c("Year",indexVarNames)

#Set up data frames
obsdat<-obsdat %>%  
  rename(Effort=!!obsEffort) %>%
  mutate_at(vars(all_of(factorNames)),factor) 
if(EstimateBycatch) {
 if(is.na(logNum))   { 
   logdat<-mutate(logdat,SampleUnits=1)
   logNum="SampleUnits"
 }
 logdat<-logdat %>%  
   rename(Effort=!!logEffort,SampleUnits=!!logNum) %>%
  mutate_at(vars(all_of(factorNames)),factor) 
 if(logEffort==sampleUnit) logdat<-mutate(logdat,Effort=SampleUnits)
 if(includeObsCatch) {
  obsdat<-obsdat %>% rename(matchColumn=!!matchColumn)
  logdat<-logdat %>% rename(matchColumn=!!matchColumn,unsampledEffort=!!logUnsampledEffort)
 }
}

# See if data set is large
BigData<-ifelse(sum(logdat$SampleUnits)>10000,TRUE,FALSE)
#Add stratum designation and check sample size in strata
if(length(requiredVarNames)>1) {
 logdat$strata<-apply( logdat[ , requiredVarNames ] , 1 , paste , collapse = "-" ) 
} else {
 logdat$strata <- pull(logdat,var=requiredVarNames)
}
if(max(tapply(logdat$SampleUnits,logdat$strata,sum))>20000) {
  print("Cannot calculate variance for large number of logbook sample units")
  VarCalc<-"None"
}  

#newDat for making index
newDat<-distinct_at(obsdat,vars(all_of(indexVarNames)),.keep_all=TRUE) %>%
 arrange(Year) %>%
 mutate(Effort=1)
temp<-allVarNames[allVarNames != "Year"]
for(i in 1:length(temp)) {
  if(!temp[i] %in% indexVarNames) {
   if(is.numeric(pull(obsdat,!!temp[i])))
    newDat[,temp[i]]<-median(pull(obsdat,!!temp[i]),na.rm=TRUE) else
    newDat[,temp[i]]<-mostfreqfunc(obsdat[,temp[i]]) 
   } 
  }

#Subtract first year if numeric to improve convergence
if(is.numeric(obsdat$Year)) {
 startYear<-min(obsdat$Year)
 obsdat$Year<-obsdat$Year-startYear
 logdat$Year<-logdat$Year-startYear
 newDat$Year<-newDat$Year-startYear
} else startYear<-min(as.numeric(as.character(obsdat$Year)))

#Set up directory for output
setwd(baseDir)
numSp<-length(sp)
if(!dir.exists(paste0("Output",runName))) dir.create(paste0("Output",runName))
outDir<-paste0(baseDir,"/output",runName)

#Make R objects to store analysis
bestmod<-NULL
predbestmod<-list()
indexbestmod<-list()
allmods<-list()
allindex<-list()
modelTable<-list()
modelSelectTable<-list()
modFits<-list()
modPredVals<-list()
modIndexVals<-list()
modelFail<-matrix("-",numSp,length(modelTry),dimnames=list(common,modelTry))
rmsetab<-list()
metab<-list()
residualTab<-list()

#Make lists to keep output, which will also be output as .pdf and .csv files for use in reports.
dirname<-list()
dat<-list()
#Loop through all species and print data summary. Note that records with NA in either catch or effort are excluded automatically
yearSum<-list()
if(NumCores>3 & numSp>1)  {
  cl<-makeCluster(NumCores-2)
  registerDoParallel(cl)
}
foreach(run= 1:numSp) %do%  {
  dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/")
  if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]])
  if(includeObsCatch) tempvars<-c(allVarNames,"Effort","Catch","matchColumn") else
    tempvars<-c(allVarNames,"Effort","Catch")
  dat[[run]]<-obsdat %>%
    rename(Catch=!!obsCatch[run])%>%
    dplyr::select_at(all_of(tempvars)) %>%
    drop_na()   %>%
    mutate(cpue=Catch/Effort,
           log.cpue=log(Catch/Effort),
           pres=ifelse(cpue>0,1,0)) 
  if(dim(dat[[run]])[1]<dim(obsdat)[1]) print(paste0("Removed ",dim(obsdat)[1]-dim(dat[[run]])[1]," rows with NA values for ",common[run]))
  yearSum[[run]]<-dat[[run]] %>% group_by(Year) %>%
    summarize(ObsCat=sum(Catch,na.rm=TRUE),
              ObsEff=sum(Effort,na.rm=TRUE),
              ObsUnits=length(Year),
              CPUE=mean(cpue,na.rm=TRUE),
              CPUEse=standard.error(cpue),
              Outlr=outlierCountFunc(cpue),
              Pos=sum(pres,na.rm=TRUE)) %>%
    mutate(PosFrac=Pos/ObsUnits)
  if(EstimateBycatch) {
   x<-logdat  %>% group_by(Year) %>%
    summarize(Effort=sum(Effort,na.rm=TRUE),Units=sum(SampleUnits)) 
   yearSum[[run]]<-merge(yearSum[[run]],x) %>% mutate(EffObsFrac=ObsEff/Effort,
                                                     UnitsObsFrac=ObsUnits/Units)
   logyear<-logdat %>% group_by(Year) %>% summarize(Effort=sum(Effort,na.rm=TRUE))
   x=ratio.func(dat[[run]]$Effort,dat[[run]]$Catch,dat[[run]]$Year,
               logyear$Effort,logyear$Effort,logyear$Year)
   yearSum[[run]]<-cbind(yearSum[[run]],CatEst=x$stratum.est,Catse=x$stratum.se) %>% 
    ungroup() %>% mutate(Year=as.numeric(as.character(Year))) %>%
     mutate(Year=ifelse(Year<startYear,Year+startYear,Year)) %>%
    dplyr::rename(!!paste0("Obs",sampleUnit):=ObsUnits,                                 ,
                  !!sampleUnit:=Units,
                  !!paste0(sampleUnit,"ObsFrac"):=UnitsObsFrac)
  } 
  write.csv(yearSum[[run]],paste0(dirname[[run]],common[run],catchType[run],"DataSummary.csv"))
}
if(NumCores>3 & numSp>1) stopCluster(cl)
save(list=c("numSp","yearSum","runName", "common", "sp"),file=paste0(outDir,"\\","sumdatR"))
rmarkdown::render("5.PrintDataSummary.rmd",
  params=list(OutDir=OutDir),
  output_file = "DataSummary.pdf",
  output_dir=outDir)
file.copy(specFile,paste0(outDir,"\\",Sys.Date(),"BycatchModelSpecification.r"))


for(run in 1:numSp) {
 residualTab[[run]]<-matrix(0,8,length(modelTry),dimnames=list(c("KS.D","KS.p", "Dispersion.ratio","Dispersion.p" , 
                            "ZeroInf.ratio" ,"ZeroInf.p","Outlier" , "Outlier.p"),
                            modelTry))
 modelTable[[run]]<-data.frame(model=modelTry,
     formula=rep("",length(modelTry)),
     RMSE=rep(NA,length(modelTry)),
     ME=rep(NA,length(modelTry)))
 modPredVals[[run]]<-rep(list(NULL),length(modelTry))
 names(modPredVals[[run]])<-modelTry
 modIndexVals[[run]]<- modPredVals[[run]]
 modFits[[run]]<- modPredVals[[run]]
 modelSelectTable[[run]]<- modPredVals[[run]]
}

