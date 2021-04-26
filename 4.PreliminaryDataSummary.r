#This sets global conditions and makes a preliminary data summary

#Set global conditions
theme_set(theme_bw())
defaultW <- 0
options(warn=defaultW)
#options(warn = -1) 
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = .8)),
  colhead = list(fg_params=list(cex = .8)),
  rowhead = list(fg_params=list(cex = .8)))
nSims=1000  #For simulating CIs where needed

# Set up variables
requiredVarNames<-as.vector(getAllTerms(simpleModel))
allVarNames<-as.vector(getAllTerms(complexModel))  
allVarNames<-allVarNames[-grep(":",allVarNames)]
if(!all(allVarNames %in% names(obsdat))) 
  print(paste0("Variable ", allVarNames[!allVarNames%in% names(obsdat) ], " not found in observer data"))
if(!all(allVarNames %in% names(logdat))) 
  print(paste0("Variable ", allVarNames[!allVarNames%in% names(logdat) ], " not found in logbook data"))
#It's all right not to see the variable name if it is a function of another variable that is present

#Set up data frames
#UnObsEffort=!!obsEffortNotSampled
obsdat<-obsdat %>%   ungroup() %>%
  rename(Effort=!!obsEffort, Year=!!yearVar) %>%
  mutate_at(vars(all_of(factorNames)),factor) 
obsdat$stratum=interaction(obsdat[,allVarNames])
logdat<-logdat %>%  ungroup() %>%
  rename(Effort=!!logEffort,Year=!!yearVar,SampleUnits=!!logNum) %>%
  mutate_at(vars(all_of(factorNames)),factor) 
if(logEffort==sampleUnit) logdat<-mutate(logdat,Effort=SampleUnits)
logdat$stratum=interaction(logdat[,allVarNames])

#if(!all(obsdat$stratum %in% logdat$stratum)) print("Not all observed trips match")

#Set up directory for output
setwd(baseDir)
source("3.BycatchFunctions.r")
numSp<-length(sp)
if(!dir.exists(paste0("Output",runName))) dir.create(paste0("Output",runName))
outDir<-paste0(baseDir,"/output",runName)

#Make lists to keep output, which will also be output as .pdf and .csv files for use in reports.
dirname<-list()
dat<-list()
#Loop through all species and print data summary. Note that records with NA in either catch or effort are excluded automatically
yearSum<-list()
fileList<-NULL
for(run in 1:numSp) {
  dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/")
  if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]])
  dat[[run]]<-obsdat %>%
    rename(Catch=!!obsCatch[run])%>%
    dplyr::select_at(all_of(c(allVarNames,"stratum","Effort","Catch"))) %>%
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
  x<-logdat  %>% group_by(Year) %>%
    summarize(Effort=sum(Effort,na.rm=TRUE),Units=sum(SampleUnits)) 
  yearSum[[run]]<-merge(yearSum[[run]],x) %>% mutate(EffObsFrac=ObsEff/Effort,
                                                     UnitsObsFrac=ObsUnits/Units)
  logyear<-logdat %>% group_by(Year) %>% summarize(Effort=sum(Effort,na.rm=TRUE))
  x=ratio.func(dat[[run]]$Effort,dat[[run]]$Catch,dat[[run]]$Year,
               logyear$Effort,logyear$Effort,logyear$Year)
  yearSum[[run]]<-cbind(yearSum[[run]],CatEst=x$stratum.est,Catse=x$stratum.se) %>% 
    ungroup() %>% mutate(Year=as.numeric(as.character(Year))) %>%
    dplyr::rename(!!paste0("Obs",sampleUnit):=ObsUnits,                                 ,
                  !!sampleUnit:=Units,
                  !!paste0(sampleUnit,"ObsFrac"):=UnitsObsFrac)
  write.csv(yearSum[[run]],paste0(dirname[[run]],common[run],catchType[run],"DataSummary.csv"))
  printTableFunc("Data summary",sp[run],yearSum[[run]],paste0(dirname[[run]],common[run],catchType[run],"DataSummary.pdf"))
  fileList<-c(fileList,paste0(dirname[[run]],common[run],catchType[run],"DataSummary.pdf"))
}
pdf_combine(fileList,paste0(outDir,"\\DataSummary.pdf"))
file.copy(specFile,paste0(outDir,"\\",Sys.Date(),"BycatchModelSpecification.r"))

