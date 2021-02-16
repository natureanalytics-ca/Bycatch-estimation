# This code runs a generic model-based bycatch estimation procedure. 
# The observer data should be aggregated to the appropriate sample unit, either trips or sets.
# Effort must be in the same units in both data sets (e.g. hook hours), but it does not matter
# how the logbook effort data is aggregated, as long has it includes data on all stratification or 
# predictor variables. FOr example, the data can be aggregated by year, East-West and season if those
# are the stratification variables, or it can be aggregated by trip.
# The user must specify the names of the databases and the variables to include in the 
# The .r file named on line 12. Everything else should work automatically.
# Contact Beth Babcock ebabcock@rsmas.miami.edu for assistance. 

############### Step 1. Enter the data specification in the file named here #############################
specFile<-"C:/Users/ebabcock/Dropbox/bycatch project/Current R code/1.BycatchModelSpecification.r"
# Either set the working directory or put the full patch in the filename.
# Complete the information in the file before continuing. You may run through specFile line by line, but it 
# will also be sourced again later. The file will be saved, with the addition of the date, to the output directory

##################################################################################
### From here on no changes should be needed.  However, there are places where you should stop and look at 
#  the output before proceeding. 

#################################################################################
#Load required libraries
library(tidyverse)
library(ggplot2)
library(MASS)
library(lme4)
library(cplm)
library(DHARMa)
library(tidyselect)
library(MuMIn)
library(gridExtra)
library(grid)
library(gtable)
library(gt)
library(pdftools)
library(tweedie)
library(reshape2)
library(glmmTMB)

########################################################################
#Read in data specification
source(specFile)
#Make data summary outputs
source(paste0(baseDir,"/4.preliminaryDataSummary.r"))
# Stop here and check data summaries to make sure each species has reasonable number
# of observations in each year for analysis.
#The file called  "Data summary all species.pdf" has all the tables. 
#There is also a .csv for each species. 

############## This is the main analysis loop ##########################
#The following runs all the models and prints all the output files. 
bestmod<-NULL
predbestmod<-list()
allmods<-list()
modelTable<-list()
modelSelectTable<-list()
modelFail<-matrix("-",numSp,length(modelTry)+1,dimnames=list(common,c("Binomial",modelTry)))
rmsetab<-list()
metab<-list()
outFiles<-list()
residualTab<-list()
for(run in 1:numSp) {
 datval<-dat[[run]]
 residualTab[[run]]<-matrix(0,8,length(modelTry)+1,dimnames=list(c("KS.D","KS.p", "Dispersion.ratio","Dispersion.p" , 
                            "ZeroInf.ratio" ,"ZeroInf.p","Outlier" , "Outlier.p"),
                            c("Binomial",modelTry)))
 modelTable[[run]]<-data.frame(model=c("Binomial",modelTry),
     formula=rep("",length(modelTry)+1),
     RMSE=rep(NA,length(modelTry)+1),
   ME=rep(NA,length(modelTry)+1))
 outVal<-paste0(dirname[[run]],common[run],catchType[run])
 outFiles[[run]]<-c(fileList[[run]],
   paste0(outVal,rep(c("Fit","Residuals","ModelSelection"),length(modelTry)+1),
     rep(dimnames(modelFail)[[2]],each=3),".pdf"),
   paste0(outVal,"AllFit.pdf"),paste0(outVal,"modelSummary.pdf"))
 if(DoCrossValidation) outFiles[[run]]<-c(outFiles[[run]],paste0(outVal,"Crossvalidation.pdf"),paste0(outVal,"BestFit.pdf"))
#Find best binomial model and print outputs. 
 varExclude<-NULL
 bin1<-findBestModelFunc(datval,"Binomial",printOutput=TRUE)
 if(!is.null(bin1)) {  #If model converged make predictions and print figures
  bintab<-bin1[[2]]
  bin1<-bin1[[1]]
  binpredvals<-makePredictionsVar(bin1,modType="Binomial",newdat=logdat,printOutput=TRUE)   
  plotFits(binpredvals,"Binomial",paste0(outVal,"FitBinomial.pdf"))
  residualTab[[run]][,"Binomial"]<-ResidualsFunc(bin1,"Binomial",paste0(outVal,"ResidualsBinomial.pdf"))
  if(residualTab[[run]]["KS.p","Binomial"]<0.05) modelFail[run,"Binomial"]<-"resid"
  if(is.null(binpredvals)) modelFail[run,"Binomial"]<-"cv"
  modelTable[[run]]$formula[1]<-paste(formula(bin1))[[3]]
  } else  {
    modelFail[run,"Binomial"]<-"fit"
    binpredvals<-NULL
  }
 predvals=rep(list(NULL),length(modelTry))
 names(predvals)=modelTry
 modfits=rep(list(NULL),length(modelTry))
 names(modfits)=modelTry
 if("Lognormal" %in% modelTry | "Gamma" %in% modelTry) {  #Delta models if requested
   posdat<-filter(dat[[run]],pres==1)
   y<-unlist(lapply(posdat[,factorNames],function(x) length(setdiff(levels(x),x)))) #See if all levels are included
   varExclude<-names(y)[y>0]
   if(length(varExclude>0)) print(paste(common[run], "excluding variable",varExclude,"from delta models for positive catch"))
   if(min(summary(posdat$Year))>0 &!is.null(bin1)) { #If all years have at least one positive observation and binomial converged, carry on with delta models
     for(mod in which(modelTry %in% c("Lognormal","Gamma"))) {
       modfits[[modelTry[mod]]]<-findBestModelFunc(posdat,modelTry[mod],printOutput=TRUE)
       if(!is.null(modfits[[modelTry[mod]]])) {
         modelSelectTable[[run]]<-modfits[[modelTry[mod]]][[2]]
         modfits[[modelTry[mod]]]<-modfits[[modelTry[mod]]][[1]]
         predvals[[modelTry[mod]]]<-makePredictionsVar(modfit1=bin1,modfit2=modfits[[modelTry[mod]]],modType=modelTry[mod],newdat=logdat,printOutput=TRUE)   
         plotFits(predvals[[modelTry[mod]]],modelTry[mod],paste0(outVal,"Fit",modelTry[mod],".pdf"))
         residualTab[[run]][,modelTry[mod]]<-ResidualsFunc(modfits[[modelTry[mod]]],modelTry[mod],paste0(outVal,"Residuals",modelTry[mod],".pdf"))
         if(residualTab[[run]]["KS.p",modelTry[mod]]<0.05) modelFail[run,modelTry[mod]]<-"resid"
         if(is.null(predvals[[modelTry[mod]]])) modelFail[run,modelTry[mod]]<-"cv"
         modelTable[[run]]$formula[mod+1]<-paste(formula(modfits[[modelTry[mod]]]))[[3]]
       } else {
         modelFail[run,modelTry[mod]]<-"fit"
       }
     }
   } else {
     print("Not all years have positive observations, skipping delta models")
     modelFail[run,c("Lognormal","Gamma")]<-"data"
   }}
 #All models except delta
 for(mod in which(!modelTry %in% c("Lognormal","Gamma"))) {
   modfits[[modelTry[mod]]]<-findBestModelFunc(datval,modelTry[mod],printOutput=TRUE)
   if(!is.null(modfits[[modelTry[mod]]])) {
     modelSelectTable[[run]]<-modfits[[modelTry[mod]]][[2]]
     modfits[[modelTry[mod]]]<-modfits[[modelTry[mod]]][[1]]
     predvals[[modelTry[mod]]]<-makePredictionsVar(modfit1=modfits[[modelTry[mod]]],modType=modelTry[mod],newdat=logdat,printOutput=TRUE)   
     plotFits(predvals[[modelTry[mod]]],modelTry[mod],paste0(outVal,"Fit",modelTry[mod],".pdf"))
     modelTable[[run]]$formula[mod+1]<-paste(formula(modfits[[modelTry[mod]]]))[[3]]
     temp<-ResidualsFunc(modfits[[modelTry[mod]]],modelTry[mod],paste0(outVal,"Residuals",modelTry[mod],".pdf"))
     if(!is.null(temp)) {
       residualTab[[run]][,modelTry[mod]]<-temp
       if(residualTab[[run]]["KS.p",modelTry[mod]]<0.05) modelFail[run,modelTry[mod]]<-"resid"
     }
     if(is.null(predvals[[modelTry[mod]]])) modelFail[run,modelTry[mod]]<-"cv"
   } else {
     modelFail[run,modelTry[mod]]<-"fit"
   }
 }
  #Compare all predictions 
  yearsumgraph<-yearSum[[run]] %>% dplyr::select(Year=Year,Total=CatEst,Total.se=Catse) %>%
    mutate(TotalVar=Total.se^2,Total.cv=Total.se/Total,TotalFixed=Total)
  allmods[[run]]<-bind_rows(c(predvals,list(Ratio=yearsumgraph)),.id="Source") 
  plotFits(allmods[[run]],modType="All",paste0(outVal,"AllFit.pdf"))
  #Show the diagnostic table
  printTableFunc("Diagnostics",sp[run],residualTab[[run]],paste0(outVal,"residualDiagnostics.pdf"),useRowNames = TRUE)
  write.csv(residualTab[[run]],paste0(outVal,"residualDiagnostics.csv"))
  ######## Cross validation 10 fold  ####################################
 if(DoCrossValidation & length(which(modelFail[run,-1]=="-"))>0) {  #Don't do unless at least one model worked
  datval$cvsample<-sample(rep(1:10,length=dim(datval)[1]),replace=FALSE)
  table(datval$cvsample,datval$Year)
  rmsetab[[run]]<-data.frame(matrix(NA,10,length(modelTry),dimnames=list(1:10,modelTry)))
  metab[[run]]<-rmsetab[[run]]
  for(i in 1:10 ) {
   datin<-datval[datval$cvsample!=i,]
   datout<-datval[datval$cvsample==i,]
   datout$SampleUnits<-rep(1,dim(datout)[1])
   bin1<-findBestModelFunc(datin,"Binomial")[[1]]
   if("Lognormal" %in% modelTry | "Gamma" %in% modelTry) { 
     posdat<-filter(dat[[run]],pres==1)
     for(mod in which(modelTry %in% c("Lognormal","Gamma"))) {
       if(modelFail[run,modelTry[mod]]=="-" & min(summary(posdat$Year))>0) {
         modfit1<-findBestModelFunc(posdat,modelTry[mod])[[1]]
         predcpue<-makePredictions(bin1,modfit1,modelTry[mod],datout)
         rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
         metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
       }
     }
   }
   for(mod in which(!modelTry %in% c("Lognormal","Gamma"))) {
     if(modelFail[run,modelTry[mod]]=="-") {
       modfit1<-findBestModelFunc(datin,modelTry[mod])[[1]]
       predcpue<-makePredictions(modfit1,modType=modelTry[mod], newdat = datout)
       rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
       metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
     }
   }
  }
# Calculate RMSE and ME
 modelTable[[run]]$RMSE[-1]<-apply(rmsetab[[run]],2,mean,na.rm=TRUE)
 modelTable[[run]]$ME[-1]<-apply(metab[[run]],2,mean,na.rm=TRUE)
 printTableFunc("Run summary",sp[run],modelTable[[run]],paste0(outVal,"modelSummary.pdf"),useRowNames = FALSE)
 write.csv(residualTab[[run]],paste0(outVal,"modelSummary.csv"))
 write.csv(rmsetab[[run]],paste0(outVal,"rmse.csv"))
 write.csv(metab[[run]],paste0(outVal,"me.csv"))
 plotCrossVal(rmsetab[[run]],metab[[run]],paste0(outVal,"CrossValidation.pdf"))
 #Select best model based on cross validation
 best<-which(!is.na( modelTable[[run]]$RMSE) & modelTable[[run]]$RMSE==min(modelTable[[run]]$RMSE,na.rm=TRUE))
 bestmod[run]<-modelTry[best]
 predbestmod[[run]]<-predvals[[modelTry[best]]]
 plotFits(predbestmod[[run]],bestmod[run],paste0(outVal,"BestFit.pdf"))
 modelTable[[run]]
 }
rm(list=c("bin1","predvals","modfits"))
print(paste(run, common[run],"complete"))
}
source(paste0(baseDir,"/5.FinalModelPrinting.r"))
###################### End of analysis loop #######################

#Save R workspace
if(saveR) save.image(file=paste0(outDir,"/R.workspace.rData"))

  