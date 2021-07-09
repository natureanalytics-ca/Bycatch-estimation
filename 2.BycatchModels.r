# The observer data should be aggregated to the appropriate sample unit, either trips or sets.
# Effort must be in the same units in both data sets (e.g. hook hours), but it does not matter
# how the logbook effort data is aggregated, as long has it includes data on all stratification or 
# predictor variables. For example, the data can be aggregated by year, East-West and season if those
# are the stratification variables, or it can be aggregated by trip.
# The user must specify the names of the databases and the variables to include in the 
# The .r file named on line 12. Everything else should work automatically.
# Contact Beth Babcock ebabcock@rsmas.miami.edu for assistance. 

############### Step 1. Enter the data specification in the file named here #############################
specFile<-"C:/Users/ebabcock/Dropbox/bycatch project/Current R code/1.BycatchModelSpecificationExample.r"
# Either set the working directory or put the full patch in the filename.
# Complete the information in the file before continuing. You may run through specFile line by line, but it 
# will also be sourced again later. The file will be saved, with the addition of the date, to the output directory

##################################################################################
### From here on no changes should be needed.  
########################################################################
#Read in data specification
source(specFile)
#Make data summary outputs
source(paste0(baseDir,"/4.preliminarySetup.r"))
# Stop here and check data summaries to make sure each species has reasonable number
# of observations in each year for analysis.
#The file called  "Data summary all species.pdf" has all the tables. 
#There is also a .csv for each species. 

############## This is the main analysis loop ##########################
StartTime<-Sys.time()   
for(run in 1:numSp) {
 datval<-dat[[run]]
 outVal<-dirname[[run]]
 varExclude<-NULL
 #Fit all models except delta
 for(mod in which(!modelTry %in% c("Lognormal","Gamma"))){
   modfit1<-findBestModelFunc(datval,modelTry[mod],printOutput=TRUE)
   modelSelectTable[[run]][[modelTry[mod]]]<-modfit1[[2]]
   modFits[[run]][[modelTry[mod]]]<-modfit1[[1]]
 }
#Fit delta models
 if("Lognormal" %in% modelTry | "Gamma" %in% modelTry) {  #Delta models if requested
   posdat<-filter(dat[[run]],pres==1)
   y<-unlist(lapply(posdat[,factorNames],function(x) length(setdiff(levels(x),x)))) #See if all levels are included
   varExclude<-names(y)[y>0]
   if(length(varExclude>0)) print(paste(common[run], "excluding variable",varExclude,"from delta models for positive catch"))
   if((min(summary(posdat$Year))>0 |  is.numeric(datval$Year)) &!is.null(modFits[[run]][["Binomial"]])) { #If all years have at least one positive observation and binomial converged, carry on with delta models
     for(mod in which(modelTry %in% c("Lognormal","Gamma")))  {
       modfit1<-findBestModelFunc(posdat,modelTry[mod],printOutput=TRUE)
       modelSelectTable[[run]][[modelTry[mod]]]<-modfit1[[2]]
       modFits[[run]][[modelTry[mod]]]<-modfit1[[1]]
     }
   } else {
     print("Not all years have positive observations, skipping delta models")
     modelFail[run,c("Lognormal","Gamma")]<-"data"
     modPredVals[[run]][[modelTry[mod]]]<-NULL
     modIndexVals[[run]][[modelTry[mod]]]<-NULL
   }
  }
 #Make predictions, residuals, etc. for all models
 for(mod in 1:length(modelTry)) {
   if(!is.null(modFits[[run]][[modelTry[mod]]])) {
     if(modelTry[mod] %in% c("Lognormal","Gamma")) {
        modfit1<-modFits[[run]][["Binomial"]]
        modfit2<-modFits[[run]][[modelTry[mod]]]
     } else {
        modfit1<-modFits[[run]][[modelTry[mod]]]
        modfit2<-NULL
     }
     if(EstimateBycatch) {
      modPredVals[[run]][[modelTry[mod]]]<-makePredictionsVar(modfit1=modfit1,modfit2=modfit2,modType=modelTry[mod],newdat=logdat,printOutput=TRUE)   
     }
     if(EstimateIndex) {
      modIndexVals[[run]][[modelTry[mod]]]<-makeIndexVar(modfit1=modfit1,modfit2=modfit2,modType=modelTry[mod],printOutput=TRUE)   
     }
     modelTable[[run]]$formula[mod]<-
       paste(formula(modFits[[run]][[modelTry[mod]]]))[[3]]
     temp<-ResidualsFunc(modFits[[run]][[modelTry[mod]]],modelTry[mod],paste0(outVal,"Residuals",modelTry[mod],".pdf"))
     if(!is.null(temp)) {
       residualTab[[run]][,modelTry[mod]]<-temp
       if(residualTab[[run]]["KS.p",modelTry[mod]]<0.01 & ResidualTest) modelFail[run,modelTry[mod]]<-"resid"
     }
     if(is.null(modPredVals[[run]][[modelTry[mod]]])) modelFail[run,modelTry[mod]]<-"cv"
   } else {
     if(modelFail[run,modelTry[mod]]=="-") modelFail[run,modelTry[mod]]<-"fit"
   }
 }
 #Combine all predictions, except Binomial
  if(EstimateBycatch) {
   yearsumgraph<-yearSum[[run]] %>% dplyr::select(Year=Year,Total=CatEst,Total.se=Catse) %>%
     mutate(TotalVar=Total.se^2,Total.cv=Total.se/Total,TotalFixed=Total)
   allmods[[run]]<-bind_rows(c(modPredVals[[run]],list(Ratio=yearsumgraph)),.id="Source") %>%
     filter(!Source=="Binomial")
   allmods[[run]]$Valid<-ifelse(modelFail[run,match(allmods[[run]]$Source,dimnames(modelFail)[[2]])]=="-" | allmods[[run]]$Source=="Ratio",1,0)
  } 
  if(EstimateIndex) {
   allindex[[run]]<-bind_rows(modIndexVals[[run]],.id="Source") %>%
     filter(!Source=="Binomial")
   allindex[[run]]$Valid<-ifelse(modelFail[run,match(allindex[[run]]$Source,dimnames(modelFail)[[2]])]=="-" ,1,0)
  } 
 #Print the diagnostic tables
  write.csv(residualTab[[run]],paste0(outVal,"residualDiagnostics.csv"))
  write.csv(modelFail,paste0(outDir,"/modelFail.csv"))
######## Cross validation 10 fold  ####################################
 rmsetab[[run]]<-data.frame(matrix(NA,10,length(modelTry),dimnames=list(1:10,modelTry)))
 rmsetab[[run]]<-rmsetab[[run]][,colnames(rmsetab[[run]])!="Binomial"]
 metab[[run]]<-rmsetab[[run]]
 if(DoCrossValidation & length(which(modelFail[run,colnames(modelFail)!="Binomial"]=="-"))>0) {  #Don't do unless at least one model worked
  if(NumCores>3 & useParallel)  {
    cl<-makeCluster(NumCores-2)
    registerDoParallel(cl)
  }
  datval$cvsample<-sample(rep(1:10,length=dim(datval)[1]),replace=FALSE)
  table(datval$cvsample,datval$Year)
  foreach(i=1:10 ) %do% {
   datin<-datval[datval$cvsample!=i,]
   datout<-datval[datval$cvsample==i,]
   datout$SampleUnits<-rep(1,dim(datout)[1])
   for(mod in which(!modelTry %in% c("Lognormal","Gamma"))) {
     if(modelFail[run,modelTry[mod]]=="-") {
       if(DredgeCrossValidation) modfit1<-findBestModelFunc(datin,modelTry[mod])[[1]] else
        modfit1<-FitModelFuncCV(formula(paste0("y~",modelTable[[run]]$formula[mod])),modType=modelTry[mod],obsdatval=datin)
       if(modelTry[mod]!="Binomial") {
        predcpue<-makePredictions(modfit1,modType=modelTry[mod], newdat = datout)
        rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
        metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
       } else {
        bin1<-modfit1    
       }
     }
   }
   if("Lognormal" %in% modelTry | "Gamma" %in% modelTry) { 
     posdat<-filter(datin,pres==1)
     for(mod in which(modelTry %in% c("Lognormal","Gamma"))) {
       if(modelFail[run,modelTry[mod]]=="-" & min(summary(posdat$Year))>0) {
         if(DredgeCrossValidation) modfit1<-findBestModelFunc(posdat,modelTry[mod])[[1]] else
          modfit1<-FitModelFuncCV(formula(paste0("y~",modelTable[[run]]$formula[mod])),modType=modelTry[mod],obsdatval=posdat)
         predcpue<-makePredictions(bin1,modfit1,modelTry[mod],datout)
         rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
         metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
       }
     }
   }
  }
 if(NumCores>3) stopCluster(cl)
# Calculate RMSE and ME
 modelTable[[run]]$RMSE[modelTable[[run]]$model!="Binomial"]<-apply(rmsetab[[run]],2,mean,na.rm=TRUE)
 modelTable[[run]]$ME[modelTable[[run]]$model!="Binomial"]<-apply(metab[[run]],2,mean,na.rm=TRUE)
 write.csv(residualTab[[run]],paste0(outVal,"modelSummary.csv"))
 write.csv(rmsetab[[run]],paste0(outVal,"rmse.csv"))
 write.csv(metab[[run]],paste0(outVal,"me.csv"))
 #Select best model based on cross validation
   best<-which(!is.na( modelTable[[run]]$RMSE) & 
     modelTable[[run]]$RMSE==min(modelTable[[run]]$RMSE,na.rm=TRUE))
  if(length(best) >0 ) {
   bestmod[run]<-modelTry[best]
   predbestmod[[run]]<-modPredVals[[run]][[modelTry[best]]]
   indexbestmod[[run]]<-modIndexVals[[run]][[modelTry[best]]]
 } else {
   bestmod[run]<-"None"
 }  
}
save(list=c("numSp","yearSum","runName", "common", "sp","bestmod",
  "predbestmod","indexbestmod","allmods","allindex","modelTable",
 "modelSelectTable","modFits","modPredVals"
 ,"modIndexVals","modelFail","rmsetab","metab",
 "residualTab" ,"run","modelTry","EstimateIndex","EstimateBycatch",
  "DoCrossValidation","indexVarNames","selectCriteria","sampleUnit",
   "modelTry","catchType","catchUnit","residualTab",
   "plotValidation","trueVals","trueCols"),
  file=paste0(outVal,"/","resultsR"))
rmarkdown::render("6.PrintResults.rmd",
   params=list(outVal=outVal))
file.rename(paste0(getwd(),"/6.PrintResults.pdf"),paste0(outVal,"/",common[run],catchType[run],"results.pdf"))
print(paste(run, common[run],"complete, ",Sys.time()))
}
###################### End of analysis loop #######################

#Save R workspace
if(saveR) save.image(file=paste0(outDir,"/R.workspace.rData"))
Sys.time()
Sys.time()-startTime
