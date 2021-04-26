# Modify this file to specify how the bycatch estimation model should be set up. 
# This will automatically be saved in the output directory as documentation.

#Specify directory where R files are found. 
baseDir<-"C:/Users/ebabcock/Dropbox/bycatch project/Current R code"
setwd(baseDir)

#Give a name to the run, which will be used to set up a directory for the the outputs
runName<-"VL catch days"

#### Read in the observer data file. This could also be an assignment to an
# R object if you have already got the data in R.
obsdat<-obsvltrip
  
#### Specify name of logbook data file. Can be aggregated or may include one line per sample unit
logdat<-logtripvl

### What is the sample unit in this data set? e.g. sets or trips. 
sampleUnit<-"trips" #Usually trips

#Specify the name of the effort variable in the observer data and logbook data. These must be
#in the same units. (e.g. 1000 hook hours). Also specify a column for effort
#that is not sampled, in trips with observers. This can be zero in all cases if observers
#sample 100% of effort in sampled trips. 
obsEffort<-"seadays" #"decDaysSampled"
#obsEffortNotSampled<-"unsampled.sets"  #Not needed yet. May be added to later version
logEffort<-"away"

# Give the name of the column in the logbook data that gives the number of sample units (trips or sets)
# that each row includes. If the logbook data is not aggregated (i.e. each row is a sample unit) 
# there should be a column with values of 1 to indicate that each column is one sample unit. 
logNum<-"trips"

#Give common and scientific names. Can be a vector of names to do multiple species at the same time
#in which case each species must have its own column in obsdat.
#This example takes the columns from a data frame to run multiple species at once.
common<-spSummary$COMMON
sp<-spSummary$SCIENTIFIC

#Give the name of the columns associated with the species. If it is a vector, 
# it must match the common and scientific names given above
obsCatch<-paste0("catch.",sp)

# Give units and type of catch to go in plot labels. Must be a vector of the same length as sp
catchUnit<-rep("kg",length(sp))
catchType<-rep("kept", length(sp))
#Specify the name of the variable defining the Years in both databases
yearVar<-"Year"

# Specify the most complex and simplest model to be considered. The code will find compare all
# intermediate models using information criteria
complexModel<-formula(y~(Year+season+EW)^2)
simpleModel<-formula(y~Year)

## The variables must have identical names and factor levels in the observer and logbook datasets
#Specify which of these variables should be interpreted as categorical to make sure they are in factor format.
#Variables not in this list will retain their original format
factorNames=c("Year","season","EW")

#Specify which observation error models to try. Options are delta-lognormal, delta-gamma, negative 
#binomial and tweedie, specified as: "Lognormal" for delta lognormal,"Gamma" for delta gamm,"NegBin" for negative binomial 
#using glm.mb in the MASS library, "Tweedie" for cpglm, and TMB nbinom1, nbinom2, and tweedie in the glmmTMB 
#library, specified with "TMB" followed by the model type. Binomial is run
#automatically as part of the delta models. 
modelTry<-c("Lognormal","Gamma","NegBin","Tweedie","TMBnbinom1","TMBnbinom2","TMBtweedie")

#Specify preferred information criteria for model selection
# Choices are AICc, AIC and BIC. 
selectCriteria<-"BIC"

#Specify whether to run a 10 fold cross-validation (TRUE or FALSE). This may not work with a small 
#or unbalanced dataset.
DoCrossValidation<-TRUE

#Specify whether to exclude models that fail the DHARMa residuals test. 
ResidualTest<-TRUE

# Specify whether to save R workspace. This should be true unless you 
# don't have space on your disk. 

saveR<-TRUE

