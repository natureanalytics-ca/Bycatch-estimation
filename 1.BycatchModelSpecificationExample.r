# Modify this file to specify how the bycatch estimation model should be set up. 
# This will automatically be saved in the output directory as documentation.

#Specify directory where R files are found. 
baseDir<-getwd()
#setwd(baseDir)

#Give a name to the run, which will be used to set up a directory for the the outputs
#and a run description to describe the run in the output files
runName<-"SimulatedExample"
runDescription<-"Example with simulated data"

# What would you like to estimate?
# You may calculate either an annual abundance index, or total bycatch, or both
# If you want total bycatch, you must have logbook data or some other source of total effort
EstimateIndex<-TRUE
EstimateBycatch<-TRUE

#### Read in the observer data file. This could also be an assignment to an
# R object if you have already got the data in R.
obsdat<-read.csv("ExampleoBs.csv",as.is=TRUE) 
 
#### Specify name of file with total effort data, if you are estimating total bycatch. 
# It can be aggregated or may include one line per sample unit (trip). 
# If estimateByatch is FALSE, this variable is not needed
logdat<-read.csv("ExampleLog.csv",as.is=TRUE) 

### What is the sample unit in this data set? e.g. sets or trips. 
sampleUnit<-"trips" #Usually trips, may be sets

#Specify the name of the effort variable in the observer data and logbook data. These must be
#in the same units. (e.g. 1000 hook hours). Also specify a column for effort
#that is not sampled, in trips with observers. This can be zero in all cases if observers
#sample 100% of effort in sampled trips. 
obsEffort<-"sampled.sets"
logEffort<-"sets" #This variable is only needed if estimating bycatch
logUnsampledEffort<-NULL  #This variable is only needed includeObsCatch is true

# Give the name of the column in the logbook data that gives the number of sample units (trips or sets)
# that each row includes. If the logbook data is not aggregated (i.e. each row is a sample unit) 
# just make this NA 
logNum<-NA

# Make includeObsCatch TRUE if (1) the observed sample units can be matched to the logbook
# sample units and (2) you want to calculate total bycatch as the observed bycatch plus the 
# predicted unobserved bycatch. This doesn't work with aggregated logbook effort. 
includeObsCatch<-FALSE
# if includeObsCatch is true, give the name of the column that matches sample
# units between the observer and logbook data. Otherwise, this can be NA
matchColumn<-NA

#Give common and scientific names. Can be a vector of names to do multiple species at the same time
#in which case each species must have its own column in obsdat.
#This example takes the columns from a data frame to run multiple species at once.
common<-"Simulated species"
sp<-"Genus species"

#Give the name of the columns associated with the species. If it is a vector, 
# it must match the common and scientific names given above
obsCatch<-"Catch"

# Give units and type of catch to go in plot labels. Must be a vector of the same length as sp
catchUnit<-rep("number",length(sp))
catchType<-rep("dead discard", length(sp))
#Specify the name of the variable defining the Years in both databases
yearVar<-"Year"

# Specify the most complex and simplest model to be considered. The code will find compare all
# intermediate models using information criteria. Use indexModel to specify which strata to 
# keep separate in calculating abundance indices. 
complexModel<-formula(y~(Year+season)^2)
simpleModel<-formula(y~Year)
indexModel<-formula(y~Year)

## The variables must have identical names and factor levels in the observer and logbook datasets
#Specify which of these variables should be interpreted as categorical to make sure they are in factor format.
#Variables not in this list will retain their original format
factorNames=c("Year","season")

#Specify which observation error models to try. Options are: "Binomial", "Normal","Lognormal",
#"Gamma","Delta-Lognormal", #"Delta-Gamma", "NegBin" for Negative binomial" using glm.mb in the MASS library,
#"Tweedie" for Tweedie GLM from the cpglm library, and "TMBnbinom1", "TMBnbinom2", and "TMBtweedie" for
# negative binomial 1, negative binomial 2 and Tweedie from the GLMMTMB library. Binomial is run
#automatically as part of the delta models if either of them are selected.

modelTry<-c("Binomial","Lognormal","Delta-Lognormal","Delta-Gamma","NegBin","TMBnbinom1","TMBnbinom2","TMBtweedie")

#Specify preferred information criteria for model selection
# Choices are AICc, AIC and BIC. 
selectCriteria<-"BIC"

#Specify whether to run a 10 fold cross-validation (TRUE or FALSE). This may not work with a small 
#or unbalanced dataset. DredgeCrossValidation specifies whether to use information criteria 
#to find the best model in cross validation, using the dredge function, or just keep the same model formula. 
#Do not use dredge for very large datasets, as the run will be slow.  
DoCrossValidation<-TRUE
DredgeCrossValidation<-FALSE

#Specify whether to exclude models that fail the DHARMa residuals test. 
ResidualTest<-FALSE

#Specify confidence interval for total bycatch estimates
CIval<-0.05 #Should be the alpha level, e.g. 0.05 for 95%

# Variance calculation method. Simulate will not work with a large number of sample units
# in the logbook data. The delta method for variance calculation is not implemented
# for the delta-lognormal or delta-gamma methods. 
VarCalc<-c("Simulate","DeltaMethod","None")[1] 

# Specify whether to save R workspace. This should be true unless you 
# don't have space on your disk. Also specify whether to use parallel processing to
# speed up calculations.
saveR<-TRUE
useParallel<-TRUE

## Validation. If you have true values of the total bycatch (for example in a simulation study)
# Make PlotValidation true and fill out the rest of the specification. 
plotValidation<-FALSE
trueVals<-NULL
trueCols<-NULL

