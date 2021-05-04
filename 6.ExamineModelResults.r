#This code allows the user to select a particular model and look at its outputs

# If the model results are not currently loaded, you may load them by 
# uncommenting the following lines.
#specFile<-"C:/Users/ebabcock/Dropbox/bycatch project/Current R code/1.BycatchModelSpecification.r"
#source(specFile)
# load("C:/Users/ebabcock/Dropbox/bycatch project/Current R code/OutputSimulated data example/R.workspace.rData")

#Select the model you want to run
#Options for model types are: c("Lognormal","Gamma","NegBin","Tweedie","TMBnbinom1","TMBnbinom2", "TMBtweedie")
ModelSelect<-"Gamma"  
#Specify a model formula for the model, or for the positive catch portion for delta models
FormulaSelect<- as.formula(y~area + fleet + Year)  
#Specify a model for the binomial portion of delta models (Not used for others)
FormulaBin<- as.formula(y~habBUM + hbf + season + 1 + area + fleet + Year)  
#Give the number of the species to pull out the correct data
SpeciesSelect<-1  

# No changes from here
datval<-dat[[SpeciesSelect]]
selectOutDir<-paste0(outDir,"/",common[SpeciesSelect],catchType[SpeciesSelect],"Selected")
if(!dir.exists(selectOutDir)) dir.create(selectOutDir)
modlist<-FitModelFunc(FormulaBin,FormulaSelect,ModelSelect,datval,selectOutDir)
predvals<-makePredictionsVar(modlist[[1]],modlist[[2]],modType=ModelSelect,newdat=logdat)   
write.csv(predvals,paste0(selectOutDir,"/",ModelSelect,"AnnualSummary.csv"))
plotFits(predvals,ModelSelect,paste0(selectOutDir,"/",ModelSelect,"SelectedModel.pdf"))
if(ModelSelect %in% c("Lognormal","Gamma")) ResidualsFunc(modlist[[1]],"Binomial",paste0(selectOutDir,"/",ModelSelect,"ResidualsBin.pdf")) else
 ResidualsFunc(modlist[[1]],ModelSelect,paste0(selectOutDir,"/",ModelSelect,"Residuals.pdf"))
if(!is.null(modlist[[2]])) ResidualsFunc(modlist[[2]],ModelSelect,paste0(selectOutDir,"/",ModelSelect,"Residuals.pdf"))

