#This file allows you to load multiple workspaces and compare the results. 
#This is mostly useful for comparing different variance estimates. 
#Alter the workspaces loaded and their names as needed

#Load multiple workspace and rename the predictions for combination
load("~/Box Sync/bycatch project (ebabcock@miami.edu)/Current R code/OutputSimulated factor year delta var/R.workspace.rData")
Delta<-allmods[[1]]
load("~/Box Sync/bycatch project (ebabcock@miami.edu)/Current R code/OutputSimulated factor year sim var/R.workspace.rData")
Simvar<-allmods[[1]]
load("~/Box Sync/bycatch project (ebabcock@miami.edu)/Current R code/OutputSimulated factor year no var/R.workspace.rData")
Novar<-allmods[[1]]

#Combine predictions
allmods<-bind_rows(DeltaMethod=Delta,Simulated=Simvar,Novar=Novar,.id="Type")
allmods$Source=paste(allmods$Type,allmods$Source)
#Adjust numeric year if needed
#allmods$Year[allmods$Year<startYear]=allmods$Year[allmods$Year<startYear]+startYear
table(allmods$Source)
names(allmods)
summary(allmods)
#Plot combinations of models
plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod NegBin","Simulated NegBin","Novar NegBin") ),"All",NULL)
plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod NegBin","Simulated NegBin","DeltaMethod TMBnbinom2","Simulated TMBnbinom2","Novar TMBnbinom2","Novar NegBin") ),"All",NULL)
plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod Tweedie","Simulated Tweedie","Novar Tweedie") ),"All",NULL)
plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod Tweedie","Simulated Tweedie", "DeltaMethod TMBtweedie", "Simulated TMBtweedie","NovarTweedie","Novar TMBtweedie","Novar Tweedie") ),"All",NULL)
plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod TMBnbinom1","Simulated TMBnbinom1","Novar TMBnbinom1") ),"All",NULL)
plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod Gamma","Simulated Gamma","Novar Gamma") ),"All",NULL)

plotSums(filter(allmods,Valid==1 & Source %in% c("Simulated Delta-Lognormal","Novar Delta-Lognormal") ),"All",NULL)
plotSums(filter(allmods,Valid==1 & Source %in% c("Simulated Delta-Gamma","Novar Delta-Gamma") ),"All",NULL)

plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod Lognormal","Simulated Lognormal","Novar Lognormal") ),"All",NULL)

plotSums(filter(allmods,Valid==1 & Source %in% c("DeltaMethod Normal","Simulated Normal","Novar Normal") ),"All",NULL)


