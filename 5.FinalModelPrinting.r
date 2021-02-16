#Print summary of which model configurations failed
write.csv(modelFail,paste0(outDir,"/modelFail.csv"))
modelFail
x=gt(data.frame(modelFail),rownames=TRUE)%>%
  tab_header(title = "Summary of model Failure mode")  
x
gtsave(x,file=paste0(outDir,"/modelFail.pdf"),zoom=1)

#Combine pdfs into single files for ease of use
for(run in 1:numSp) {
  outFile<-paste0(dirname[[run]],common[run],catchType[run],"CompleteOutput.pdf")
  if(file.exists(outFile)) file.remove(outFile)
  pdf_combine(outFiles[[run]][file.exists(outFiles[[run]])],outFile)
  outFile2<-paste0(outDir,"/",common[run],catchType[run],"CompleteOutput.pdf")
  if(file.exists(outFile2)) file.remove(outFile2)
  file.copy(outFile,outFile2)
}

outFile<-paste0(outDir,"/BestFits.pdf")
if(file.exists(outFile)) file.remove(outFile)
bestfits<-NULL
for(run in 1:numSp) {
  bestfits=c(bestfits,outFiles[[run]][grep("BestFit",outFiles[[run]])])
}  
pdf_combine(bestfits[file.exists(bestfits)],outFile)
outFile<-paste0(outDir,"/CrossVal.pdf")
if(file.exists(outFile)) file.remove(outFile)
crossvals<-NULL
for(run in 1:numSp) {
  crossvals=c(crossvals,outFiles[[run]][grep("Crossvalidation",outFiles[[run]])])
}  
pdf_combine(crossvals[file.exists(crossvals)],outFile)

# #If you are comparing to known values of total catch use this to plot
# for(run in 1:numSp) {
#   logValidate<-logyearsum %>% filter(Year %in% 2009:2017) %>%
#                               rename(Total=!!paste0("catch.",sp[run]))  %>%
#                               dplyr::select(Year,Total)
#   plotFitsValidate(predbestmod[[run]],logValidate,paste0(outDir,"/Model Summary All Species Validate.pdf"))
# }
