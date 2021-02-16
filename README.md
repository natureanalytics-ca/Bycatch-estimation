

## Bycatch estimation tool

These R files run a generic bycatch estimation tool. See the User Guide (UserGuide.pdf or UserGuide.Rmd) for details. The .r files are:

1.BycatchModelSpecification.r
2.BycatchModels.r
3.BycatchFunctions.r
4.PreliminaryDataSummary.r
5.FinalModelPrinting.r
6.ExamineModelResults.r

To use the tool, give the names of your data files, and the desired settings in
1.BycatchModelSpecification.r, and then run 2.BycatchModels.r in its entirety.
You should not need to modify anything except in 1.BycatchModelSpecification.r. The files called examplelog.csv and exampleobs.csv are simulated data examples loosely based on the Gulf of Mexico reef fish longline fishery, which will work with the example. 

The file labeled A.ReadDataLL.r sets up the data for analysis from the actual reef longline observer and logbook data, not included due to privacy concerns. 
