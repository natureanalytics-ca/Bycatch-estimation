# !diagnostics off  

# Basic ratio estimator with variance (Cochran)
# x, y and g are vectors giving the effort/catch, bycatch and stratum of 
# each observed sample unit. X is the total effort/catch by stratum 
# N is the total number of sample units by stratum, if available, othewise total effort 
# Output is mean and standard error of bycatch by stratum and the total 
# bycatch with SE. Assumes unobserved strata have zero catch
ratio.func= function(x,y,g,X,N,G) {
  N=N[G %in% g]
  X=X[G %in% g]
  G=G[G %in% g]
  a=order(G)
  X=X[a]
  N=N[a]
  G=G[a]
  if(sum(X-N)!=0)  n=tapply(x,g,length)  else  n=tapply(x,g,sum)
  f=n/N
  xmean=tapply(x,g,mean)
  ymean=tapply(y,g,mean)
  Rhat=ymean/xmean
  sx2=var(x)
  sy2=var(y)
  sxy=cov(x,y)
  stratum.est=X*Rhat
  stratum.var=X^2*(1-f)/(n*xmean^2)*(sy2+Rhat^2*sx2-2*Rhat*sxy)
  total.est=sum(stratum.est)
  total.var=sum(stratum.var)
  list(stratum.est=stratum.est,stratum.se=sqrt(stratum.var),coverage=f,total.est=total.est,total.se=sqrt(total.var))
}


# Statified mean estimator. y and g are observed values of the variable and 
# associated stratum identifier. N is the total sample size in each stratum (Cochran)
stratified.func=function(y,g,N) {
  ymean=tapply(y,g,mean)
  n=tapply(x,g,length) 
  yvar=tapply(y,g,var)
  stratum.est=ymean*N
  stratum.var=N*(N-n)*yvar/n
  total.est=sum(stratum.est)
  total.var=sum(stratum.var)
  list(stratum.est=stratum.est,stratum.se=sqrt(stratum.var),total.est=sum(stratum.est))
}

# Function to convert data in excel format with date and time separated by a blank
# into an R format date
getdatefunc=function(x,dateformat="%m/%d/%Y") {
  require(reshape2)
  y=colsplit(as.character(x)," ",names=c("date","time"))
  as.Date(y[,1],format=dateformat)
}  
getdatefunc2=function(x,dateformat="%d%b%Y") {
  require(reshape2)
  y=colsplit(as.character(x),":",names=c("date","hour","sec"))
  as.Date(y[,1],format=dateformat)
}  
gettimefunc2=function(x) {
  require(reshape2)
  y=colsplit(as.character(x),":",names=c("date","hour","min","sec"))
  timeval=y[,2]+y[,3]/60+y[,4]/3600
  timeval
}  

seasonfunc<-function(month,numseason=4) {
  if(!is.numeric(month)) month=as.numeric(as.character(month))
  seasons=rep(1:numseason,each=12/numseason)
  seasons[month]
}

# Function to count the number of unique levels in a vector
length.unique=function(x) length(unique(x))

# Stratum designations from Scott-Denton paper, and shrimp observer manual for GOM shrimp areas
areafunc=function(x) {
  x$StatZone[is.na(x$StatZone)]=-1
  area=rep(-1,dim(x)[1])  
  area[x$StatZone>=1 & x$StatZone<=9]=1 #"WFL"
  area[x$StatZone>=10 & x$StatZone<=12]=2 #"AL-MS"
  area[x$StatZone>=13 & x$StatZone<=17]=3 #"LA"
  area[x$StatZone>=18 &x$StatZone<=21]=4 #"TX"
  x$LatInS[is.na(x$LatInS)]=0
  x$LatInM[is.na(x$LatInM)]=0
  x$LatOutS[is.na(x$LatOutS)]=0
  x$LatOutM[is.na(x$LatOutM)]=0
  x$lat=x$LatInD+x$LatInM/60+x$LatInS/3600
  x$lat[is.na(x$lat)]=(x$LatOutD+x$LatOutM/60+x$LatOutS/3600)[is.na(x$lat)]
  area[x$StatZone>=24 & x$StatZone<=29 ]=5 #"EFL"
  area[x$StatZone==30 & x$lat<30.708]=5 #"EFL"
  area[x$StatZone==30 & is.na(x$lat)]=5 #"EFL"
  area[x$StatZone>=30 & x$lat>=30.708 ]=6 #"GA"
  area[x$StatZone==31 ]=6 #"GA"
  area[x$StatZone==32 ]=7 #"SC"
  area[x$StatZone==33 & x$lat<33.86]=7 #"SC"
  area[x$StatZone==33 & x$lat>=33.86]=8 #"NC"
  area[x$StatZone>=34]=8 #"NC"
  area[is.na(x$StatZone) & x$LonInD>88 & x$LonInD<89]=2
  area
}

# Function to find mode of a categorical variable
mostfreqfunc<-function(x) {
  x=x[!is.na(x)]
  if(length(x)>0) {
   y=aggregate(x,list(x),length)
   temp=y$Group.1[y$x==max(y$x)][1]
  } else temp=NA
  temp
}

# Function to find range of a numerical variable
getRange<-function(x) {
  max(x,na.rm=TRUE)-min(x,na.rm=TRUE)
}

## Function to convert new areas to old areas from Kevin McCarthy
areaGOM=function(x) {
  x[x>=2383& x<= 2384]=2
  x[x >= 2483 & x <= 2485]=2
  x[x >= 2581 & x <= 2585]=3
  x[x >= 2681 & x <= 2685]=4
  x[x >= 2782 & x <= 2785]=5
  x[x >= 2882 & x <= 2884]=6
  x[x >= 2982 & x <= 2984]=7
  x[x >= 3083 & x <= 3084]=7
  x[x==3085 | x==2985 | x==2885]=8
  x[x %in% c(3086,2986,2886,2786,2686,2586,2486)] = 9
  x[x %in% c(3087,2987,2887,2787,2687,2587)] = 10
  x[x %in% c(3088,2988,2888,2788,2688,2588)] = 11
  x[x %in% c(3089,3090)] = 12
  x[x %in% c(2989,2889,2789,2689,2589)] = 13
  x[x %in% c(2990,2890,2790,2690,2590)] = 14
  x[x %in% c(2991,2891,2791,2691,2591)] = 15
  x[x %in% c(2992,2892,2792,2692,2592)] = 16
  x[x %in% c(2993,2893,2793,2693,2593)] = 17
  x[x %in% c(2994,2894,2794,2694,2594)] = 18
  x[x %in% c(2995,2895,2896)] = 19
  x[x %in% c(2797,2796,2795)] = 20
  x[x %in% c(2697,2696,2695)] = 21
  x[x>1000]=NA
  x
}

## Function to count outliers, defined as more than numSD standard deviations from the mean.
outlierCountFunc=function(x,numSD=8) {
  length(which(x>mean(x,na.rm=TRUE)+numSD*sd(x,na.rm=TRUE) |  x<mean(x,na.rm=TRUE)-numSD*sd(x,na.rm=TRUE)))
}

##Function to calculate number of hours fished in each day counted as first set to last haul
fishTimeFunc<-function(timeout,timein,prop.sampled=1) { 
  timeout=ifelse(timeout>timein & !is.na(timein+timeout),timeout,timeout+24)
  maxtimeday=(max(timeout,na.rm=TRUE)-min(timein,na.rm=TRUE))
  if(maxtimeday<0.5) maxtimeday=0.5   #Arbitrary minimum
  meanprop=mean(prop.sampled,na.rm=TRUE)
  if(is.na(meanprop)) meanprop=1  #based on median proportion sampled=1
  maxtimeday*meanprop
}

fishTimeFunc<-function(timeout,timein,prop.sampled=1) { 
  timeout=ifelse(timeout>timein & !is.na(timein+timeout),timeout,24)
  maxtimeday=(max(timeout,na.rm=TRUE)-min(timein,na.rm=TRUE))
  if(maxtimeday<0.5) maxtimeday=0.5   #Arbitrary minimum
  meanprop=mean(prop.sampled,na.rm=TRUE)
  if(is.na(meanprop)) meanprop=1  #based on median proportion sampled=1
  maxtimeday*meanprop
}


#######Delta lognormal functions
#Variance of a product from Lo et al. (1992)
lo.se=function(x1,x1e,x2,x2e) {
  cor12=cor(x1[!is.na(x1) &!is.na(x2)],x2[!is.na(x1) &!is.na(x2)])
  (x1e^2 * x2^2 +x2e^2*x1^2+2*x1*x2*cor12*x1e*x2e)^0.5
}

#Calculate normal mean and standard error from lognormal mean and se
norm.mean=function(x1,x1e) {
  log(x1)-0.5*log(1+(x1/x1e)^2)
}
norm.se=function(x1,x1e) {
  sqrt(log(1+(x1/x1e)^2))
}
#Calculate lognormal mean and standard error from normal mean and se
lnorm.mean=function(x1,x1e) {
  exp(x1+0.5*x1e^2)
}
lnorm.se=function(x1,x1e) {
  ((exp(x1e^2)-1)*exp(2*x1+x1e^2))^0.5
}

#Inverse logit
ilogit=function(x) {
  1/(1+exp(-x))
}

#Exact variance of the product of two independent variables, from Goodman (1960)
goodman.var<-function(x,y) {
  var(x)*y+var(y)*x-var(x)*var(y)
}

#Standard error of a mean
standard.error<-function(x) {
  x=x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#Generate standard errors of predictions from simulation from regression coefficients and their var/covar matrix
getSimSE<-function(modfit, df1, transFunc="none",offsetval=NULL, nsim=nSims) {
  if(length(coef(modfit))==dim(vcov(modfit))[1]) {
    b=t(mvrnorm(nsim,coef(modfit),vcov(modfit)))
    yvar=sub( " ", " ",formula(modfit) )[2]
    df1<-cbind(y=rep(1,dim(df1)[1]),df1)
    names(df1)[1]<-yvar
    a=model.matrix(formula(modfit),data=df1)
    d<- a %*% b
    if(transFunc=="exp") d<-exp(d)
    if(transFunc=="ilogit") d<-ilogit(d)
    e<-apply(d,1,sd)
    if(!is.null(offsetval)) e<-e*data.frame(newdat)[,offsetval]  #For negbin only
  } else {
    e<-rep(NA,dim(df1)[1])
  }
  e
}

#Generate simulations to use as input to DHARMa residual calculations for the Tweedie from cpglm. 
simulateTweedie <- function(modfit1, nsims=nSims){
  pred = predict(modfit1, type = "response")
  nObs = length(pred)
  sim = replicate(nsims,rtweedie(nObs,xi=modfit1$p, mu=pred,phi=modfit1$phi))
  return(sim)
}

#Generate mean and standard error of predictions for delta lognormal by simulation
getSimDeltaLN<-function(modfitBin,modfitLnorm, df1, nsim=10000) {
  b1=t(mvrnorm(nsim,coef(modfitBin),vcov(modfitBin)))
  yvar=sub( " ", " ",formula(modfitBin) )[2]
  df11<-cbind(y=rep(1,dim(df1)[1]),df1)
  names(df11)[1]<-yvar
  a1=model.matrix(formula(modfitBin),data=df11)
  d1<- a1 %*% b1
  b2=t(mvrnorm(nsim,coef(modfitLnorm),vcov(modfitLnorm)))
  yvar=sub( " ", " ",formula(modfitLnorm) )[2]
  df12<-cbind(y=rep(1,dim(df1)[1]),df1)
  names(df12)[1]<-yvar
  a2=model.matrix(formula(modfitLnorm),data=df12)
  d2<- a2 %*% b2
  lnormmean<-apply(exp(d2),1,mean)
  lnormse<-apply(exp(d2),1,sd)
  deltamean<-apply(ilogit(d1)*exp(d2),1,mean)
  deltase<-apply(ilogit(d1)*exp(d2),1,sd)
  data.frame(delta.mean=deltamean,delta.se=deltase,lnorm.mean=lnormmean,lnorm.se=lnormse)
}


# Make tables with headers using gt and print to a graphics file
printTableFunc=function(headerText,subheaderText="",df1,filename,useRowNames=FALSE) {
  rowNames=rownames(df1)
  if(is.matrix(df1) | class(df1)[1]=="model.selection") df1<-as.data.frame(df1)
  numvars=names(select_if(df1,is.numeric))
  numvars=numvars[!numvars %in% c("Year","year")]
  if(useRowNames) df1<-cbind(" "=rowNames,df1)
  x=gt(df1)%>%
    tab_header(title = headerText,
      subtitle=subheaderText) %>%
    fmt_number(
      columns=one_of(numvars),
      n_sigfig = 3) %>%
    tab_options(table.width=720,
      table.font.size=10,
      column_labels.font.size=10,
       heading.title.font.size=12  )
  gtsave(x,file=filename, vwidth=720,vheight=864,zoom=1)
}

# Make pretty table with headers using kable
printTableFunc=function(headerText,subheaderText="",df1,filename,useRowNames=FALSE) {
  if(is.matrix(df1) | class(df1)[1]=="model.selection") 
    df1<-as.data.frame(df1)
  df1<-df1 %>%
  mutate_if(
    is.numeric,
    ~ ifelse(abs(.x) > 1, round(.x), round(.x, 2))
  )
  kable(df1,caption = paste(headerText,subheaderText),format="html") %>%
  kable_styling(full_width=FALSE) %>%
  save_kable(filename)
}

# Make ANOVA tables with headers using GridExtra library and print to a graphics file
printTableFunc=function(headerText,df1,rowNames=FALSE,speciesRun=TRUE) {
 padding <- unit(0.5,"line")
 if(speciesRun) title <- textGrob(paste(headerText,common[run],catchType[run]),gp=gpar(fontsize=12)) else
   title <- textGrob(headerText,gp=gpar(fontsize=12))
 if(rowNames)
  table <-tableGrob(df1,theme=mytheme) else
    table <-tableGrob(df1,theme=mytheme,rows=NULL)
  table <- gtable_add_rows(table, 
                         heights = grobHeight(title) + padding,
                         pos = 0)
 table<- gtable_add_grob(table,title,
                        t=1, l=1, 
                        r=ncol(table))
 grid.newpage()
 grid.draw(table)
}

# Make pretty anova tables with kable
printANOVAFunc2=function(headerText,subheaderText="",df1,filename,useRowNames=TRUE) {
  numcols=dim(df1)[2]
  dust(df1,  keep_rownames=useRowNames) %>%
  sprinkle(col=2:numcols,round=1) %>%
  sprinkle(col=numcols+1,fn=quote(pvalString(value))) %>%
  sprinkle_colnames(.rownames="Term") %>%
  kable(caption = paste0(headerText," ",subheaderText)) %>%
  kable_styling(full_width=FALSE) %>%
  save_kable(file=filename)
}

# Make  anova tables with gt
printANOVAFunc2=function(headerText,subheaderText="",df1,filename,useRowNames=TRUE) {
  numcols=dim(df1)[2]
  dust(df1,  keep_rownames=useRowNames) %>%
  sprinkle(col=2:numcols,round=1) %>%
  sprinkle(col=numcols+1,fn=quote(pvalString(value))) %>%
  sprinkle_colnames(.rownames="Term") %>%
  kable(caption = paste0(headerText," ",subheaderText)) %>%
  kable_styling(full_width=FALSE) %>%
  save_kable(file=filename)
}


#Function to find best model by information criteria, of model types 
#Binomial, Lognormal, Gamma, Negative Binomial and Tweedie.
findBestModelFunc<-function(obsdatval,modType,printOutput=FALSE) {
  keepVars=requiredVarNames
  offset=NULL
  extras=c("AICc","AIC", "BIC")
  if(modType %in% c("NegBin","TMBnbinom1","TMBnbinom2"))     
    obsdatval$y=round(obsdatval$Catch)
  if(modType %in% c("Tweedie","TMBtweedie")) 
    obsdatval$y=obsdatval$cpue
  if(modType =="Binomial")   {  
    obsdatval$y=obsdatval$pres
    funcName="glm"
    args=list(formula="",data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=100),na.action=na.fail)
  }
  if(modType=="Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    funcName="lm"
    args=list(formula="",data=obsdatval,na.action=na.fail)
  }
  if(modType=="Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    funcName="glm"
    args=list(formula="",data=obsdatval,family=Gamma(link="log"),na.action=na.fail)
  }
  if(modType == "NegBin") {
    funcName="glm.nb"
    offset="+offset(log(Effort))"
    keepVars=c(requiredVarNames,"offset(log(Effort))")
    args=list(formula="",data=obsdatval,control=glm.control(epsilon=1E-6,maxit=30),na.action=na.fail)
  }
  if(modType=="Tweedie") {
    funcName="cpglm"
    args=list(formula="",data=obsdatval,na.action=na.fail)
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    funcName="glmmTMB"
    TMBfamily=gsub("TMB","",modType)
    offset="+offset(log(Effort))"
    keepVars=paste0("cond(",c(requiredVarNames,"offset(log(Effort))"),")")
    args=list(formula="",family=TMBfamily,data=obsdatval,na.action=na.fail)
  }
  if(modType %in% c("TMBtweedie") ){
    funcName="glmmTMB"
    TMBfamily=gsub("TMB","",modType)
    keepVars=paste0("cond(",requiredVarNames,")")
    args=list(formula="",family=TMBfamily,data=obsdatval,na.action=na.fail)
  }
  formulaList<-list(as.formula(paste("y~",paste(getAllTerms(complexModel),collapse="+"),offset)),
                    as.formula(paste("y~",paste(allVarNames,collapse="+"),offset)),
                    as.formula(paste("y~",paste(allVarNames[!allVarNames %in% varExclude],collapse="+"),offset)),
                    as.formula(paste("y~",paste(requiredVarNames,collapse="+"),offset)), NA)
  if(length(varExclude)>0 & modType %in% c("Lognormal","Gamma")) {
    formulaList[[1]]<-NULL
    formulaList[[1]]<-NULL
  }  
  args$formula=formulaList[[1]]
  modfit1<-try(do.call(funcName,args))
  for(i in 2:(length(formulaList))-1) {
    if(class(modfit1)[1] %in% c("glm","lm","glm.nb")) {
      if(modfit1$rank<length(coef(modfit1))) class(modfit1)<-"try-error" 
    } 
    if(class(modfit1)[1] == "cpglm")  {
      if(length(coef(modfit1))!=dim(vcov(modfit1))[1])  class(modfit1)<-"try-error" 
    }
    if(class(modfit1)[1]=="try-error")   {
      args$formula=formulaList[[i]]
      modfit1<-try(do.call(funcName,args))
    }
  }
  if(class(modfit1)[1]=="try-error")   {
    returnval=NULL
    print(paste(common[run],modType,"failed to converge"))
  } else {
    if(modType=="Binomial") modfit1<-glm(formula(modfit1),data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=100),na.action=na.fail)
    if(modType=="Lognormal") modfit1<-lm(formula(modfit1),data=obsdatval,na.action=na.fail)
    if(modType=="Gamma") modfit1<-glm(formula(modfit1),data=obsdatval,family=Gamma(link="log"),na.action=na.fail)
    if(modType=="NegBin") modfit1<-glm.nb(formula(modfit1),data=obsdatval,control=glm.control(epsilon=1E-6,maxit=30),na.action=na.fail)
    if(modType=="Tweedie") modfit1<-cpglm(formula(modfit1),data=obsdatval,na.action=na.fail)
    if(modType %in% c("TMBnbinom1","TMBnbinom2","TMBtweedie") ) 
      modfit1<-glmmTMB(formula(modfit1),family=TMBfamily,data=obsdatval,na.action=na.fail)
   if(useParallel) {
      cl2<-makeCluster(NumCores-2)
      registerDoParallel(cl2)
      clusterExport(cl2,varlist=c("obsdatval","modfit1","keepVars"),envir=environment())
     clusterEvalQ(cl2, {library(glmmTMB)
                        library(cplm)    
                        library(MASS) } )
     modfit2<-try(MuMIn::pdredge(modfit1,rank=selectCriteria,fixed=keepVars,extra=extras,cl2))
     stopCluster(cl2)
   } else {
     modfit2<-try(dredge(modfit1,rank=selectCriteria,fixed=keepVars,extra=extras))
   }
     if(class(modfit2)[1]!="try-error") {
      modfit3<-get.models(modfit2,1)[[1]] 
    } else {
      modfit2<-NULL
      modfit3<-NULL
      }
    if(printOutput & !is.null(modfit2)) {
     write.csv(modfit2,paste0(dirname[[run]],common[run],catchType[run],"ModelSelection",modType,".csv"))
#     printTableFunc(paste(modType," Model Selection"),df1=mutate(modfit2,weight=round(weight,4)),filename=paste0(dirname[[run]],common[run],catchType[run],"ModelSelection",modType,".pdf"),useRowNames=TRUE)
     if(modType %in% c("Binomial","NegBin")) anova1=anova(modfit3,test="Chi")
     if(modType %in% c("Tweedie","TMBnbinom1","TMBnbinom2","TMBtweedie")) anova1=NULL
     if(modType %in% c("Lognormal","Gamma")) anova1=anova(modfit3,test="F")
     if(!is.null(anova1)) {
      write.csv(anova1,paste0(dirname[[run]],common[run],catchType[run],modType,"Anova.csv"))
#      printANOVAFunc(paste(modType,"Anova"),df1=anova1,filename=paste0(dirname[[run]],common[run],catchType[run],modType,"Anova.pdf"),useRowNames=TRUE)
     }
    }
    returnval=list(modfit3,modfit2)
  }
  returnval  
}

#Function to predict total positive trips or total catches with prediction SE
#NOte variance calculations are wrong. Did not include covariance. Use makePredictiosnSimVar instead
# makePredictionsVar<-function(modfit1,modfit2=NULL,modType,newdat,obsdatval,printOutput=FALSE,nsims=nSims) {
#   newdat$Effort=newdat$Effort/newdat$SampleUnits
#   newdat=uncount(newdat,SampleUnits)
#   nObs=dim(newdat)[1]
#   if(!is.null(modfit1)) {
#     response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
#     if(dim(response1)[2]==1) {
#       names(response1)="fit"
#       if(modType=="Tweedie")
#         response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp",offsetval=NULL, nsim=nSims) else
#           response1$se.fit=rep(NA,dim(response1)[2])
#     }
#     if(!is.null(modfit2))  {
#       response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
#       names(response2)=paste0(names(response2),"2")
#     }    
#     if(modType== "Binomial") {
#       allpred<-cbind(newdat,response1) %>%
#         mutate(Total=fit,
#                TotalVar=se.fit^2+fit*(1-fit))
#     }
#     if(modType == "Lognormal" ){
#       allpred<-cbind(newdat,response1,response2) %>% 
#         mutate(pos.cpue=lnorm.mean(fit2,sqrt(se.fit2^2+sigma(modfit2)^2)),
#                pos.cpue.se=lnorm.se(fit2,sqrt(se.fit2^2+sigma(modfit2)^2)),
#                prob.se=sqrt(se.fit^2+fit*(1-fit))) %>% 
#         mutate(Total=Effort*fit*pos.cpue,
#                TotalVar=Effort^2*lo.se(fit,prob.se,pos.cpue,pos.cpue.se)^2) 
#     }
#     if(modType == "Gamma"){
#       allpred<-cbind(newdat,response1,response2) %>%
#         mutate(pos.cpue.se=sqrt(se.fit2^2+fit2^2*gamma.shape(modfit2)[[1]]),
#                prob.se=sqrt(se.fit^2+fit*(1-fit))) %>% 
#         mutate(Total=Effort*fit*fit2,
#                TotalVar=Effort^2*lo.se(fit,prob.se,fit2,pos.cpue.se)^2) 
#     }
#     if(modType =="NegBin") {
#       allpred<-cbind(newdat,response1)   %>% 
#         mutate(Total=fit,
#                TotalVar=se.fit^2+fit+fit^2/modfit1$theta)              
#     }
#     if(modType =="Tweedie") {
#       allpred<-cbind(newdat,response1)   %>% 
#         mutate(Total=Effort*fit,
#                TotalVar=Effort^2*(se.fit^2+modfit1$phi*fit^modfit1$p)) 
#     }
#     if(modType == "TMBnbinom1") {
#       allpred<-cbind(newdat,response1)  %>% 
#         mutate(Total=fit,
#                TotalVar=se.fit^2+fit+fit*sigma(modfit1))              
#     }
#     if(modType == "TMBnbinom2") {
#       allpred<-cbind(newdat,response1)  %>% 
#         mutate(Total=fit,
#                TotalVar=se.fit^2+fit+fit^2/sigma(modfit1))              
#     }
#     if(modType == "TMBtweedie") {
#       allpred<-cbind(newdat,response1)  %>% 
#         mutate(Total=Effort*fit,
#                TotalVar=Effort^2*(se.fit^2+sigma(modfit1)*fit^(glmmTMB:::.tweedie_power(modfit1))))              
#     }
#     stratapred<-allpred %>% 
#       group_by_at(all_of(requiredVarNames)) %>% 
#       summarize(Total=sum(Total),
#                 TotalVar=sum(TotalVar))  %>% 
#       mutate(Total.se=sqrt(TotalVar))    
#     yearpred<-stratapred%>% group_by(Year) %>% 
#       summarize(Total=sum(Total,na.rm=TRUE),TotalVar=sum(TotalVar)) %>%
#       mutate(Year=as.numeric(as.character(Year)),Total.se=sqrt(TotalVar))%>%
#       mutate(Total.cv=Total.se/Total)
#     if(is.na(max(yearpred$Total.cv)) | max(yearpred$Total.cv,na.rm=TRUE)>10) {
#       print(paste(common[run],modType," CV >10 or NA variance")) 
#       yearpred<-NULL
#       stratapred<-NULL
#       allpred<-NULL
#       returnval=NULL
#     }  else {
#       returnval=yearpred
#       if(printOutput) {
#         write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modType,"StratumSummary.csv"))
#         write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modType,"AnnualSummary.csv"))
#       }
#     } 
#   } else {
#     returnval=NULL
#   }
#   returnval
# }

#Function to predict without variances to get predictions quickly for cross validation
makePredictions<-function(modfit1,modfit2=NULL,modType,newdat) {
  if(!is.null(modfit1)) {
    predval1<-try(data.frame(predict(modfit1,newdata=newdat,se.fit=TRUE,type="response")))
    if(class(predval1)[[1]]!="try-error") {
     if(!is.null(modfit2))  {
      predval2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      names(predval2)=paste0(names(predval2),"2")
     }
     if(modType=="Lognormal") {
        allpred<-cbind(newdat,predval1,predval2)  %>% 
          mutate(est.cpue=fit*lnorm.mean(fit2,sqrt(se.fit2^2+sigma(modfit2)^2))) 
     }
     if(modType=="Gamma") {
      allpred<-cbind(newdat,predval1,predval2)   %>% 
        mutate(est.cpue=fit*fit2) 
    }
    if(modType %in% c("NegBin","TMBnbinom1","TMBnbinom2")) {
      allpred<-cbind(newdat,predval1)   %>% 
        mutate(est.cpue=fit/Effort) 
    }
    if(modType=="TMBtweedie"){
      allpred<-cbind(newdat,predval1)   %>% 
        mutate(est.cpue=fit) 
    }
    if(modType =="Tweedie") {  
      allpred<-data.frame(est.cpue=predict(modfit1,newdata=newdat,type="response"))  
    }
    returnval=allpred
    } else returnval=NULL
  } else {
    returnval=NULL
  }
  returnval
}

## Function to get an abundance index with SE, Does not yet have year interactions as random effects
makeIndexVar<-function(modfit1,modfit2=NULL,modType,obsdatval,newdat=newDat,printOutput=FALSE,nsims=nSims) {
    returnval=NULL
    if(!is.null(modfit1)) {
    response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modType=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp",offsetval=NULL, nsim=nSims) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      names(response2)=paste0(names(response2),"2")
    }    
    if(modType== "Binomial") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Index=fit,
               SE=se.fit)
    }
    if(modType == "Lognormal" ){
      allpred<-cbind(newdat,response1,response2) %>% 
        mutate(pos.cpue=lnorm.mean(fit2,se.fit2),
               pos.cpue.se=lnorm.se(fit2,se.fit2),
               prob.se=se.fit) %>% 
        mutate(Index=fit*pos.cpue,
               SE=lo.se(fit,prob.se,pos.cpue,pos.cpue.se)) 
    }
    if(modType == "Gamma"){
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue.se=se.fit2,
               prob.se=se.fit) %>% 
        mutate(Index=fit*fit2,
               SE=lo.se(fit,prob.se,fit2,pos.cpue.se)) 
    }
    if(modType =="NegBin") {
      allpred<-cbind(newdat,response1)   %>% 
        mutate(Index=fit,
               SE=se.fit)              
    }
    if(modType =="Tweedie") {
      allpred<-cbind(newdat,response1)   %>% 
        mutate(Index=fit,
               SE=se.fit) 
    }
    if(modType == "TMBnbinom1") {
      allpred<-cbind(newdat,response1)  %>% 
        mutate(Index=fit,
               SE=se.fit)              
    }
    if(modType == "TMBnbinom2") {
      allpred<-cbind(newdat,response1)  %>% 
        mutate(Index=fit,
               SE=se.fit)              
    }
    if(modType == "TMBtweedie") {
      allpred<-cbind(newdat,response1)  %>% 
        mutate(Index=fit,
               SE=se.fit)              
    }
    allpred=allpred %>% mutate(ymin=Index-SE,ymax=Index+SE)  %>%
     mutate(ymin=ifelse(ymin<0,0,ymin))
    returnval=allpred
    if(printOutput) {
        write.csv(allpred,paste0(dirname[[run]],common[run],catchType[run],modType,"Index.csv"))
      }
  }
  returnval
}

#Generate standard errors and confidence intervals of predictions from simulation from regression coefficients and their var/covar matrix
makePredictionsSimVar<-function(modfit1,modfit2=NULL,newdat, modtype,  nsim=nSims, printOutput=TRUE) {
  #Separate out sample units
  newdat$Effort=newdat$Effort/newdat$SampleUnits
  newdat=uncount(newdat,SampleUnits)
  nObs=dim(newdat)[1]
  #Get predictions
  response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modtype=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp",offsetval=NULL, nsim=nSims) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      names(response2)=paste0(names(response2),"2")
    }
  #Set up model matrices for simulation
  yvar=sub( " ", " ",formula(modfit1) )[2]
  newdat<-cbind(y=rep(1,nObs),newdat)
  names(newdat)[1]=yvar
  a=model.matrix(formula(modfit1),data=newdat)
  if(!is.null(modfit2)) {
   yvar=sub( " ", " ",formula(modfit2) )[2] 
   if(! yvar %in% names(newdat)) {
    newdat<-cbind(y=rep(1,nObs),newdat)
    names(newdat)[1]=yvar
   }  
   b=model.matrix(formula(modfit2),data=newdat)
  }
  #Get predictions
  if(modtype == "Binomial") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=fit,TotalVar=se.fit^2+fit*(1-fit))
      sim=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )  
     }    
  if(modtype == "Lognormal") {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(fit2,sqrt(se.fit2^2+sigma(modfit2)^2)),
               pos.cpue.se=lnorm.se(fit2,sqrt(se.fit2^2+sigma(modfit2)^2)),
               prob.se=sqrt(se.fit^2+fit*(1-fit))) %>%
        mutate(Total=Effort*fit*pos.cpue,
               TotalVar=Effort^2*lo.se(fit,prob.se,pos.cpue,pos.cpue.se)^2)
      sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )  
      sim2=replicate(nsim,newdat$Effort*exp(rnorm(nObs,b %*% 
           mvrnorm(1,coef(modfit2),vcov(modfit2)),sigma(modfit2)) ) )  
      sim=sim1*sim2
  }
  if(modtype == "Gamma") {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue.se=sqrt(se.fit2^2+fit2^2*gamma.shape(modfit2)[[1]]),
               prob.se=sqrt(se.fit^2+fit*(1-fit))) %>%
        mutate(Total=Effort*fit*fit2,
               TotalVar=Effort^2*lo.se(fit,prob.se,fit2,pos.cpue.se)^2)
      sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )  
      sim2=replicate(nsim,newdat$Effort*simulateGammaDraw(modfit2,nObs,b) )  
      sim=sim1*sim2
  }
  if(modtype=="NegBin") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=fit,TotalVar=se.fit^2+fit+fit^2/modfit1$theta)
      sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))*newdat$Effort, 
        size=modfit1$theta))  #Simulate negative binomial data
  }
  if(modtype=="Tweedie") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=Effort*fit,
               TotalVar=Effort^2*(se.fit^2+modfit1$phi*fit^modfit1$p))
      sim=replicate(nsim,rtweedie(nObs,power=modfit1$p,
        mu=as.vector(exp(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))),phi=modfit1$phi))*newdat$Effort  
  }
  if(modtype=="TMBnbinom1") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=fit,
               TotalVar=se.fit^2+fit+fit*sigma(modfit1))
      sim = replicate(nsim,simulateNegBin1Draw(modfit1,nObs,a,newdat$Effort))  
  }
  if(modtype=="TMBnbinom2") {
       allpred<-cbind(newdat,response1)  %>%
        mutate(Total=fit,
               TotalVar=se.fit^2+fit+fit^2/sigma(modfit1))
      sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,fixef(modfit1)[[1]],
        vcov(modfit1)[[1]]))*newdat$Effort, size=sigma(modfit1)))  #Simulate negative binomial data
  }
  if(modtype=="TMBtweedie") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=Effort*fit,
               TotalVar=Effort^2*(se.fit^2+sigma(modfit1)*fit^(glmmTMB:::.tweedie_power(modfit1))))
      sim=replicate(nsim,rtweedie(nObs,power=glmmTMB:::.tweedie_power(modfit1),
        mu=as.vector(exp(a %*% mvrnorm(1,fixef(modfit1)[[1]],vcov(modfit1)[[1]]))),
        phi=sigma(modfit1)))*newdat$Effort  
  }
  stratatotal<-allpred %>%
      group_by_at(all_of(requiredVarNames)) %>%
      summarize(Total=sum(Total,na.rm=TRUE))
  yeartotal<-allpred%>% group_by(Year) %>%
      summarize(Total=sum(Total,na.rm=TRUE))
  stratapred<-cbind(newdat,sim) %>% 
       group_by_at(all_of(requiredVarNames)) %>% 
       summarize_at(.vars=as.character(1:nsim),.fun=sum,na.rm=TRUE) %>%
       rowwise() %>%
       mutate(Total.mean=mean(c_across(as.character(1:nsim))),
         TotalVar=var(c_across(as.character(1:nsim))),
         TotalLCI=quantile(c_across(as.character(1:nsim)),p=CIval/2),
         TotalUCI=quantile(c_across(as.character(1:nsim)),p=1-CIval/2)) %>% 
       mutate(Total.se=sqrt(TotalVar))  %>%
       mutate(Total.cv=Total.se/Total.mean)  %>%
        dplyr::select(-one_of(as.character(1:nsim)))
  stratapred$Total=stratatotal$Total
  yearpred<-cbind(newdat,sim) %>% 
       group_by(Year) %>% 
       summarize_at(.vars=as.character(1:nsim),.fun=sum,na.rm=TRUE) %>%
       rowwise() %>%
       mutate(Total.mean=mean(c_across(as.character(1:nsim))),
         TotalVar=var(c_across(as.character(1:nsim))),
         TotalLCI=quantile(c_across(as.character(1:nsim)),p=CIval/2),
         TotalUCI=quantile(c_across(as.character(1:nsim)),p=1-CIval/2)) %>% 
       mutate(Total.se=sqrt(TotalVar))  %>%
       mutate(Total.cv=Total.se/Total.mean)  %>%
        dplyr::select(-one_of(as.character(1:nsim)))
  yearpred$Total<-yeartotal$Total 
  if(is.na(max(yearpred$Total.cv)) | max(yearpred$Total.cv,na.rm=TRUE)>10) {
       print(paste(common[run],modtype," CV >10 or NA variance")) 
       yearpred<-NULL
       stratapred<-NULL
       returnval=NULL
  }  else  {     returnval=yearpred  }
  if(printOutput) {
       write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modtype,"StratumSummary.csv"))
       write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modtype,"AnnualSummary.csv"))
  }
  returnval
}

simulateGammaDraw<-function(modfit,nObs,b) {
    muval<-exp(b %*% mvrnorm(1,coef(modfit),vcov(modfit)))
    shapeval<-gamma.shape(modfit)[[1]]
    scaleval<-muval/shapeval
    rgamma(nObs,shape=shapeval,scale=scaleval) 
}
simulateNegBin1Draw<-function(modfit,nObs,b,Effort) {
    muval<-exp(b %*% mvrnorm(1,fixef(modfit)[[1]],vcov(modfit)[[1]]))*Effort
    thetaval<-sigma(modfit)
    rnbinom(nObs,mu=muval,size=muval/thetaval) 
}


#Function to plot either total positive trips (binomial) or total catch/bycatch (all other models)
plotSums<-function(yearpred,modType,fileName,subtext="") {
  if(!is.null(yearpred)) {
    if(modType=="Binomial") ytitle=paste0(common[run]," ","predicted total positive trips") else
      ytitle=paste0("Total",common[run]," ",catchType[run]," (",catchUnit[run],")")
    if(modType %in% c("Lognormal","Gamma")) modType=paste("Delta",modType)
    yearpred<-yearpred %>% 
      mutate(Year=as.numeric(as.character(Year)),ymin=Total-Total.se,ymax=Total+Total.se) %>%
      mutate(ymin=ifelse(ymin>0,ymin,0))
    if(modType=="All") {
      g<-ggplot(yearpred,aes(x=Year,y=Total,ymin=TotalLCI,ymax=TotalUCI,fill=Source))+
        geom_line(aes(col=Source))+ geom_ribbon(alpha=0.3)+xlab("Year")+
#        geom_line(aes(y=Total.mean,col=Source),lty=2,lwd=2)+
        ylab(ytitle)
    } else {
      g<-ggplot(yearpred,aes(x=Year,y=Total,ymin=TotalLCI,ymax=TotalUCI))+
        geom_line(aes(y=Total.mean),lty=2)+
        geom_line()+ geom_ribbon(alpha=0.3)+xlab("Year")+
        ylab(ytitle)
    }
    print(g)
    if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
}
}

#Function to plot abundance index plus minus standard error
plotIndex<-function(yearpred,modType,fileName,subtext="") {
  if(!is.null(yearpred)) {
    if(modType=="Binomial") ytitle=paste0(common[run]," ","Positive trip index") else
      ytitle=paste0("Index ", common[run]," ",catchType[run]," (",catchUnit[run],")")
    if(modType %in% c("Lognormal","Gamma")) modType=paste("Delta",modType)
    yearpred<-yearpred %>% mutate(Year=as.numeric(as.character(Year)),ymin=Index-SE,ymax=Index+SE) %>%
      mutate(ymin=ifelse(ymin>0,ymin,0))
    if(modType=="All") {
      g<-ggplot(yearpred,aes(x=Year,y=Index,ymin=ymin,ymax=ymax,fill=Source))+
        geom_line(aes(col=Source))+ geom_ribbon(alpha=0.3)+xlab("Year")+
        ylab(ytitle)
    } else {
      g<-ggplot(yearpred,aes(x=Year,y=Index,ymin=ymin,ymax=ymax))+
        geom_line()+ geom_ribbon(alpha=0.3)+xlab("Year")+
        ylab(ytitle)
    }
    if(length(indexVarNames)>1) {
        varplot=as.formula(paste0("~",paste(grep("Year",indexVarNames,invert=TRUE,value=TRUE),sep="+")))
        g=g+facet_wrap(varplot)
    }
    print(g)
    if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
  }
}

#Function plots residuals with both R and Dharma library and calculate residual diagnostics.
ResidualsFunc<-function(modfit1,modType,fileName,nsim=250) {
  require(quantreg)
  if(!is.null(fileName))   pdf(fileName,height=5,width=7)
  if(!is.null(modfit1)) {
    dfcheck<-data.frame(Expected=predict(modfit1),Residuals=residuals(modfit1))
    g1<-ggplot(dfcheck,aes(x=Expected,y=Residuals))+geom_point()+
        geom_abline(intercept=0,slope=0)+
        ggtitle(paste("a. ",modType,"ordinary residuals"))
    g2<-ggplot(dfcheck,aes(sample=Residuals))+geom_qq()+geom_qq_line()+
        ggtitle(paste("b. QQ normal of residuals"))
    if(class(modfit1)[1] =="cpglm") {  #Extra step to simulate DHARMa for cpglm or mgcv
      simvals=simulateTweedie(modfit1,nsim)
      simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                          fittedPredictedResponse = predict(modfit1,type="response")))
      if(class(simulationOutput)[1]=="try-error") {
        simvals=simulateTweedie(modfit1,nsim*4)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                            fittedPredictedResponse = predict(modfit1,type="response")))
      }
      if(class(simulationOutput)[1]=="try-error") {
        simvals=simulateTweedie(modfit1,nsim*10)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                            fittedPredictedResponse = predict(modfit1,type="response")))
      }
    }  
    if(class(modfit1)[1]=="gam" & modType=="NegBin") {
      simvals=simulateNegBinGam(modfit1,nsim)
      simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                          fittedPredictedResponse = predict(modfit1,type="response")))
      if(class(simulationOutput)[1]=="try-error") {
        simvals=simulateNegBinGam(modfit1,nsim*4)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                            fittedPredictedResponse = predict(modfit1,type="response")))
      }
      if(class(simulationOutput)[1]=="try-error") {
        simvals=simulateNegBinGam(modfit1,nsim*10)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                            fittedPredictedResponse = predict(modfit1,type="response")))
      }
      if(class(simulationOutput)[1]=="try-error") {
        simvals=simulateNegBinGam(modfit1,nsim*15)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y , 
                                            fittedPredictedResponse = predict(modfit1,type="response")))
      }
    }  
    if( class(modfit1)[1] !="cpglm" & !(class(modfit1)[1]=="gam" & modType=="NegBin"))     {  #Regular DHARMa residuals work for everything but cpglm
      simulationOutput <- try(simulateResiduals(fittedModel = modfit1, n = nsim))
      if(class(simulationOutput)[1]=="try-error") simulationOutput <- try(simulateResiduals(fittedModel = modfit1, n = nsim*4))
      if(class(simulationOutput)[1]=="try-error") simulationOutput <- try(simulateResiduals(fittedModel = modfit1, n = nsim*10))
    }
    if(class(simulationOutput)[1]!="try-error") {
#      plot(simulationOutput, quantreg = F)
#      title(modType,outer=2,line=-1)
      df1<-data.frame(Residual=simulationOutput$scaledResiduals,Predictor=simulationOutput$fittedPredictedResponse) %>%
        arrange(Residual) %>%
        mutate(Empirical.Quantile=(1:n())/(n()+1)) %>%
        mutate(Expected.Quantile=qunif(Empirical.Quantile)) %>%
        mutate(Rank.Predictor= rank(Predictor, ties.method = "average")) %>%
        mutate(Rank.Predictor = Rank.Predictor/max(Rank.Predictor))
      g3<-ggplot(df1,aes(x=Expected.Quantile,y=Residual))+geom_point()+
        geom_abline()+ylab("DHARMa scaled residuals")+ xlab("Expected quantile")+
        ggtitle("c. QQ uniform scaled residuals")
      g4<-ggplot(df1,aes(x=Rank.Predictor,y=Residual))+
        geom_point()+xlab("Model predictions (rank transformed)")+
        ylab("DHARMa scaled residuals")+ggtitle("d. Scaled residual vs. predicted")+
        geom_hline(aes(yintercept=0.5),lty=2)+geom_hline(aes(yintercept=0.75),lty=2)+geom_hline(aes(yintercept=0.25),lty=2)+
        geom_quantile(method = "rqss",col="red", formula=y ~ qss(x, lambda = 2))
      grid.arrange(g1,g2,g3,g4,ncol=2)
      test1=testUniformity(simulationOutput,plot=FALSE)
      test2=testDispersion(simulationOutput,plot=FALSE)
      test3=testZeroInflation(simulationOutput,plot=FALSE)
      test4=testOutliers(simulationOutput,plot=FALSE)
      returnval=c(test1$statistic,
                  test1$p.value,
                  test2$statistic,
                  test2$p.value,
                  test3$statistic,
                  test3$p.value,
                  test4$statistic,
                  test4$p.value)
      names(returnval)=c("KS.D","KS.p","Dispersion.ratio","Dispersion.p","ZeroInf.ratio","ZeroInf.p","Outlier","Outlier.p")
    } else returnval=NULL    
  } else returnval=NULL
 if(!is.null(fileName))  dev.off()
  returnval
}

#Function to plot total catch by all models plus a validation number 
plotSumsValidate<-function(yearpred,trueval,fileName,colName) {
  yearpred<-yearpred %>% mutate(ymin=Total-Total.se,ymax=Total+Total.se)%>%
    mutate(ymin=ifelse(ymin>0,ymin,0)) %>%
    dplyr::filter(Valid==1)
  trueval<-trueval %>% rename(Total=!!colName) %>%
    mutate(ymin=rep(NA,dim(trueval)[1]),ymax=rep(NA,dim(trueval)[1]),
                  Source=rep("Validation",dim(trueval)[1]))
  yearpred<-rbind(yearpred[,c("Year","Total","ymin","ymax","Source")],trueval[,c("Year","Total","ymin","ymax","Source")])
  g<-ggplot(yearpred,aes(x=Year,y=Total,ymin=ymin,ymax=ymax,color=Source,fill=Source))+
      geom_line()+ geom_ribbon(alpha=0.3)+ylab("Total catch (kg)")+xlab("Year")+
      ylab(paste0(common[run]," ",catchType[run]," (",catchUnit[run],")"))+
      geom_point(data=yearpred[yearpred$Source=="Validation",],aes(x=Year,y=Total),size=2)
  print(g)
  if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
 }

#Function to plot boxplots of RMSE and ME across folds
plotCrossVal<-function(rmse,me,fileName) {
 RMSE<-pivot_longer(rmse,everything(),names_to="Model")
 ME<-pivot_longer(me,everything(),names_to="Model")
 df<-bind_rows(list(RMSE=RMSE,ME=ME),.id="Metric")
 g<-ggplot(df)+geom_boxplot(aes(x=Model,y=value),fill="lightgrey")+
   facet_wrap(Metric~.,ncol=1,scales="free")+
   xlab("Model")+ylab("Cross validation metrics")
 print(g)
 if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
}

#Calculate RMSE
getRMSE<-function(yhat,y) {
  if(!is.null(yhat))
    rmse=sqrt(sum((yhat-y)^2)/length(y)) else
  rmse=NA
  rmse
}

#Calculate mean error
getME<-function(yhat,y) {
  if(!is.null(yhat))
    me=sum(yhat-y)/length(y) else
      me=NA
    me
}

#Function to fit a specified model formula and print outputs
FitModelFunc<-function(formula1,formula2,modType,obsdatval,outputDir) {
  modfit2=NULL
  formula3=update(formula2,~.+offset(log(Effort)))  
  if(modType %in% c("Binomial","Lognormal","Gamma") )  {  
    obsdatval$y=obsdatval$pres
    modfit1<-try(glm(formula1,data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=1000)))
  }
  if(modType=="Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit2=try(lm(formula2,data=obsdatval))
  }
  if(modType=="Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit2=try(glm(formula2,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="NegBin") {
    obsdatval$y=round(obsdatval$Catch)
    modfit1=try(glm.nb(formula3,data=obsdatval,control=glm.control(epsilon=1E-6,maxit=45),na.action=na.fail))
  }
  if(modType=="Tweedie") {
    obsdatval$y=obsdatval$cpue
    modfit1=try(cpglm(formula2,data=obsdatval))
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    obsdatval$y=round(obsdatval$Catch)
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula3,family=TMBfamily,data=obsdatval))
  }
  if(modType =="TMBtweedie"){
    obsdatval$y=obsdatval$cpue
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula2,family=TMBfamily,data=obsdatval))
  }
  if(class(modfit1)[1]=="try-error") modfit1=NULL
  if(class(modfit2)[1]=="try-error") modfit2=NULL
  if(!is.null(modfit1)) {
    if(modType %in% c("Binomial","Lognormal","Gamma"))  #for delta models write binomial anova
      write.csv(anova(modfit1,test="Chi"),file=paste0(outputDir,"/BinomialAnova.csv"))
    if(modType %in% c("NegBin")) anova1=anova(modfit1,test="Chi")
    if(modType %in% c("Tweedie","TMBtweedie","TMBnbinom1","TMBnbinom2")) anova1=NULL
    if(modType %in% c("Lognormal","Gamma")) anova1=anova(modfit2,test="F")
    if(!is.null(anova1)) {
      write.csv(anova1,paste0(outputDir,"/anova",modType,".csv"))
    }
    }
  list(modfit1=modfit1,modfit2=modfit2)
}

#Function to fit a specified model formula and print outputs for Cross validation
FitModelFuncCV<-function(formula1,modType,obsdatval) {
  if(modType %in% c("Binomial") )  {  
    obsdatval$y=obsdatval$pres
    modfit1<-try(glm(formula1,data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=1000)))
  }
  if(modType=="Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(lm(formula1,data=obsdatval))
  }
  if(modType=="Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(glm(formula1,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="NegBin") {
    obsdatval$y=round(obsdatval$Catch)
    modfit1=try(glm.nb(formula1,data=obsdatval,control=glm.control(epsilon=1E-6,maxit=45),na.action=na.fail))
  }
  if(modType=="Tweedie") {
    obsdatval$y=obsdatval$cpue
    modfit1=try(cpglm(formula1,data=obsdatval))
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    obsdatval$y=round(obsdatval$Catch)
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula1,family=TMBfamily,data=obsdatval))
  }
  if(modType =="TMBtweedie"){
    obsdatval$y=obsdatval$cpue
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula1,family=TMBfamily,data=obsdatval))
  }
  if(class(modfit1)[1]=="try-error") modfit1=NULL
  modfit1
}

#Function to simulate DHARMa residuals from a negative binomial GAM or GLM
simulateNegBinGam <- function(modfit, nsims=250, offsetval=1){
  muval = predict(modfit, type = "response")*offsetval  #Get the mean with offset
  nObs = length(muval)  
  thetaval = modfit$family$getTheta(trans=TRUE)  #Get theta not log transforme
  sim = replicate(nsims,rnbinom(nObs,mu=muval, size=thetaval))  #Simulate negative binomial data
  sim
}

