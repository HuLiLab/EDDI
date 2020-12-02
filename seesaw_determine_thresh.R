options(stringsAsFactors=F)
library(dplyr)
library(randomForest)
library(caret)
library(doParallel)

# Static control variables
numFolds=10

## Input: Dependency and expression for a gene dep exp pair,
##        bound considered to generate Brier score
## numFolds folds of cross validated classifier models (logistic regression)
## Aggregates Brier score across each fold to generate an average BS across folds
## Return: Average Brier score for the considered bound
find_BS=function(dep, exp, bound){
  dCells=ifelse(dep<=bound,yes=1,no=0)
  foldI=createFolds(as.factor(dCells), k=numFolds, list=FALSE)
  MeanBS=c()
  for (n in 1:numFolds){
    d_test=dCells[foldI==n]
    d_train=dCells[foldI!=n]
    exp_test=exp[foldI==n]
    exp_train=exp[foldI!=n]
    dataTrain=data.frame(d=d_train, exp=exp_train)
    classify=glm(d~exp, family=binomial,data=dataTrain)
    predictions<-predict(classify,data.frame(exp=exp_test),type='response')
    MeanBS=c(MeanBS,mean((predictions-d_test)^2))
  }
  mean(MeanBS)
}

## Input: Dependency and expression for a gene dep exp pair,
##        minimal amount of positive cell lines possible
## Increments by 0.1 along range of bounds from min amount to 0
## Generates a Brier score for each bound
## Return: Data frame containing all considered bounds and Brier scores for each
find_allBound=function(dep, exp, cellThresh){
  depThresh=dep[order(dep)][1:cellThresh][cellThresh]
  BSlist=data.frame()
  for (bound in seq(0,depThresh,by=-0.1)){
    BS=find_BS(dep, exp, bound)
    BSlist=rbind(BSlist,data.frame("BS"=BS,
                                     row.names=c(bound)))
  }
  BSlist
}
## Input: Data frame containing bounds as rownames and
##        a column BS to measure the Brier Score of each bound
## Uses the change in Brier score to determine improvement of the classifier
## After a certain threshold, the classifier will cease to improve
## The threshold is currently set as the median of the improvement function
## Return: Optimal bound in decimal,BSscores alongside improvement function, improvement boundary
find_optimal_fullstat=function(BSlist){
  BSordered=BSlist[order(-as.numeric(rownames(BSlist))),,drop=FALSE]
  improve=with(BSordered, ave(BS,FUN=function(x) c(0,diff(x))))
  BS_improve=cbind(BSordered,improve)
  ## improveBound is how the function knows at which point the function ceases to improve
  improveBound=quantile(BS_improve$improve)[3]
  for (bound in rownames(BS_improve)){
    if (BS_improve[bound,"improve"]<0 & BS_improve[bound,"improve"]>improveBound){
      return(list(as.numeric(bound),BS_improve,improveBound))
    }
  }
}

find_optimal=function(BSlist){
  BSordered=BSlist[order(-as.numeric(rownames(BSlist))),,drop=FALSE]
  improve=with(BSordered, ave(BS,FUN=function(x) c(0,diff(x))))
  BS_improve=cbind(BSordered,improve)
  ## improveBound is how the function knows at which point the function ceases to improve
  improveBound=quantile(BS_improve$improve)[3]
  for (bound in rownames(BS_improve)){
    if (BS_improve[bound,"improve"]<0 & BS_improve[bound,"improve"]>improveBound){
      return(as.numeric(c(bound,BS_improve[bound,"BS"]))) #includes the optimal brier score and the the bound
    }
  }
}

## Input: Dependency and expression for a gene dep exp pair,
##        data frame containing bounds as rownames and
##        a column BS to measure the Brier Score of each bound
## Plots all the bounds with respect to their BS score
## Colorizes the best regarded bound in magenta
plot_optimal=function(dep,exp,BSlist){
  plot(dep,exp)
  for (bound in rownames(BSlist)){
    pig=(BSlist[bound,1]-range(BSlist[,1])[1])/(range(BSlist[,1])[2]-range(BSlist[,1])[1])
    abline(v=bound,col=rgb(0,1-pig,0,1))
  }
  fI=find_optimal(BSlist)
  abline(v=fI,col="magenta")
}

plot_improve=function(BSlist){
  stat=find_optimal_fullstat(BSlist)
  plot(x=rownames(stat[[2]]),y=stat[[2]]$improve)
  abline(v=stat[[1]],col="magenta")
}

#--------------------------------------------------------TEST FUNCTIONS--------------------------------------------------

#test_find_AllBound=function(){
 # source("init_dataLoad.R")
 # tempCellDep="TXNDC8"
 # tempCellExp="RPL23AP21"
 # dep=score_filtered[,tempCellDep]
 # exp=m_filtered[,tempCellExp]
 # cellThresh=20
 # BSlist=find_allBound(dep,exp,cellThresh)
 # plot_optimal(dep,exp,BSlist)
 # plot_improve(BSlist)
#}
