ts=proc.time()[3]
options(stringsAsFactors=F)
args<-commandArgs(trailingOnly=TRUE)
load(args[1])

library(dplyr)
library(caret)
library(doParallel)
library(scales)

source("init_threshSet.R")
source("seesaw_determine_thresh.R")
source("seesaw_data_analysis.R")
cat('time used: ', proc.time()[3]-ts, '\n')


## Control variables
thresh_pStat=0.01

## Minimum amount of cells deemed positive
cellThresh=20
registerDoParallel(cores=10)

all_gene_imp=list()
run.date<-Sys.Date()
batch_size=12
splits=ceiling(length(topWithinThresh)/batch_size)
#split the dependencies into batches
###IF BATCH ALREADY EXISTS, i.e, OLD BATCH SET IS LOADED IN, 
cat('new stuff')
if(!exists("batch")){
cat('making new batch')
batch<-createFolds(topWithinThresh, k=splits)}
if(exists("fold")){
cat('here is a fold')
start=fold
rm(fold)
}else
{
start=1
cat('no folds')
}
cat(start)
## For each gene dependency
for (fold in start:length(batch)){
workload<-topWithinThresh[batch[[fold]]]
cat(fold)
for (depGene in names(workload)){
 cat(depGene)
  ## Gene importance for one gene
  gene_imp=data.frame()
  ## For each gene dependency-expression pair
  result=foreach(expGene = colnames(m_filtered), .combine=cbind, .multicombine=TRUE) %dopar%{
    # Filters the NAs in dependencies and expressions and also scales expressions
    depExp=filt_dep_exp(score_filtered[,depGene],m_filtered[,expGene])
    dp=depExp[[1]]
    ex=depExp[[2]]
    
    ## Determining boundary between positive and negative cell lines
    ## Data frame containing BS (seesaw_determine_thresh.R) for multiple thresholds
    boundList=find_allBound(dp,ex,cellThresh)
    ## Best boundary as determined 
    bound=find_optimal(boundList)
    pred<-dp<=bound[1]
    coefs<-glm(pred~ex, family='binomial')[[1]]
    

    
    ## combines results by generating seeScore,sawScore,coeff,pStat and appending to bound
    ret=c(generate_quant_sep(dp,ex,bound),generate_coeff(dp,ex,bound),bound[1],bound[2],coefs[1],coefs[2])
  }
  cat('end of foreach time: ', proc.time()[3]-ts, '\n')
  ## Creates a dataframe gene_importance that contains 5 parameters:
  ## Seescore, Sawscore, Coefficient, pStat, Bound, Brier Score
  gene_imp=data.frame("seeScore"=matrix(unlist(result[1,])),
                      "sawScore"=matrix(unlist(result[2,])),
                      "coefficient"=matrix(unlist(result[3,])),
                      "pStat"=matrix(unlist(result[4,])),
                      "bound"=matrix(unlist(result[5,])),
			"BS"=matrix(unlist(result[6,])),
			"b0"=matrix(unlist(result[7,])),
			"b1"=matrix(unlist(result[8,])))
  rownames(gene_imp)=colnames(m_filtered)
  ## Filters gene_imp for multiple things
  ## See documentation on seesaw_data_analysis
  gene_imp_filt=filt_gene_imp(gene_imp,thresh_pStat) 

  ## Adds them all to a big list
  all_gene_imp[[depGene]]=gene_imp_filt
  save(batch,all_gene_imp,fold,file=paste('mainframe_results_',run.date,'.RData',sep="")) #will help keep track of results after each dependency is run, ss_vectors may not be built but that's no issue;the spaced formatting indicates stopped in middle of a fold
}
 ##construct vectors
 ss_vect<-build_SeeSaw_vect(m_filtered, all_gene_imp)
cat('time used: ', proc.time()[3]-ts, '\n')

#save(all_gene_imp,ss_vect,file=paste('mainframe_results_11.1.19',fold,Sys.Date(),'.RData'))
#saves the data split and updates the data after moving through each fold so problems can be diagnosed, NOTE REDUNDANCIES, TAKE THIS INTO ACCOUNT WHEN ANALYZING
fold<-fold+1
save(batch,all_gene_imp,ss_vect,fold,file=paste('mainframe_results_',run.date,'.RData',sep=""))
fold<-fold-1
cat(paste('Finished batch',fold))
}
