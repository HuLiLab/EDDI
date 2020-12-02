## Input: dependency and expression for a gene dep exp pair, and optimal pos neg bound
## Takes the 25th quantile of Pos See Cells and 75th quantile of Neg See Cells
## Takes the 75th quantile of Pos Saw Cells and 25th quantile of Neg Saw Cells
## Return: List of 2 values representing degree of separation, positive is more separated
#this ought to be standardized if it isn't
generate_quant_sep=function(dep,exp,bound){
  quantPos=quantile(exp[dep<=thresh])
  quantNeg=quantile(exp[dep>thresh])
  seeScore=quantPos[2]-quantNeg[4]
  sawScore=-(quantPos[4]-quantNeg[2])
  
  list(seeScore,sawScore)
}

## Input: dependency and expression for a gene dep exp pair, and optimal pos neg bound
## Creates a linear regression for the positive cell lines and generates pStat and coeff
## Return: List of 2 values: coefficients, pStat
generate_coeff=function(dep,exp,bound){
  linear=lm(exp[dep<=bound]~dep[dep<=bound])
  dataM=summary(linear)$coefficients
  coef=dataM[2,1]
  pStat=dataM[2,4]
  
  list(coef,pStat)
}

###------------------------------FILTERING----------------------------------

## Input: a database with 6 parameters in order: seeScore, sawScore, coefficient, pStat, bound, brier score
## Filters through p values, see scores, saw scores
## Return: a filtered database with the same 5 parameters
filt_gene_imp=function(gene_imp,thresh_pStat){
  gene_imp=gene_imp[(gene_imp[,"seeScore"]>0|gene_imp[,"sawScore"]>0),,drop=FALSE]
  gene_imp=gene_imp[gene_imp[,"pStat"]<thresh_pStat,,drop=FALSE]
  return(gene_imp)
}

##Input: full gene expression database, filtered gene importance data packaged with Seesaw scores, linear regression slope and pstat,
##dependency boundary and brier score, and logistic regression coefficients on the optimal bounds
##By dependency, each cell line has a vector of seesaw probabilities of each gene
##Returns: a list of dataframes for each dependency; rows are cell lines, columns are genes, elements are see-saw probability
build_SeeSaw_vect<-function(expressions,geneimp){
	#IMPLEMENT BS threshold
  output<-foreach(dep=names(geneimp)) %dopar%{
    data=geneimp[[dep]]
    cl= rownames(expressions)
    #maybe a brier filtering step; parameters refer to seesaw score cutoffs
    b1=data[,'b1'] #vector of b1 for each seesaw gene
    b0=data[,'b0'] #vector of b0 for each seesaw gene
    #now generate expression val for all see-saw genes =>big matrix
    exp=expressions[,rownames(data)]
    #double check logistic regression forumla
    ex=1/(exp(-t(t(exp)*b1+b0))+1)
    good_saw<-data[,'sawScore']>0
    ex[,good_saw]=ex[,good_saw]*(-1)
    ex
  }
    names(output)=names(geneimp)
    output
}
## Input: dependency and expression data for one gene dep-exp pair
## Filters the NAs in the data, also scales the expression values
## Return: list with two elements: filtered dep,filtered and scaled exp
filt_dep_exp=function(dep,exp){
  exNA=is.na(exp)
  ex=exp[!exNA]
  dp=dep[!exNA]
  dpNA=is.na(dp)
  ex=ex[!dpNA]
  dp=dp[!dpNA]
  
  if (isTRUE(all.equal(sd(ex),0))){
    ex[ex!=0]=0
  }
  else{
    ex=scale(ex, center=TRUE, scale=TRUE)
  }
  
  list(dp,ex)
}
