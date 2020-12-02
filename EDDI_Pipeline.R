#assumes I've run init_dataLoad, this is the final pipeline
library(caret)
library(e1071)
library(randomForest)
library(pROC)

source('seesaw_result_analysis.R')
protoClusterEval<-function(depID,dep){
  genes<-names(Cluster2D[[dep]])
  scores<-score_filtered[names(Cluster2D[[dep]][[1]]),dep]
  depGroup<-scores[depID]
  nonDepGroup<-scores[-depID]
  
  
  #general details
  template<-rep(0,length(scores))
  template[depID]<-1
  expInfo<-sapply(genes,function(gene){
    depExp<-m_filtered[names(scores)[depID],gene]
    ndExp<-m_filtered[names(scores)[-depID],gene]
    dquant<-unname(quantile(depExp));ndquant<-unname(quantile(ndExp))
    dIQR=IQR(depExp);dMed<-median(depExp)
    ndIQR=IQR(ndExp);ndMed<-median(ndExp)
    
    fit<-cor(template,m_filtered[names(scores),gene])
    eVar<-var(depExp)
    eRange<-c(lowDepExp=max(dquant[2]-dIQR,min(depExp)),
     highDepExp=dquant[4]+dIQR) 
    neRange<-c(lowNDExp=max(ndquant[2]-ndIQR,min(ndExp)),
     highNDExp=ndquant[4]+ndIQR)
    correlation=cor(depGroup,m_filtered[names(scores)[depID],gene])
    c(templateScore=fit,variance=eVar,correlation=correlation,eRange,neRange)
  })
  
  #obtain score infor
  scoreMax<-max(depGroup)
  scoreVar<-var(depGroup)
  rbind(expInfo,scoreMax=scoreMax,scoreVar=scoreVar)
  
}

f1Models<-function(dep,features=NULL,trueCL=NULL){
  if(is.null(features)){
    features<-dependencies[dependencies$dependency==dep,2]
  }
  featureDat<-m_filtered[,features,drop=FALSE]
  labels<-as.numeric(names(score_filtered[,dep])%in%trueCL)
  #for debugging, assign the labels list to Labels
  keep<-1:length(labels)
  if(anyNA(labels)){keep<-which(!is.na(labels))}
  data<-data.frame(labels=labels[keep],featureDat[keep,,drop=FALSE])
  
  F1<-vector()
  AUC<-vector()
  for(i in 1:5){
  folds<-createFolds(as.factor(data$labels),k=10)
  for(fold in folds){
    ktrain<-data[-fold,,drop=FALSE]
    ktrain<-upSample(ktrain[,-1,drop=FALSE],as.factor(ktrain[,1]),yname="label")
    ktest<-data[fold,,drop=FALSE]
    tune<-tuneRF(x=ktrain[,-ncol(ktrain),drop=FALSE],y=as.factor(ktrain[,ncol(ktrain)]),ntreeTry = 500,trace=FALSE,plot=FALSE) # get a sense of grid size
    krf<-randomForest(x=ktrain[,-ncol(ktrain),drop=FALSE],y=as.factor(ktrain[,ncol(ktrain)]),mtry=tune[which.min(tune[,2]),1])
    predict<-predict(krf,newdata=ktest[,-1,drop=FALSE])
    confusMat<-table(ktest[,1],predict)
    precision<-confusMat[2,2]/sum(confusMat[,2])
    recall<-confusMat[2,2]/sum(confusMat[2,])
    roc=roc(response=as.numeric(ktest[,1]), predictor=as.numeric(predict), auc=T)$auc 
    F1<-c(F1,2*(precision*recall)/(precision+recall))
    AUC<-c(AUC,roc)
  }
  }
  F1[which(is.na(F1))]<-0
  tune<-tuneRF(data[,-1,drop=FALSE],as.factor(data$labels),ntreeTry = 500,trace=FALSE,plot=FALSE) # get a sense of grid size
  rf<-randomForest(x=data[,-1,drop=FALSE],y=as.factor(data$labels),
                   mtry=tune[which.min(tune[,2]),1],importance=TRUE)
  importance<-varImpPlot(rf)
  varImpPlot(rf)
  
  
  list(cvF1=c(mean(F1),sd(F1)),varImp=importance,cvAUC=c(mean(AUC),sd(AUC)),cv=AUC,rf)
}

##----------------------------Generating Clusters-------------------------------
#cluster each pair based on dep-exp
Cluster2D<-vector(mode='list',length=length(unique(dependencies$dependency)))
names(Cluster2D)<-unique(dependencies$dependency)

for(i in 1:nrow(dependencies)){
  dep<-dependencies$dependency[i]
  exp<-dependencies$seesaw.gene[i]
  d<-na.omit(score_filtered[,dep])
  e<-m_filtered[names(d),exp]
  temp<-cbind(d,e)
  clusters<-kmeans(temp,centers=2)$cluster
  Cluster2D[[dep]][[exp]]<-clusters
}

#union
Union2D<-vector(mode='list',length=length(unique(dependencies$dependency)))
names(Union2D)<-unique(dependencies$dependency)
for(dep in names(Cluster2D)){
  depCluster<-lapply(Cluster2D[[dep]],function(x){
    cl<-names(x)
    scores<-na.omit(score_filtered[cl,dep])
    clusta<-cl[x==1];clustb<-cl[x==2]
    true<-which.max(c(sum(scores[clusta]<=-2),sum(scores[clustb]<=-2)))
    cl<-cl[x==true]
    common<-intersect(cl,names(x)[scores<=-2])
    (names(x)%in% common)+1
  })
  depClusterNatural<-lapply(Cluster2D[[dep]],function(x){
    cl<-names(x)
    scores<-na.omit(score_filtered[cl,dep])
    clusta<-cl[x==1];clustb<-cl[x==2]
    true<-which.max(c(sum(scores[clusta]<=-2),sum(scores[clustb]<=-2)))
    cl<-cl[x==true]
    (names(x)%in% cl)+1
  })
  unionPrep<-lapply(depCluster,function(x){which(x==2)})
  unionPrepNat<-lapply(depClusterNatural,function(x){which(x==2)})
  
  Union2D[[dep]][['-2 cut']]<-unique(do.call('c',unionPrep))
  Union2D[[dep]][['natural']]<-unique(do.call('c',unionPrepNat))
  
}
#run eval
options(warn=0)
eval2Dn2<-list()
for(dep in names(Union2D)){
  eval2Dn2[[dep]]<-protoClusterEval(Union2D[[dep]][[1]],dep)
}

labels<-lapply(names(Union2D),function(x){
  names(Cluster2D[[x]][[1]])[Union2D[[x]][[1]]]
})
names(labels)<-names(Union2D)
save(labels,file='labels.RData')


##--------------------------------Build Classification Models--------------------------
load('labels.RData')
set.seed(1)
F1res<-list()
tryDep<-depCount$dependency
for(name in tryDep[-74]){
r<-NULL
while(is.null(r)){
  print(paste('attempting',name))  
  try(
  r<-f1Models(name,trueCL = labels[[name]])
  )
}
F1res[[name]]<-r
}
auc<-sapply(F1res,function(x){x[[3]][1]})
hist(auc)
sum(auc>0.75) 
auc[auc>0.7]


#compile final list of well-performing models
finalModList<-F1res[names(auc)[auc>0.7]]
save(finalModList,file='models.RData') 


##-------------------Differential Expression-----------------------------------------
source('DESeq.R')
resDESeq2fold<-lapply(names(resDESeqFull),function(x){
  dat<-resDESeqFull[[x]]
  diff<-which(abs(dat$log2FoldChange)>=1.0 & dat$padj<=0.05)
  dat<-dat[diff,]
  dat<-dat[order(dat$padj,decreasing=FALSE),]
})
names(resDESeq2fold)<-names(resDESeqFull)


resDESeqScreen<-lapply(names(resDESeqFull),function(x){
  dat<-resDESeq2fold[[x]]
  keep<-intersect(rownames(dat),colnames(eval2Dn2[[x]]))
  dat<-as.data.frame(dat[keep,])
})
names(resDESeqScreen)<-names(resDESeqFull)
##-----------------------------RF Feature Classifications------------------------
fastplot<-function(dep,set){
  exp<-colnames(set[[dep]])
  for(e in exp){
    scores<-na.omit(score_filtered[,dep])
    plot(scores,m_filtered[names(scores),e],xlab=dep,ylab=e,
         col=rownames(m_filtered)%in%labels[[dep]]+1)
  }
}

#using DESeq results
evalGoodDiff<-lapply(names(finalModList),function(x){
  features<-rownames(resDESeqScreen[[x]])
  eval2Dn2[[x]][,features,drop=FALSE]
})
names(evalGoodDiff)<-names(finalModList)
#gives the set of differentially expressed features, the extremely strong ones
#are determined by the template matching

#find genes that are essentially silenced
ko<-list()
for(dep in names(evalGoodDiff)){
  keep<-sapply(colnames(evalGoodDiff[[dep]]),function(x){
    median(m_filtered[labels[[dep]],x])<=0.1
  })
  ko[[dep]]<-evalGoodDiff[[dep]][,keep,drop=FALSE]
}

##---------------------------Relaxed Lasso Models-----------------------------------
###using glmnet's relaxed lasso
relaxo<-list()
for(dep in names(finalModList)){
  x<-m_filtered[labels[[dep]],colnames(eval2Dn2[[dep]]),drop=FALSE]
  x<-predict(preProcess(x,method=c('center','scale')),x)
  y<-score_filtered[labels[[dep]],dep,drop=FALSE]
  if(ncol(x)>1){
    mod<-cv.glmnet(x,y,relax=TRUE,nfolds=nrow(x))
    relFit<-mod$relaxed
    optimal<-which(mod$lambda==mod$relaxed$lambda.min) #get the optimal lambda
    relaxo.r2<-mod$glmnet.fit$relaxed$dev.ratio[optimal]
    relaxo.rmse<-sqrt(mod$cvm[optimal])
    relaxo[[dep]]<-list(mod,relaxo.r2,relaxo.rmse,coef(mod))
  }
  else
    relaxo[[dep]]<-lm(y~x)
}

##---------------------------------Extract See-saw-----------------------------------
bigoCorMat<-lapply(c('FOXA1','SOX10','PAX8'),function(dep){
  features<-colnames(eval2Dn2[[dep]])
  cl<-labels[[dep]]
  corMat<-matrix(rep(1,length(features)^2),nrow=length(features),
                 dimnames = list(features,features))
  if(length(features)>1){
    tries<-combn(features,2)
    for(i in 1:ncol(tries)){
      a<-tries[1,i];b<-tries[2,i]
      cor.ab<-cor(m_filtered[cl,a],m_filtered[cl,b])
      corMat[a,b]<-cor.ab;corMat[b,a]<-cor.ab
    }
  }
  return(corMat)
})
names(bigoCorMat)<-c('FOXA1','SOX10','PAX8')

seesaw<-list()

for(dep in c("FOXA1","SOX10","PAX8")){
  feature_index<-relaxo[[dep]][[4]]@i+1
  features<-rownames(relaxo[[dep]][[4]])[feature_index[-1]]
  mat<-bigoCorMat[[dep]][features,]
  mat<-mat[,apply(mat,2,function(y){sum(abs(y)>=0.75)>0})]
  features<-union(features,colnames(mat))
  features<-features[sapply(features,function(x){
    abs(cor(score_filtered[labels[[dep]],dep],m_filtered[labels[[dep]],x]))>0.4
  })]
  seesaw[[dep]]<-eval2Dn2[[dep]][,features]
}

save(relaxo,file='relaxedLasso.RData')
save(seesaw,file='seesaw.RData')
