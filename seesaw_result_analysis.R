library(CARS)
files<-list.files()
files<-files[grep('mainframe',files)]

dependencies<-vector(mode='list',length=ncol(score_filtered))
names(dependencies)<-colnames(score_filtered)

for(file in files){
load(file)
c<-vector()
for(i in 1:length(all_gene_imp)){
  if(dim(all_gene_imp[[i]])[1]>0)
    c<-c(c,i)
}
##finds index in the result of dependencies that has a found a potential see-saw pair

dat<-c(c,names(all_gene_imp)[c]) 
##gets the names of the aforementioned dependencies
dependencies[dat]<-all_gene_imp[dat]
}


dependencies<-dependencies[sapply(dependencies,function(x){length(x)>0})]
dependencies<-lapply(1:length(dependencies),function(i){
  data.frame(dependency=rep(names(dependencies[i]),nrow(dependencies[[i]])),
                       seesaw.gene=rownames(dependencies[[i]]),dependencies[[i]])
})
dependencies<-do.call('rbind',dependencies)


##----------------------------Re-verify linearity------------------------------------
#gather the r2 for each model I decided to fit
seesaw.r2<-sapply(1:nrow(dependencies),function(i){
  pair<-dependencies[i,]
  scores<-score_filtered[,pair[,1]]
  kept<-which(scores<=pair[,7])
  (cor(score_filtered[kept,pair[,1]],m_filtered[kept,pair[,2]]))^2
})
dependencies$r2<-seesaw.r2
#most R2 are weak; good ones are >0.3, >0.2 is also fine

##--------------------------Automatic Plot Formatting-------------------------------
row<-dependencies['',]
scores<-score_filtered[,"SMARCA2"]
expression<-m_filtered[,"SMARCA4"]
#true.dep<-which(scores<=row[[7]])
true.dep<-which(scores<=(-0.8))
line<-lm(expression[true.dep]~scores[true.dep])
plot(scores,expression,xlab="SOX10",ylab="SDC3")
plot(scores[true.dep],expression[true.dep],xlab="SOX10",ylab="SDC3")
abline(v=-2)
abline(line)


##--------------------------Hunt for Good Features----------------------------------
#which genes have many see-saw candidates?
depCount<-data.frame(dependency=unique(dependencies[,1]),count=rep(0,length(unique(dependencies[,1]))))
rownames(depCount)<-depCount[,1]
for(dep in unique(dependencies[,1])){
  depCount[dep,2]<-length(unique(dependencies[dependencies[,1]==dep,2]))
}

depCand1<-rownames(depCount)[depCount$count>1]
#57 dependencies have at least 2 dependency candidates
stop()

##-------------------------Filter By See-Saw Scores--------------------------------
generate_quant_sep=function(dep,exp,thresh){
  quantPos=quantile(exp[dep<=thresh])
  quantNeg=quantile(exp[dep>thresh])
  seeScore=quantPos[2]-quantNeg[4]
  sawScore=-(quantPos[4]-quantNeg[2])
  
  list(seeScore,sawScore)
}
new<-dependencies[,1:2]
new$see<-NA
new$saw<-NA
for(i in 1:nrow(dependencies)){
  cl<-rownames(score_filtered)[which(!is.na(score_filtered[,dependencies$dependency[i]]))]
  dep<-score_filtered[cl,dependencies$dependency[i]];exp<-m_filtered[cl,dependencies$seesaw.gene[i]]
  a<-generate_quant_sep(dep,exp,-2)
  new$see[i]<-a[[1]]
  new$saw[i]<-a[[2]]
}
ss_filt<-dependencies[new$see>=0.5|new$saw>=0.5,] 

##-------------------------R2 Filter------------------------------------------------
ss_filt2<-ss_filt[ss_filt$neg2r2>=0.25,]

##--------------------------------------Linear models with split at dep score=-2-----------------------
#becuase brier scoring churned out BS
neg2seesaw.r2<-sapply(1:nrow(dependencies),function(i){
  pair<-dependencies[i,]
  scores<-score_filtered[,pair[,1]]
  kept<-which(scores<=-2)
  (cor(score_filtered[kept,pair[,1]],m_filtered[kept,pair[,2]]))^2
})
dependencies$neg2r2<-neg2seesaw.r2

#now select those with -2 R2>0.2

depneg2R2<-dependencies[dependencies$neg2r2>=0.2,] #186 pairs left
unique(depneg2R2$dependency) 

