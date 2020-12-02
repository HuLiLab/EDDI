library(dplyr)
library(DESeq2)
library(CePa)
library(biomaRt)
library(EnsDb.Hsapiens.v79)
#make sure Ive run init dataload
##---------------------------Generate Counts Matrix-----------------
counts<-read.gct('CCLE_DepMap_18q3_RNAseq_reads_20180718.gct')
#add a pseudo count 1 to make sure size estimation won't break
counts<-counts+1

common<-sapply(rownames(m_filtered),function(x){
  grep(x,colnames(counts))
})
common[["EN_ENDOMETRIUM"]]<-170
common<-unlist(common)

counts<-counts[,common]
colnames(counts)<-rownames(m_filtered)

mart<-useMart('ENSEMBL_MART_ENSEMBL')
mart<-useDataset('hsapiens_gene_ensembl',mart)

names<-gsub('\\..*',"",rownames(counts))
annot<-getBM(mart=mart,filter='ensembl_gene_id',values=names,
             attributes=c("ensembl_gene_id",
                          "gene_biotype", "external_gene_name"))
annotMatch<-data.frame(rownames(counts)[match(annot$ensembl_gene_id,names)],annot)
#THIS STILL ISN'T THE RIGHT AMOUNT OF THINGS AJADJFSDKLAFF;needs fixing
#Just make sure my diff dep genes are there
geneID<-sapply(rownames(counts),function(x){
  annotMatch[annotMatch[,1]==x,4]
})
geneID<-unlist(geneID)

#find which m_filtered genes went unmatched
unmatched<-setdiff(colnames(m_filtered),geneID)

#try to recover some unmatched
recover <- ensembldb::select(EnsDb.Hsapiens.v79, keys= unmatched, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
match<-sapply(recover$GENEID,grep,rownames(counts)) #find count matrix ids that fit
match<-match[sapply(match,function(x){length(x)>0})] #remove the ones without a fit
match<-unlist(match)
#get the count ID for each recovered symbol and add to geneID
geneIDadd<-recover[recover$GENEID%in%names(match),1]
names(geneIDadd)<-rownames(counts)[match]
geneID<-c(geneID,geneIDadd)

unmatched2<-setdiff(colnames(m_filtered),geneID) #71 left to go



#format it like the names that were matched
countsFilt<-counts[names(geneID),] #select the counts
rownames(countsFilt)<-geneID #replace the ensemblid with the gene ID
countsFilt<-countsFilt[intersect(colnames(m_filtered),geneID),] #get those that are in
#m_filtered



##------------------------Generate Sample Tables----------------------
#for each dependency, make a table stating whether a cell is dep vs nd
#and specify the genes I wish to look at
sampleTables<-lapply(names(finalModList),function(i){
  depStatus<-rep('NonDependent',440)
  depStatus[match(labels[[i]],colnames(countsFilt))]<-'dependent'
  annot<-data.frame(colnames(countsFilt),depStatus)
  rownames(annot)<-colnames(countsFilt)
  colnames(annot)<-c('sample','condition')
  annot
})
names(sampleTables)<-names(finalModList)

##-------------------------Construct DESeq Dataset|Ignore--------------------
dsData<-lapply(sampleTables,function(x){
  dds<-DESeqDataSetFromMatrix(countData=countsFilt[x[[1]],,drop=FALSE],colData=x[[2]],
                              design=~condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
})


##-------------------------Run Differential Expression Analysis|Ignore------
deseqDat<-list()
for(dep in names(evalGoodDiff)){
  tryCatch({
  deseqDat[[dep]]<-DESeq(dsData[[dep]],fitType = 'mean')
  },error=function(e){cat("Error:",conditionMessage(e),"\n")})
}
#apparently, some of the over-under expressed genes are within 2 mag
manual<-setdiff(names(evalGoodDiff),names(deseqDat))
for(dep in manual){
  dds<-dsData[[dep]]
  dds<-estimateSizeFactors(dds)
dds<-estimateDispersionsGeneEst(dds)
dispersions(dds)<-mcols(dds)$dispGeneEst
deseqDat[[dep]]<-nbinomWaldTest(dds)
}

resDESeq<-lapply(deseqDat,results,
                 contrast = c('condition','dependent','NonDependent'))

##---------------------DESeq except I use the whole counts matrix each time----------

dsDataFull<-lapply(sampleTables,function(x){
  dds<-DESeqDataSetFromMatrix(countData=countsFilt,colData=x,
                              design=~condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
})
resDESeqFull<-list()
for(dep in names(finalModList)){
  tryCatch({
    deseq<-DESeq(dsDataFull[[dep]],fitType = 'mean')
    resDESeqFull[[dep]]<-results(deseq,contrast = c('condition','dependent','NonDependent'))
  },error=function(e){cat("Error:",conditionMessage(e),"\n")})
}

#resDESeqFull<-lapply(deseqDatFull,results,
                 #contrast = c('condition','dependent','NonDependent'))

save(resDESeqFull,file='resDESeq.RData')

resDESeq2fold<-lapply(names(resDESeqFull),function(x){
  dat<-resDESeqFull[[x]]
  diff<-which(abs(dat$log2FoldChange)>=1.0 & dat$padj<=0.05)
  dat<-dat[diff,]
  dat<-dat[order(dat$padj,decreasing=FALSE),]
})
names(resDESeq2fold)<-names(resDESeqFull)
#each of these has waay more than I need, so reduce down to what was 
#included in their dep set, but make note of any potential predictive features

resDESeqScreen<-lapply(names(resDESeqFull),function(x){
  dat<-resDESeq2fold[[x]]
  keep<-intersect(rownames(dat),colnames(eval2Dn2[[x]]))
  dat<-as.data.frame(dat[keep,])
})
names(resDESeqScreen)<-names(resDESeqFull)
#okay so this says about everything is differentially expressed
#there are only some instances where the gene list is differentf

View(resDESeqScreen[['FOXA1']][colnames(evalGoodDiff[['FOXA1']]),])
#the log2 fold changes are on the higher side here

