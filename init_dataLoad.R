ts=proc.time()[3]
options(stringsAsFactors=F)
library(dplyr)
library(randomForest)
library(caret)
library(doParallel)
library(scales)
##-------------------------------Load Data---------------------------------------------------------
fnge='CCLE_DepMap_18q3_RNAseq_RPKM_20180718.hgnc_sym_log2_rpkm_plus_1.rdata' #gene expressions file (format is .rdata)
fnscore='ExpandedGeneZSolsCleaned.csv' ## dependency scores file (format is csv)

m=get(load(fnge)) ## load gene expressions
score=read.table(fnscore, sep=',' ,header=T, row.names=1, check.names=F) #read the score csv file it
m_sd<-get(load("CCLE_DepMap_18q3_RNAseq_RPKM_20180718.hgnc_sym_log2_rpkm_plus_1.mean_sd.rdata"))# load r data w/ the gene expression standard deviations
## LOADS THE CUT DATA THAT IS GENERATED FROM seesaw_tscore_gene_importance_filter
#load("seesaw_results_gene_imp_cut.rdata")


cl1=colnames(score)## names of cell lines in scores
cl2=colnames(m)## names of cell lines in gene expression data
##> cl2[1:10]
## [1] "22RV1_PROSTATE (ACH-000956)"                        
## [2] "2313287_STOMACH (ACH-000948)"                       
## [3] "253JBV_URINARY_TRACT (ACH-000026)"                  
## [4] "253J_URINARY_TRACT (ACH-000011)"                    
## [5] "42MGBA_CENTRAL_NERVOUS_SYSTEM (ACH-000323)"         
## [6] "5637_URINARY_TRACT (ACH-000905)"                    
## [7] "59M_OVARY (ACH-000520)"                             
## [8] "639V_URINARY_TRACT (ACH-000973)"                    
## [9] "647V_URINARY_TRACT (ACH-000896)"                    
##[10] "697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE (ACH-000070)"

cl2_=stringr::str_split(cl2, ' ', simplify=T)[,1] ## splits the gene expression cell line names from its ID in the original database
##> str(cl1)
## chr [1:501] "22RV1_PROSTATE" "697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" ...
##> str(cl2)
## chr [1:1156] "22RV1_PROSTATE (ACH-000956)" ...
##> str(cl2_)
## chr [1:1156] "22RV1_PROSTATE" "2313287_STOMACH" ...

common=intersect(cl1, cl2_) ## checks for cell lines that exists in both datasets
##> str(common)
## chr [1:487] "22RV1_PROSTATE" "697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" ...

colnames(m)=cl2_ ## remove the database ID from m

score=score[, common] ## filter cell lines to what's common
m=m[, common] 
m<-t(m) ## transpose matrix
score<-t(score) ## transpose matrix
##-------------------------------------------------Data Processing, SDDI Project------------------------------------------------------------------
genes<-rownames(m_sd) ## names of all the genes
m_sd<-as.data.frame(m_sd, row.names=TRUE) ## convert standard deviation data to a data frame
m_sd$gene<-genes ## create a column with gene names
m_sd<-m_sd%>%arrange(desc(sd)) ## sort the genes my standard deviation largest-> smallest; done with the dplyr package
gene_most_var<-m_sd$gene[1:200] ## select names of the 200 genes with highest SD
NO_HAEM<-!grepl("HAE",rownames(m)) ## logical vector, haematopoietic=true, not=false
## sets the filter
filt=NO_HAEM
## select cells according to filter
m_filtered<-m[filt,]
score_filtered<-score[filt,]
