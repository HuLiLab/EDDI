source("init_dataLoad.R")
## set threshold (-2 as detailed in CELL paper) and apply threshold
thresh=-2
useThreshold = function(x){
  sum(na.omit(x)<=thresh)
}
withinThresh=apply(score_filtered, 2, useThreshold)
topWithinThresh=withinThresh[order(-withinThresh)][1:2000]  #2000 standard

