#Version 2 - IBD filtering now turned to a function 
# Filtering for low call rate samples
#Format: R --vanilla --slave --args 1 $PLINK_FILE < $PATH/QC.R 
args<-commandArgs(trailingOnly=TRUE)
state <- args[1]
plinkoutput <- args[2]

#current program mode 0 - unassigned;
#1 - list call rate outliers; 
#2 - list sample relatedness (IBD >= 0.9) 
#3 - IBD >= 0.185

#Defining IBD filtering function
IBDMode <- function(pihat, state){
  writeLines(sprintf("\nSelected for MODE %s - Sample Relatedness (IBD>=%.3f)\n",state,pihat))
  filename <- plinkoutput
  check_ibd <- read.table(filename,header=T)
  ibd_sorted <- check_ibd[order(-check_ibd$PI_HAT), ]
  write.csv(ibd_sorted,file=paste("list_ibd_",pihat,"_sorted.csv",sep=""),row.names=F)
  writeLines(sprintf("Generated a list of sorted IBD values: list_ibd%s_sorted.csv",pihat))
  ibdfiltered = check_ibd[(check_ibd$PI_HAT >= as.numeric(pihat)), ]
  ibdfiltered_sorted = ibdfiltered[order(-ibdfiltered$PI_HAT), ]
  write.table(ibdfiltered_sorted, file = paste("list.ibd",pihat,".filtered.txt",sep=""), row.names = F)
  writeLines(sprintf("Generated a list of filtered IBD values: list.ibd%s.filtered.txt",pihat))
}
if (state == 1) {
  writeLines("\nSelected for MODE 1 - Call Rate Outliers\n
Input file required as .imiss file from PLINK\n")
  filename <- plinkoutput
  miss <- read.table(filename,header=T)
  miss$call0.98 <- (1-miss$F_MISS)
  call.outliers = miss[miss$call0.98 < 0.98, ]
  call.outliers
  write.table(call.outliers, file = "list.calloutliers0.98.txt")
  print("Generated output as 'list.calloutliers0.98.txt'")
} else if (state == 2){
  IBDMode(0.9,state)
} else if (state == 3){
  IBDMode(0.185,state)
} else {
  print(sprintf("Invalid mode input. Retrieved %s when 1, 2 or 3 expected",state))
}


