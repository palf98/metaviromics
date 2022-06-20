# 1. finds VirFinder csv output files (VFcsv) under path provided in first argument
# 2. processes VFcsv into new csv files with new columns:
#    a) fdr estimated by BH method with p.adjust()
#    b) logical vector with TRUE meaning fdr < 0.05
#    c) logical vector with TRUE meaning that the row has one of the 30% lowest p-values

args = commandArgs(trailingOnly=TRUE)

findVFcsvFiles <- function(path){
  setwd(path)
  current_dir <- paste(c(getwd(), '/'), collapse = "")
  csv <- dir(pattern = glob2rx("*VirFinder.csv"), recursive = T)
  abs.paths <- c()
  j <- 0
  for (i in 1:length(csv)){
    j <- j +  1
    abs.paths[j] <- paste(c(current_dir, csv[i]), collapse = "")
  }
  return(abs.paths)
}

virfinderpost <- function(vfcsv){
  name <- paste(strsplit(x = vfcsv, split = ".", fixed = T)[[1]][1], "postprocess.csv", sep = "_")
  #open data.frame from csv file
  df <- na.omit(read.csv(vfcsv, header = T))
  fdrs <- p.adjust(p = df$pvalue, method = "fdr")
  percentil30 <- df$pvalue < quantile(df$pvalue, probs = 0.3)
  df <- cbind(df, data.frame("fdr" = fdrs, "fdr.less0.05" = fdrs < 0.05, "percentil30" = percentil30))
  write.table(df, file = name, row.names = FALSE, col.names = TRUE, sep = ",")
  print(paste(name, " WRITTEN!", sep = ""))
}

lapply(findVFcsvFiles(args[1]), virfinderpost)
