require(Biostrings)
require(parallel)
require(VirFinder)

#### parVF.run ####
# parVF.run checks sequence length to predict with the appropiate model
parVF.run <- function (seqFaIn){
  # vector of individual nucleotides
  seqFa <- strsplit(x=as.character(seqFaIn),split="",fixed=T)[[1]]
  predResult <- NULL
  
  # kmer occurence count (k = 8)
  featureOut <- countSeqFeatureCpp(seqFa, modEPV) 
  featureOut_kmerCount <- featureOut$kmerCount # same result using VF.trainMod8mer
  # sequence length
  seqLength <- length(seqFa)
  
  # select previously trained logistic regression model based on sequence length. 
  # 3 models for lengths < 1kb, 1kb - 3kb, > 3kb are part of modEPV or default VF.trainMod8mer
  if (seqLength < 1 * 1000) {
    lasso.mod <- attr(modEPV, "lasso.mod_0.5k")
    nullDis <- attr(modEPV, "nullDis_0.5k")
  }
  else if (seqLength < 3 * 1000) {
    lasso.mod <- attr(modEPV, "lasso.mod_1k")
    nullDis <- attr(modEPV, "nullDis_1k")
  }
  else {
    lasso.mod <- attr(modEPV, "lasso.mod_3k")
    nullDis <- attr(modEPV, "nullDis_3k")
  }
  # predict with chosen lasso model using kmer counts as variables to maximize AUC
  lasso.pred <- predict(lasso.mod, t(as.matrix(featureOut_kmerCount)),
                        type = "response")
  pvalue <- mean(nullDis > as.numeric(lasso.pred))
  print(paste("len", seqLength,
              "score", round(lasso.pred, 4), "pvalue", round(pvalue,
                                                             4)))
  predResult <- rbind(predResult, c(seqLength, lasso.pred, pvalue))
  colnames(predResult) <- c("length", "score", "pvalue")
  predResult_df <- as.data.frame(predResult)
  predResult_df$length <- as.numeric(as.character(predResult_df$length))
  predResult_df$score <- as.numeric(as.character(predResult_df$score))
  predResult_df$pvalue <- as.numeric(as.character(predResult_df$pvalue))
  return(predResult_df)
}

#### parVF.pred ####
# parVF.pred calls parVF.run for each sequence and returns 
# the data.frame for the whole assembly fasta file
parVF.pred <- function(inFaFile, cores=2) {
  predResult_df <- NULL
  # readDNAStringSet extracts individual sequences from FASTA file
  dnastringset <- readDNAStringSet(inFaFile)
  # parVF.run for each sequence, parallelized
  predResult_df <- do.call(rbind, mclapply(dnastringset, parVF.run, mc.preschedule=F, mc.cores=cores))
  predResult_df$name <- rownames(predResult_df)
  try.res <- try(predResult_df$qvalue <- VF.qvalue(predResult_df$pvalue), silent = T)
  
  # qvalue can't be calculated for smaller fasta files where uniform distribution 
  # of p-values is not explored enough, i.e. values near boundaries (0,1) are missing
  if (class(try.res) != "try-error") {
    predResult_df <- predResult_df[,c("name","length","score","pvalue", "qvalue")]
  }
  else {
    predResult_df <- predResult_df[,c("name","length","score","pvalue")]
  }
  predResult_df <- predResult_df[names(dnastringset),]
  return(predResult_df)
}

#### findFaFiles ####
# finds Trinity FASTA assembly files under path
# example: PHYLUM superfolder as path
findFaFiles <- function(path){
  setwd(path)
  current_dir <- paste(c(getwd(), '/'), collapse = "")
  assembly <- dir(pattern = glob2rx("*.Trinity.fasta"), recursive = T)
  abs.paths <- c()
  #ignore *Trinity.fasta files inside read_partitions/ directory
  j <- 0
  for (i in 1:length(assembly)){
    if ("read_partitions" %in% strsplit(assembly[i], "/")[[1]] == FALSE){
      j <- j +  1
      abs.paths[j] <- paste(c(current_dir, assembly[i]), collapse = "")
    }
  }
  return(abs.paths)
}



#### parVF.batch ####
# parVF.batch takes a vector of absolute paths to FASTA files and calls parVF.pred for each of them, 
# creates directories and writes data.frame in csv with VirFinder results
parVF.batch <- function(FaFilePath, cores){
  # getting directory where assembly is found
  path.end <- strsplit(FaFilePath, "/")[[1]][length(strsplit(FaFilePath, "/")[[1]])]
  subfolder <- paste(c(paste(strsplit(FaFilePath, "/")[[1]][1:length(strsplit(FaFilePath, "/")[[1]])-1], collapse = "/"), "/"), collapse = "")
  # accession
  name <- strsplit(path.end, ".Trinity.fasta")[[1]]

  preds <- parVF.pred(FaFilePath, cores = cores)
  
  #to output file
  suppressWarnings(dir.create(paste(c(subfolder, "VirFinder"), collapse = "")))
  write.table(preds, file = paste(c(subfolder, "VirFinder/", name, "_VirFinder.csv"), collapse = ""), row.names = FALSE, col.names = TRUE, sep = ",")
}





#### Code to run ####
# load model trained with prokaryotic and eukaryotic viruses
load("/home/paualico/.conda/envs/R_env/VirFinder/VF.modEPV_k8.rda")

# provide path above assembly files to findFaFiles to generate vector of paths and provide number of cores to run
lapply(findFaFiles("/home/paualico/samples/pau/PHYLUM/"), parVF.batch, cores = 25)


# important modifications from https://github.com/rec3141/VirFinder/blob/master/linux/VirFinder/R/parVF.pred.R

# deleted:
# rmWordID <- attr(VF.trainMod8mer, "rmWordID_0.5k") #modEPV lacks rmWordID and is not needed for correct dimensionality

# modified:
# model load was put upstream instead of being loaded for every sequence as in original script
# generalization for our file structure


# added:
# qvalue column added only if enough data to calculate
