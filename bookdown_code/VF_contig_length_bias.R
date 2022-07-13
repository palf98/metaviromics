
tables_per_method <- function(final_substring, path){
  lista = list()
  setwd(path)
  for (i in 1:length(dir(pattern = glob2rx(final_substring), recursive = T))){
    lista[[i]] <- read.table(dir(pattern = glob2rx(final_substring), recursive = T)[[i]], header = T, sep =
                               ",")[sample(dim(
                                 read.table(dir(pattern = glob2rx(final_substring), recursive = T)[[i]], header = T, sep =
                                              ","))[1], 5000),c(2,6)]
  }
  df <- plyr::rbind.fill(lista)
  df[is.na(df)] <- 0
  return(df)
}

VFpvals <- tables_per_method("*postprocess.csv", "/home/pau/TFM/projecte/")
write.csv(VFpvals, file = "data/sample50/VFpvals.csv", row.names = F)

par(mfrow = c(2,1), mar = c(1,3,1,1))
hist(VFpvals$fdr, breaks = 50)
plot(VFpvals$length, VFpvals$fdr, cex = 0.8, pch = 16)
abline(h = 0.05, lty = 2)
