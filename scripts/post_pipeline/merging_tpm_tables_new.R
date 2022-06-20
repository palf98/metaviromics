### EN TERCER LUGAR

args = commandArgs(trailingOnly=TRUE)
# args[1]: path donde buscar las tablas individuales por muestra, puede estar por encima
# args[2]: donde se quiera escribir las tablas

tables_per_method <- function(final_substring, path){
  lista = list()
  setwd(path)
  for (i in 1:length(dir(pattern = glob2rx(final_substring), recursive = T))){
    lista[[i]] <- read.table(dir(pattern = glob2rx(final_substring), recursive = T)[[i]], header = T, sep =
"\t")
  }
  df <- plyr::rbind.fill(lista)
  df[is.na(df)] <- 0
  return(df)
}

diam <- tables_per_method("*accession_ID.tsv", args[1])

setwd(args[2])
write.csv(diam, "merged_diamond_tpm.csv", row.names = F)
