abun <- read.csv("~/Desktop/CoDA-virus-family-abun/absabun_merged_001.csv", header = T)
df <- read.delim("data/db_def_con_DNA.tsv", header = T)

subset <- dplyr::filter(df, run %in% abun$accession)

mezcolanza <- cbind(subset[order(subset$run),], abun[order(abun$accession),])
mezcolanza$bases <- mezcolanza$bases / 1e9 #Gbases

diff_fams <- function(comp.vec){
  length(comp.vec[which(comp.vec > 0)])
}
shannon <- function(comp.vec){
  obj <- t(comp.vec / sum(comp.vec) * log(comp.vec / sum(comp.vec)))
  obj[is.nan(obj)] <- 0
  return(-1 * sum(obj))
}

shannons <- apply(mezcolanza[,25:101], 1, shannon)

#evennes = shannon/log(richness)

nfams <- apply(mezcolanza[,25:101], 1, diff_fams)

plot(mezcolanza$bases, nfams) # family richness
plot(mezcolanza$bases, shannons) # Shannon diversity index
abline(v = 1:10)
plot(mezcolanza$bases, shannons / log(nfams)) # Evenness

plot(mezcolanza$bases, mezcolanza$size)




dat <- read.delim("data/db_def_solo_RNA.tsv", header = T)

sample_example <- data.frame()

n <- 1
(seed <- runif(1)*100000000)
set.seed(seed)
for (i in 1:length(unique(dat$species))){
  nacc_not_asmbld <- length(dplyr::filter(dat, bases > 3e9, size < 20, assembly_suyo == F, species == unique(dat$species)[i])[,1])
  nacc <- length(dplyr::filter(dat, bases > 3e9, size < 20, species == unique(dat$species)[i])[,1])
  # si el numero de accessions para la especie descontando los assemblies hechos es mayor que n
  if ((nacc - (nacc - nacc_not_asmbld)) > n - 1){
    acc_smpl <- sample(1:nacc_not_asmbld, n)
    sample_example <- rbind(sample_example, dplyr::filter(dat, bases > 3e9, size < 20, species == unique(dat$species)[i])[sample(1:length(dplyr::filter(dat, bases > 3e9, size < 20, species == unique(dat$species)[i], assembly_suyo == F)[,1]), n),])
  }
  # si es positivo y no mayor que n
  else if ((nacc - (nacc - nacc_not_asmbld)) > 0){
    sample_example <- rbind(sample_example, dplyr::filter(dat, bases > 3e9, size < 20, species == unique(dat$species)[i], assembly_suyo == F))
  }
}

table(sample_example$phylum)
round(table(sample_example$phylum) / length(sample_example$phylum) , 4)
length(unique(sample_example$run)) # 1834 run acc
mean(sample_example$size)


write.table(sample_example, file = "data/seed21338385_todo.tsv", row.names = F, sep = "\t")
dat <- read.delim("data/seed21338385_todo.tsv", header = T)


raros <- dplyr::filter(dat, phylum %in% c("Brachiopoda", "Bryozoa", "Chaetognatha", "Ctenophora", "Cycliophora", "Dicyemida", "Entoprocta", "Gastrotricha",
                                            "Gnathostomulida", "Hemichordata", "Kinorhyncha", "Loricifera", "Mollusca", "Nematomorpha", "Nemertea", "Onychophora",
                                            "Orthonectida", "Phoronida", "Placozoa", "Priapulida", "Rotifera", "Tardigrada", "Urochordata", "Xenacoelomorpha"))
resto <- dplyr::filter(dat, !(run %in% raros$run))

unique(resto$phylum)
unique(raros$phylum)
write.table(raros, file = "data/raros.tsv", row.names = F, sep = "\t")

plot(dat$size)

set.seed(1)
resto_5xphyl <- rbind(dplyr::filter(resto, phylum == "Annelida")[sample(dim(dplyr::filter(resto, phylum == "Annelida"))[1], 5),],
                      dplyr::filter(resto, phylum == "Arthropoda")[sample(dim(dplyr::filter(resto, phylum == "Arthropoda"))[1], 5),],
                      dplyr::filter(resto, phylum == "Cnidaria")[sample(dim(dplyr::filter(resto, phylum == "Cnidaria"))[1], 5),],
                      dplyr::filter(resto, phylum == "Echinodermata")[sample(dim(dplyr::filter(resto, phylum == "Echinodermata"))[1], 5),],
                      dplyr::filter(resto, phylum == "Nematoda")[sample(dim(dplyr::filter(resto, phylum == "Nematoda"))[1], 5),],
                      dplyr::filter(resto, phylum == "Platyhelminthes")[sample(dim(dplyr::filter(resto, phylum == "Platyhelminthes"))[1], 5),],
                      dplyr::filter(resto, phylum == "Porifera")[sample(dim(dplyr::filter(resto, phylum == "Porifera"))[1], 5),])
write.table(resto_5xphyl, file = "data/avanzadilla_resto_filos.tsv", row.names = F, sep = "\t")



# que sea logaritmico tambien mola porque no aumentan los falsos positivos con el numero de lecturas





library(zCompositions)
X <- abun[,6:82][-which(apply(abun[,6:82], 1, sum) == 0),]
imputed <- cmultRepl(X, method = "CZM", output = "prop", z.warning = 0.999)


Y <- clr(X)
Y
aitch <- adist(X)
eucl <- dist(as.matrix(Y))
eucl

