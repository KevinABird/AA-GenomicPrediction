### Export sigmas from GFBLUP for plotting

load("../data/input_files/gfblup_input.Rdata")
# use command line arguments
args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
traits <- colnames(phenotypes[,2:28])
pheno <- traits[index]

result <- NULL

for (i in 1:1000){
    load(paste0("../gfblup_results/", pheno, "/fitGF_", pheno, i, ".RData"))
    print(paste(pheno, i, sep=" "))
    df <- as.data.frame(fitGF$sigmas)
    print(colMeans(df))
    result <- rbind(result, colMeans(df))	
    write.table(result, paste0("../gfblup_results/", pheno, "/sigmas_", pheno, ".txt"), 
    row.names=FALSE, col.names=TRUE)
}

df <- read.table(paste0("../gfblup_results/", pheno, "/sigmas_", pheno, ".txt"))
print(paste(pheno, "sigmas"))
print(df)
