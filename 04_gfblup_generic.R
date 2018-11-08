###############################################################################################################
##  GFBLUP for Arabidopsis amino acid panel
##  This script runs the a genomic feature model for amino acid traits 
##  S. Turner
##  Jan. 31 2018; Updated Oct 16 2018
###############################################################################################################

# load libraries
library(qgg)
library(regress)

###############################################################################################################
## Load and prepare data
###############################################################################################################

load("data/input_files/gfblup_input.Rdata")
# use command line arguments
args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
traits <- colnames(phenotypes[,2:55])
pheno <- traits[index]
print(pheno)
path <- args[2]
print(path)

dir <- paste0("gfblup_results/", unlist(strsplit(path, "/"))[2])

dir.create(paste0("gfblup_results/", unlist(strsplit(path, "/"))[2]))

# load principal components
pcs <- read.table("data/input_files/principal_components.txt", header=TRUE)

### clean up missing data
p_missing <- phenotypes$Taxa[is.na(phenotypes[,pheno])]

phenotypes <- phenotypes[!phenotypes$Taxa %in% p_missing,]

pcs <- pcs[!rownames(pcs) %in% p_missing,]

W <- W[!rownames(W) %in% p_missing,]

input_data <- cbind(phenotypes, pcs)

# create character vector of markers for specified pathway 
markers <- read.table(path)[,2]
markers <- as.character(markers)

###############################################################################################################
# set up cross validation
###############################################################################################################

n <- nrow(W)
fold <- 10
nsets <- 50

# ids for validation sets
validate <- replicate(nsets, sample(1:n, as.integer(n / fold))) # matrix input

###############################################################################################################
### Run GFBLUP (genomic prediction with feature set)
###############################################################################################################

# formula
fm <- if (pheno[1] %in% c("glu", "gly", "val", "BCAA", "gly_t", "val_t")) {
  paste0(pheno[1], "~ PC1 + PC2")
} else if (pheno[1] %in% "met"){
  paste0(pheno[1], "~ PC1")
} else {
  paste0(pheno[1], "~1")
}
  
fm <- as.formula(fm)
print(fm)

# specify marker sets
setsGF <- list(S = markers, G = colnames(W)[!colnames(W) %in% markers])

# fit GFBLUP model
fitGF <- gfm(fm = fm, W = W, sets = setsGF, data = input_data, validate = validate, mkplots = FALSE)
traceback()

# export variances
df_gf <- as.data.frame(fitGF$sigmas)
write.table(fitGF$sigmas, paste0(dir, "/gf_sigmas_", pheno, ".txt"), col.names = c("S", "G", "e"))

# export log likelihood
gf_llik <- fitGF$fit$llik
write.table(gf_llik, paste0(dir, "/gf_llik_", pheno, ".txt"))

# export predictive accuracy
gf_pa <- as.data.frame(fitGF$pa)
write.table(gf_pa, paste0(dir, "/gf_pa_", pheno, ".txt")) 
