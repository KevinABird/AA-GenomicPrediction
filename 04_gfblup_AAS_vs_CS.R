###############################################################################################################
##  GFBLUP for Arabidopsis amino acid panel
##  This script uses GFBLUP to compare the amino acid SNP (AAS) subset to the control SNP (CS) subset
##  S. Turner
##  Jan. 31 2018; Updated Aug. 16 2018
###############################################################################################################

# load libraries
library(qgg)

###############################################################################################################
## Load and prepare data

load("./data/input_files/gfblup_input.Rdata")
# use command line arguments
args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
traits <- colnames(phenotypes[,2:54])
pheno <- traits[index]
print(pheno)

# load principal components
pcs <- read.table("./data/input_files/principal_components.txt", header=TRUE)

###############################################################################################################
### Run GBLUP (infintesimal model)
# remove individuals with missing phenotypes
p_missing <- phenotypes$Taxa[is.na(phenotypes[,pheno])]
phenotypes <- phenotypes[!phenotypes$Taxa %in% p_missing,]
# G <- G[!rownames(G) %in% p_missing, !colnames(G) %in% p_missing]
W <- W[!rownames(W) %in% p_missing,]
# L <- ncol(G)
pcs <- pcs[!rownames(pcs) %in% p_missing,]

# create lists of subsets for different models
# setsGB <- list(F=colnames(W)) # define SNP set for gblup model - not really necessary?

# formula
# fm <- as.formula(paste(pheno[1], "~", "PC1", "+", "PC2", sep=""))
fm <- ifelse(pheno[1]=="val", paste0(pheno[1], "~ PC1 + PC2"), paste0(pheno[1], "~ 1"))
fm <- as.formula(fm)
print(fm)

# parameters for REML analyses and cross validation
n <- nrow(W)
fold <- 10
nsets <- 1000

# ids for validation sets
validate <- replicate(nsets, sample(1:n, as.integer(n / fold))) # matrix input

input_data <- cbind(phenotypes, pcs)

# standard gblup model
# dir.create(paste0("./gfblup_results/", pheno, "/"))

fitGB <- gfm(fm=fm, W=W, data=input_data, validate=validate, mkplots=FALSE)
# save(fitGB, file=paste0("./gfblup_results/", pheno, "/gblup_", pheno, ".RData"))
fitGB$sigmas
df_gs <- as.data.frame(fitGB$sigmas)
write.table(df_gs, paste0("./gfblup_results/gs_sigmas_", pheno, ".txt"))
gs_llik <- fitGB$fit$llik
write.table(gs_llik, paste0("./gfblup_results/gs_llik_", pheno, ".txt"))
gs_pa <- as.data.frame(fitGB$pa)
write.table(gs_pa, paste0("./gfblup_results/gs_pa_", pheno, ".txt"))

## GFBLUP - AA set
setsAA <- list(AA=m_aas)
fitAA <- gfm(fm=fm, W=W, sets=setsAA, data=input_data, validate=validate, mkplots=FALSE)
# save(fitAA, file=paste0("./gfblup_results/", pheno, "/fitAA_", pheno, ".RData"))
fitAA$sigmas
df_aa <- as.data.frame(fitAA$sigmas)
write.table(df_aa, paste0("./gfblup_results/aa_sigmas_", pheno, ".txt"))
aa_llik <- fitAA$fit$llik
write.table(aa_llik, paste0("./gfblup_results/aa_llik_", pheno, ".txt"))
aa_pa <- as.data.frame(fitAA$pa)
write.table(aa_pa, paste0("./gfblup_results/aa_pa_", pheno, ".txt"))

## GFBLUP - control sets
# gfblup model comparing AAS and CS subsets
# exports input for binomial test
cs_results <- NULL
cs_llik <- NULL
cs_pa <- NULL

for (i in 1:1000){
  print(i)
  # read in control marker data
  m_control <- read.table(paste0("./data/control_markers/plink_control_snplist_", i, ".txt"))
  m_control <- as.character(m_control[,1])
  # define marker sets
  setsC <- list(C=m_control) # define SNP set for gfblup model using random snps
  # fit the model!
  fitC <- gfm(fm=fm, W=W, sets=setsC, data=input_data, validate=validate, mkplots=FALSE)
  # save(fitC, file=paste("./gfblup_results/", pheno, "/fitC_", pheno, i, ".RData", sep=""))
  # export sigmas and log likelihoods
  df_c <- as.data.frame(fitC$sigmas)
  llik <- fitC$fit$llik
  pa <- mean(fitC$pa)
  cs_results <- rbind(cs_results, df_c)
  cs_llik <- rbind(cs_llik, llik)
  cs_pa <- rbind(cs_pa, pa)
  write.table(cs_results, paste0("./gfblup_results/cs_sigmas_", pheno, ".txt"))
  write.table(cs_llik, paste0("./gfblup_results/cs_llik_", pheno, ".txt"))
  write.table(cs_pa, paste0("./gfblup_results/cs_pa_", pheno, ".txt"))
  cat("Done!\n")
}
