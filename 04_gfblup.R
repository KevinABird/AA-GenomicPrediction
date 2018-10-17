###############################################################################################################
##  GFBLUP for Arabidopsis amino acid panel
##  This script runs the a genomic feature model for amino acid traits 
##  S. Turner
##  Jan. 31 2018; Updated Oct 16 2018
###############################################################################################################

###############################################################################################################
# TODO
# add principal components (see previous model selection results from GAPIT)
# generic input for marker classes (to test multiple pathways)
###############################################################################################################

# load libraries
library(qgg)

###############################################################################################################
## Load and prepare data
###############################################################################################################

load("data/input_files/gfblup_input.Rdata")
# use command line arguments
args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
traits <- colnames(phenotypes[,2:54])
pheno <- traits[index]
print(pheno)

###############################################################################################################
### TODO: EDIT THIS SECTION IN PREP DATA SCRIPT
###############################################################################################################

# load principal components
pcs <- read.table("data/input_files/principal_components.txt", header=TRUE)

### clean up missing data
p_missing <- phenotypes$Taxa[is.na(phenotypes[,pheno])]

phenotypes <- phenotypes[!phenotypes$Taxa %in% p_missing,]

pcs <- pcs[!rownames(pcs) %in% p_missing,]

W <- W[!rownames(W) %in% p_missing,]

input_data <- cbind(phenotypes, pcs)


###############################################################################################################
# set up cross validation
###############################################################################################################

n <- nrow(W)
fold <- 10
nsets <- 50

# ids for validation sets
validate <- replicate(nsets, sample(1:n, as.integer(n / fold))) # matrix input

###############################################################################################################
### Run GBLUP (infintesimal model)
###############################################################################################################

setsGB <- list(G = colnames(W)) # define SNP set for gblup model

# formula
fm <- ifelse(pheno[1]=="val", paste0(pheno[1], "~ PC1 + PC2"), paste0(pheno[1], "~ 1"))
fm <- as.formula(fm)
print(fm)

fitGB <- gfm(fm = fm, W = W, sets = setsGB, data = input_data, validate = validate, mkplots = FALSE)

# export variances
# df_gs <- as.data.frame(fitGB$sigmas)
write.table(fitGB$sigmas, paste0("gfblup_results2/gs_sigmas_", pheno, ".txt"))

# export log likelihood
gs_llik <- fitGB$fit$llik
write.table(gs_llik, paste0("gfblup_results2/gs_llik_", pheno, ".txt"))

# export predictive accuracy
gs_pa <- as.data.frame(fitGB$pa)
write.table(gs_pa, paste0("gfblup_results2/gs_pa_", pheno, ".txt"))

###############################################################################################################
## GFBLUP - AA set
###############################################################################################################

setsAA <- list(AA = m_aas, G = colnames(W)[!colnames(W) %in% m_aas])

fitAA <- gfm(fm = fm, W = W, sets = setsAA, data = input_data, validate = validate, mkplots = FALSE)

# export variances
# df_aa <- as.data.frame(fitAA$sigmas)
write.table(fitAA$sigmas, paste0("gfblup_results2/aa_sigmas_", pheno, ".txt"))

# export log likelihood
aa_llik <- fitAA$fit$llik
write.table(aa_llik, paste0("./gfblup_results2/aa_llik_", pheno, ".txt"))

# export predictive accuracy
aa_pa <- as.data.frame(fitAA$pa)
write.table(aa_pa, paste0("./gfblup_results2/aa_pa_", pheno, ".txt"))

###############################################################################################################
## GFBLUP - control sets
###############################################################################################################

# set up empty objects to store results
cs_results <- NULL
cs_llik <- NULL
cs_pa <- NULL

# loop through subsets (probably a more efficient way to do this)
for (i in 1:1000){
  print(i)
  
  # read in control marker data
  m <- m_control[,i]
  gs <- colnames(W)[!colnames(W) %in% m] # all other snps
  
  # define marker sets
  setsC <- list(C = m, G = gs) # define SNP set for gfblup model using random snps
  
  # fit the model!
  fitC <- gfm(fm = fm, W = W, sets = setsC, data = input_data, validate = validate, mkplots = FALSE)
  
  # export variances
  # df_c <- as.data.frame(fitC$sigmas)
  df_c <- fitC$sigmas
  cs_results <- rbind(cs_results, df_c)
  write.table(cs_results, paste0("gfblup_results2/cs_sigmas_", pheno, ".txt"))
  
  # export log likelihood
  llik <- fitC$fit$llik
  cs_llik <- rbind(cs_llik, llik)
  write.table(cs_llik, paste0("gfblup_results2/cs_llik_", pheno, ".txt"))
  
  # export predictive accuracy 
  pa <- mean(fitC$pa)
  cs_pa <- rbind(cs_pa, pa)
  write.table(cs_pa, paste0("gfblup_results2/cs_pa_", pheno, ".txt"))
  
  cat("Done!\n")
}

