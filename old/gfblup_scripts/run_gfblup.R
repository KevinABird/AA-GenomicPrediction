###############################################################################################################
##  Compare predictive abilities for GBLUP vs GFBLUP w/amino acid SNP subset
##  S. Turner
##  Jan. 31 2018
###############################################################################################################

# load libraries
library(qgg)
library(plyr)

###############################################################################################################
## Load and prepare data

load("/home/sarahturner/Documents/Bird_et_al/input_files/gfblup_input.Rdata")
# use command line arguments
traits <- colnames(phenotypes[,2:28])
print(traits)

###############################################################################################################
### Run GBLUP and GFBLUP 

fitGB <- list()
fitGF <- list()

for(name in traits) {
  # remove individuals with missing phenotypes
  p_missing <- phenotypes$Taxa[is.na(phenotypes[,name])]
  pheno_tmp <- phenotypes[!phenotypes$Taxa %in% p_missing,]
  G_tmp <- G[!rownames(G) %in% p_missing, !colnames(G) %in% p_missing]
  W_tmp <- W[!rownames(W) %in% p_missing,]
  L <- ncol(G_tmp)
  
  # create lists of subsets for different models
  # setsGB <- list(F=colnames(W)) # define SNP set for gblup model - not really necessary?
  
  # formula
  fm <- as.formula(paste(name, "~", 1, sep=""))
  print(fm)
  
  # parameters for REML analyses and cross validation
  n <- ncol(G_tmp)
  fold <- 10
  nsets <- 50
  
  # ids for validation sets
  validate <- replicate(nsets, sample(1:n, as.integer(n / fold))) # matrix input
  
  # gblup model
  fitGB[[name]] <- gfm(fm=fm, W=W_tmp, G=list(G_tmp), data=pheno_tmp, validate=validate)
  
  # define marker sets
  setsGF <- list(AAS=m_aas) # define SNP set for gfblup model using AA snps
  # fit the model! 
  fitGF[[name]] <- gfm(fm=fm, W=W_tmp, G=list(G_tmp), sets=setsGF, data=pheno_tmp, validate=validate)
}

# save predictive abilities and sigmas
gblup_pa <- data.frame(method="gblup")
gblup_sigma <- list()

for (name in traits){
  df_pa <- as.data.frame(fitGB[[name]]$pa)
  colnames(df_pa) <- name
  # df <- colnames(df, c(name))
  gblup_pa <- cbind(gblup_pa, df_pa)
  df_sigma <- as.data.frame(fitGB[[name]]$sigmas)
  gblup_sigma[[name]] <- df_sigma
}

gfblup_pa <- data.frame(method="gfblup")
gfblup_sigma <- list()

for (name in traits){
  df_pa <- as.data.frame(fitGF[[name]]$pa)
  colnames(df_pa) <- name
  # df <- colnames(df, c(name))
  gfblup_pa <- cbind(gfblup_pa, df_pa)
  df_sigma <- as.data.frame(fitGF[[name]]$sigmas)
  gfblup_sigma[[name]] <- df_sigma
}

# save output
save(gblup_pa, gblup_sigma, gfblup_pa, gfblup_sigma, file="/home/sarahturner/Documents/Bird_et_al/GBLUP_vs_GFBLUP.RData")


### Run for BCAT2 gene

fitGB <- list()
fitGF <- list()

for(name in traits) {
  # remove individuals with missing phenotypes
  p_missing <- phenotypes$Taxa[is.na(phenotypes[,name])]
  pheno_tmp <- phenotypes[!phenotypes$Taxa %in% p_missing,]
  G_tmp <- G[!rownames(G) %in% p_missing, !colnames(G) %in% p_missing]
  W_tmp <- W[!rownames(W) %in% p_missing,]
  L <- ncol(G_tmp)
  
  # create lists of subsets for different models
  # setsGB <- list(F=colnames(W)) # define SNP set for gblup model - not really necessary?
  
  # formula
  fm <- as.formula(paste(name, "~", 1, sep=""))
  print(fm)
  
  # parameters for REML analyses and cross validation
  n <- ncol(G_tmp)
  fold <- 10
  nsets <- 50
  
  # ids for validation sets
  validate <- replicate(nsets, sample(1:n, as.integer(n / fold))) # matrix input
  
  # gblup model
  fitGB[[name]] <- gfm(fm=fm, W=W_tmp, G=list(G_tmp), data=pheno_tmp, validate=validate)
  
  # define marker sets
  setsGF <- list(BCAT2 = c("S5372", "S5373", "S5374", "S5375", "S5376", "S5377", "S5378", "S5379", "S5380", "S5381", 
                            "S5382", "S5383", "S5384", "S5385", "S5386"))
  # fit the model! 
  fitGF[[name]] <- gfm(fm=fm, W=W_tmp, G=list(G_tmp), sets=setsGF, data=pheno_tmp, validate=validate)
}

for (name in traits){
  df_pa <- as.data.frame(fitGF[[name]]$pa)
  colnames(df_pa) <- name
  # df <- colnames(df, c(name))
  gfblup_pa <- cbind(gfblup_pa, df_pa)
  df_sigma <- as.data.frame(fitGF[[name]]$sigmas)
  gfblup_sigma[[name]] <- df_sigma
}


for (name in traits){
  print(mean(gfblup_pa[[name]]) > mean(gblup_pa[[name]]))
}
