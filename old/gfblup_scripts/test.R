library(qgg)
# source("https://bioconductor.org/biocLite.R")
# biocLite("snpStats")
library(snpStats)

phenotypes <- read.csv("./data/input_files/AfterBlupsAbsLevels04222014.csv")
phenotypes$Taxa <- paste0("A", phenotypes$Taxa) # append 'A' to taxa names
phenotypes <- phenotypes[-78,] # remove individual with no genotypic data
ratios <- read.csv("~/Documents/AA-GenomicPrediction/data/input_files/Blup-transformedforGWAS.csv", header=TRUE)
ratios$Taxa <- paste0("A", ratios$Taxa)
phenotypes <- merge(phenotypes, ratios, by = "Taxa")
phenotypes <- replace(phenotypes, is.na(phenotypes), NA)

# Construct the genotype matrix (W)
W_full <- read.plink(bed="./data/GS_subset/plinkGeneOmeSubset.bed", bim="./data/GS_subset/plinkGeneOmeSubset.bim", fam="./data/GS_subset/plinkGeneOmeSubset.fam", na.strings=c("-9"), sep=".")

W <- W_full$genotypes
rm(W_full)
rownames(W) <- paste("A", rownames(W), sep="")
W <- as.matrix(as.data.frame(W))
class(W) <- "numeric"
W <- scale(W, center=TRUE, scale=TRUE)
dim(W)
W[1:5,1:5]

W <- W[rownames(W)%in%phenotypes$Taxa,] # remove individuals with missing phenotypic data
dim(W)

# extract names of SNPs with missing data and remove from W
m_missing <- colnames(W)[colSums(is.na(W)) > 0] # markers with missing data in full set
W <- W[,!colSums(is.na(W)) > 0]
dim(W)

write.matrix(W, "./data/input_files/W_mat.txt")

## EDIT - this no longer works
## this appears to be different as of 14 August 2018
Glist <- prepG(bedfiles="./data/GS_subset/plinkGeneOmeSubset.bed",
               bimfiles="./data/GS_subset/plinkGeneOmeSubset.bim",
               famfiles="./data/GS_subset/plinkGeneOmeSubset.fam", study="AA-gfblup",
               fnRAW="./data/GS_subset/plinkGeneOmeSubset.out")

Glist$ids <- paste("A", Glist$ids, sep="")

W2 <- getW(Glist)