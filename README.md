# AA-GenomicPrediction

Code and figures for manuscript "Subset-based genomic prediction provides insights into the genetic architecture of free amino acid levels in dry <i>Arabidopsis thaliana</i> seeds"

**snp_subset_pca.Rmd** - Runs PCA to compare structure across amino acid (AAS), control (CS), and genomic (GS) SNP sets

## GFBLUP Scripts
**01_extract_control_SNP_names.sh** - simple loop to extract marker names for control SNPs. Note that input data is not available on the git, but will be available as a supplementary file with the manuscript.

**02_prep_data_gfblup.Rmd** - cleans phenotypic and marker data and prepares inputs for GFBLUP, including genotype matrix (W) and genomic relationship matrix (G). Output is saved as .RData for future use.

**03_gfblup_AAS_only.Rmd** - runs standard GBLUP and GFBLUP for each amino acid trait; in GFBLUP only the AAS subset weighted. 

**Note:** scripts 4-6 are designed to run on a computing cluster using SLURM for job scheduling

**04_gfblup_AAS_vs_CS.R** - runs GFBLUP with both AAS and CS subsets weighted. Allows comparison for the relative proportion of genomic variance explained by the AAS and CS subsets. Exports TRUE/FALSE results testing for whether or not the proportion of genomic variance explained by AAS is greater than that of the CS subset. 
**04_gfblup_AAS_vs_CS_sub.sh** - submit file to run script 04 for all amino acid traits using sbatch

**05_gfblup_binom_test.R** - performs binomial test for the output of script 04
**05_gfblup_binom_sub.sh** - submit file for script 05

**06_gfblup_export_sigmas.R** - extracts values for the proportion of genomic variance explained for each run of GFBLUP with the AAS and CS subsets.
**06_gfblup_export_sigmas_sub.sh** - submit file for script 06

**07_plot_gfblup_results.Rmd** - compiles output from GFBLUP comparing the AAS and CS subsets. Plots results from standard GBLUP comparing predictive accuracy for the AAS and CS subsets and results from GFBLUP comparing the proportion of genomic variance explained by the AAS and CS subsets. See Figure 2 in the manuscript. 

## Helpful GFBLUP references:
**qgg R package**: https://github.com/psoerensen/qgg

Edwards SM, Thomsen B, Madsen P, Sørensen P. (2015) Partitioning of Genomic Variance Reveals Biological Pathways Associated with Udder Health and Milk Production Traits in Dairy Cattle. GSE https://doi.org/10.1186/s12711-015-0132-6

Edwards SM, Sørensen IF, Sarup P, Mackay TF, Sørensen P. (2016) Mapping Variants to Gene Ontology Categories Improves Genomic Prediction for Quantitative Traits in Drosophila melanogaster. Genetics https://doi.org/10.1534/genetics.116.187161

Sarup P, Jensen J, OstersenT, Henryon M, Sørensen P. (2016) Increased Prediction Accuracy using a Genomic Feature Model Including Prior Information on Quantitative Trait Locus Regions in Purebred Danish Duroc Pigs. BMC Genetics https://doi.org/10.1186/s12863-015-0322-9
