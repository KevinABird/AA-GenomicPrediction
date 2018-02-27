# AA-GenomicPrediction

Code and figures for manuscript "Subset-based genomic prediction provides insights into the genetic architecture of free amino acid levels in dry <i>Arabidopsis thaliana</i> seeds"

**BioMart_GeneList_Construction.Rmd** - Takes biosynthesis and catabolism csv files from GeneList_files directory to get start and stop position of genes usng biomaRt. Output should resemble At_aa_Biosynthesis_BiomaRt_DataSet.txt and At_aa_Catab_BiomaRt_Dataset 

**snp_subset_pca.Rmd** - Runs PCA to compare structure across amino acid (AAS), control (CS), and genomic (GS) SNP sets

## GBLUP Scripts

**01a_Make_AAS.Rmd** - use start and stop position from combined Biosynthesis and Catabolsm files output from BioMart_GeneList_Construction.Rmd to select SNPs from the whole SNP sets available in supplementary file.

**01b_Make_GS.Rmd** - Use start and stop positions from TAIR10 annotated genes to select sto select SNPs from the whole SNP sets available in supplementary file

**01c_MakeCS_1to50.R**- Recommended to run on a computing cluster and can potentially take over 24 hours to complete: Example script that iteratively makes 50 subsets by randomly selecting 335 genes to call SNPs, and comparing SNP count of randomly generated subset to AAS to ensure they're the same size. Recommend running several scripts in parallele to generate large numbers, e.g. 20 scripts making 50 subsets to generate 1000 subsets. Change value of *n* on line 23 and range in line 25.

**02a_MakeAAS-GS_Kin.R** - General GAPIT script to make a VanRaden kinship matrix from AAS or GS GAPIT input files.

**02b_MakeCSKin.R** - Script to make a directory for each CS subset and create a VanRaden kinship matrix, outputting in the newly made directory.

**Note:** Scripts 03a-03c are recommended to run on a computing cluster and can potentially take over 24 hours to complete:
**03a_AAS_GAPIT_1_50.R**- Sample script to run GAPIT 50 times on AA Subset for a single AA trait. Recommended to run several scripts in parallel changing the range in line 12: for(number in c(1:50)). To run on differen AA traits change line 20: AAabs[,c(1,2)] and line 57 colnames(AAabs)[2]) 

**03b_GS_GAPIT.R** -Sample script to run a single GAPIT run with 1000 replicates for a single AA trait. To run on differen AA traits change line 20: AAabs[,c(1,2)] and line 57 colnames(AAabs)[2]) 

**03c_CS_GAPIT_1_50.R** Sample script to run GAPIT on 50 different CS subsets for a single AA trait. Recommended to run several scripts in parallel changing the range in line 12: for(number in c(1:50)). To run on differen AA traits change line 20: AAabs[,c(1,2)] and line 57 colnames(AAabs)[2]) 

**04a_AAS_CSBinom.Rmd** Take GAPIT output from supplementary file and run a binomial test.

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
