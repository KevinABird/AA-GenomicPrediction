#!/bin/bash -l
#SBATCH -D /home/turnersa/data/AA-GenomicPrediction/gfblup_scripts
#SBATCH -J export_sigmas
#SBATCH -o slurm-log/sigmas%j.out
#SBATCH -e slurm-log/sigmas%j.err
#SBATCH --mem 10000
#SBATCH --array=1-27

set -e
set -u

Rscript 06_gfblup_export_sigmas.R $SLURM_ARRAY_TASK_ID
