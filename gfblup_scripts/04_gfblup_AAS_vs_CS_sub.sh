#!/bin/bash -l
#SBATCH -D /home/turnersa/data/AA-GenomicPrediction/gfblup_scripts
#SBATCH -J gfblup
#SBATCH -o slurm-log/gfblup%j.out
#SBATCH -e slurm-log/gfblup%j.err
#SBATCH --mem 10000
#SBATCH --array=1-27

set -e
set -u

Rscript 04_gfblup_AAS_vs_CS.R $SLURM_ARRAY_TASK_ID
