#!/bin/bash -l
#SBATCH -D /home/turnersa/data/AA-GenomicPrediction
#SBATCH -J gfblup
#SBATCH -o slurm-logs/gfblup%A_%a.out
#SBATCH -e slurm-logs/gfblup%A_%a.err
#SBATCH --array=1-53

set -e
set -u

Rscript 04_gfblup.R $SLURM_ARRAY_TASK_ID
