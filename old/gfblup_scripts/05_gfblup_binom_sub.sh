#!/bin/bash -l
#SBATCH -D /home/turnersa/data/AA-GenomicPrediction/gfblup_scripts
#SBATCH -J binom_test
#SBATCH -o slurm-log/binom%j.out
#SBATCH -e slurm-log/binom%j.err
#SBATCH --mem 10000
#SBATCH --array=1-27

set -e
set -u

Rscript 05_gfblup_binom_test.R $SLURM_ARRAY_TASK_ID
