#!/bin/bash
#SBATCH -p ghpc_v3
#SBATCH -N 1
#SBATCH --mem=120G
#SBATCH -t 4:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime feature-table filter-samples \
  --i-table ./table.qza \
  --p-min-frequency 6500 \
  --o-filtered-table ./table-6.5k.qza
  
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID