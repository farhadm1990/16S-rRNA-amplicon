#!/bin/bash
#SBATCH -p ghpc_v3
#SBATCH -N 1
#SBATCH --mem=120G
#SBATCH -t 4:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime taxa barplot \
--i-table ./table-6.5k.qza \
--i-taxonomy ./taxonomy-dss.qza \
--m-metadata-file ./metadata.tsv \
--o-visualization ./taxa-barplot.qzv

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID