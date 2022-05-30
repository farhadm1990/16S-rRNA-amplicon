#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1  
#SBATCH --mem=64G
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8 
qiime diversity alpha-rarefaction --i-table ~/data/dss/tableBac.qza \
--m-metadata-file ~/data/dss/metadataFilt.tsv \
--o-visualization ~/data/dss/alpha-rarefactionBac.qzv \
--p-min-depth 10 \
--p-max-depth 34400

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
