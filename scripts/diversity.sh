#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=120G
#SBATCH -t 2:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime diversity core-metrics-phylogenetic --i-table ~/data/dss/tableBac.qza --i-phylogeny ~/data/dss/treeBac.qza \
--m-metadata-file ~/data/dss/metadataFilt.tsv --p-sampling-depth 16000 --p-n-jobs-or-threads 10 \
--output-dir ~/data/dss/core_metricsBac16000

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID