#!/bin/bash
#SBATCH -p ghpc_v3
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -t 1:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime feature-table tabulate-seqs --i-data ~/data/dss/repseqs.qza \
--o-visualization ~/data/dss/repseqs.qzv && qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

cd $SLURM_SUBMIT_DIR 
rm -rf /scratch/$USER/$SLURM_JOBID
