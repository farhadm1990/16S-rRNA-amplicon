#!/bin/bash
#SBATCH -p ghpc_v3
#SBATCH -N 1
#SBATCH -n 10  
#SBATCH --mem=64G
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
cd ~/data/dss && mkdir -p classifier && cd ./classifier && qiime feature-classifier classify-sklearn \
--i-reads ~/data/dss/repseqs.qza --i-classifier ~/data/dss/silva138-classifier-341f-805r.qza --o-classification \
 ~/data/dss/classifier/taxonomy-dss.qza --p-n-jobs 10 && qiime metadata tabulate --m-input-file ~/data/dss/classifier/taxonomy-dss.qza \
 --o-visualization ~/data/dss/classifier/taxonomy-dss.qzv

cd $SLURM_SUBMIT_DIR 
rm -rf /scratch/$USER/$SLURM_JOBID