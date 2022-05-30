#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 10  
#SBATCH --mem=64G
#SBATCH -t 2:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime feature-classifier classify-sklearn \
--i-reads ~/data/dss/repseqsNoFilt.qza --i-classifier ~/data/dss/classifier/silva138-classifier-341f-805r.qza --o-classification \
 ~/data/dss/Taxonomy/taxonomyNoFilt.qza --p-n-jobs 10 && qiime metadata tabulate --m-input-file ~/data/dss/Taxonomy/taxonomyNoFilt.qza \
 --o-visualization ~/data/dss/Taxonomy/taxonomyNoFilt.qzv

cd $SLURM_SUBMIT_DIR 
rm -rf /scratch/$USER/$SLURM_JOBID