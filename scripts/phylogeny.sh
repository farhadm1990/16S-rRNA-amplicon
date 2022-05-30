#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=120G
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPTDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime fragment-insertion sepp --i-representative-sequences ~/data/dss/repseqsNoFilt.qza \
 --i-reference-database ~/data/dss/sepp-ref-gg-13-8.qza --o-tree ~/data/dss/treeNoFilt.qza \
 --o-placements ~/data/dss/tree-placementsNoFilt.qza --p-threads 10

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
