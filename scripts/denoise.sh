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
qiime dada2 denoise-paired --i-demultiplexed-seqs ~/data/dss/demuxed-dss.qza --p-trim-left-f 17 --p-trim-left-r 21 --p-trunc-len-f 260 --p-trunc-len-r 220 --o-table ~/data/dss/tableNoFilt.qza --o-representative-sequences ~/data/dss/repseqsNoFilt.qza --o-denoising-stats ~/data/dss/denoising-statsNoFilt.qza --p-n-threads 10

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
