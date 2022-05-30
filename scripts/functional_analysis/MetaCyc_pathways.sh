#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=120G
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPTDIR
mkdir -p $TMPDIR

source activate picrust2
pathway_pipeline.py -i ~/data/dss/functional_analysis/EC_metagenome/pred_metagenome_contrib.tsv -o ~/data/dss/functional_analysis/EC_metagenome/functioMetaCyc_Pathways -p 10

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
