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
metagenome_pipeline.py -i ~/data/dss/functional_analysis/asv.dss.biom -m ~/data/dss/functional_analysis/16s_copy_numbers_and_nsti.tsv -f ~/data/dss/functional_analysis/EC_predicted_genom_DSS.tsv --max_nsti 1.99 -o ~/data/dss/functional_analysis/EC_metagenome1.99 --strat_out

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
