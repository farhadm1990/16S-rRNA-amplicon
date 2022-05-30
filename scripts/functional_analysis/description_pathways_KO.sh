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
add_descriptions.py -i ./functional_analysis/KO_metagenome/functionKEGG_Pathways/path_abun_unstrat.tsv --custom_map_table ./functional_analysis/KO_metagenome/KEGG_pathways_info.tsv -o ./functional_analysis/KO_metagenome/functionKEGG_Pathways/path_abun_unstrat_descrip.tsv

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
