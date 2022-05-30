#!/bin/bash
#SBATCH -p ghpc_v3
#SBATCH -N 1
#SBATCH --mem=120G
#SBATCH -t 4:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime diversity beta-group-significance \
--i-distance-matrix ./core-metrics/unweighted_unifrac_distance_matrix.qza --m-metadata-file ./metadata.tsv --m-metadata-column treatment \
--o-visualization ./core-metrics/unw-unifrac-treatment.qzv --p-pairwise && qiime diversity beta-group-significance \
--i-distance-matrix ./core-metrics/weighted_unifrac_distance_matrix.qza --m-metadata-file ./metadata.tsv --m-metadata-column treatment \
--o-visualization ./core-metrics/w-unifrac-treatment.qzv --p-pairwise

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID