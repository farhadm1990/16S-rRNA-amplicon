# Microbiome Analysis

## Analysis of (gut) microbiome in Qiime2 and R.

In this module, I will walk you through the necessary steps involved in the analysis of 16S rRNA microbiome amplicone data from raw sequences to publication quality visualisations and statistical analysis. 
Note: all this workflow has been done on Jupyter notebook on a cluster node with 120 GB processer from Aarhus University, Denmark. In order to multitask in different nodes, tasks on Qiime2 have been submited to the cluster by seperate bash scripts.

This module includes the following steps:

## [1. Importing raw data into Qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Import.sh)
Here we import the raw sequences into a qiime artifact. But first you need to download the amplicon data [here](https://www.dropbox.com/scl/fo/o2j5uwiaynh6owsom9kf0/h?dl=0&rlkey=4qvl191j9zfx4332tfm4pul1k) and save it on your local drive. In order to bring the sequences from the local path (here they are saved on a folder in my cluster called 'amplicons', you must use this [bash command](https://github.com/farhadm1990/Microbiome_analysis/blob/main/scripts/import.sh). The input path must be defined by a [manifest file](https://github.com/farhadm1990/Microbiome_analysis/blob/main/manifestArranged.tsv), which includes the name of the sample and a path to each sample sequence for both forward (in one column) and reverse reads (in another column). Sequence data are paired end in the format of FASTA with quality score (Fastaq); therefore, in qiime2 the type will be "SampleData[PairedEndSequencesWithQuality]" and their imput format asigned as PairedEndFastqManifestPhred33V2.  
Here is the command:
```python
#!/bin/bash
#SBATCH -p ghpc_v3
#SBATCH -N 1
#SBATCH --mem=32G 
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path ~/data/dss/amplicons/manifestArranged.tsv \
  --output-path ~/data/dss/demuxed-dss.qza

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
```

## 2. Quality control and merging paired-end sequences (in Qiime2)
In this step we use DADA2 package to create the ASV tables and representative sequences (repseqs).

## 3. Training a primer region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)
In this step RESCRIPr will be used for creating more region specific, more sensitive based on our primerset.

## 4. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method
Using ASV table and repseqs we create a phylogenetic tree using SEPP package in Qiime2

## 5. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
Using qiime2R package, we can bring all generated artifacts from qiime2 into R and integrate them into one phyloseq object by qiime_to_phyloseq()

## 6. R-based analysis of microbiome data
When we imported all the artifacts from **qiime2** into **R**, we can use different packages and costume functions to render different *preprocessing* and *analitycal* steps.
