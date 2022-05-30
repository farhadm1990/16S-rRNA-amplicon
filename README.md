# Microbiome Analysis

## Analysis of (gut) microbiome in **Qiime2** and **R**.

In this module, I will walk you through the necessary steps involved in the analysis of 16S rRNA microbiome amplicone data from raw sequences to publication quality visualisations and statistical analysis. 
Note: all this workflow has been done on Jupyter notebook on a cluster node with 120 GB processer from Aarhus University, Denmark. In order to multitask in different nodes, tasks on Qiime2 have been submited to the cluster by seperate bash scripts.

This module includes the following steps

# Steps in **Qiime**


Before we continue on this jurney, I highly recommend you to check out the doc vignette of the [qiime tutorials](https://docs.qiime2.org/2022.2/tutorials/) related to the [installation](https://docs.qiime2.org/2022.2/install/native/) and case studies, such as [Parkinson's Mouse Tutorial](https://docs.qiime2.org/2022.2/tutorials/pd-mice/) for instance. 
## [1. Importing raw data into Qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Import.sh)
Here we import the raw sequences into a qiime artifact. But first you need to download the amplicon data [here](https://www.dropbox.com/scl/fo/o2j5uwiaynh6owsom9kf0/h?dl=0&rlkey=4qvl191j9zfx4332tfm4pul1k) and save it on your local drive. In order to bring the sequences from the local path (here they are saved on a folder in my cluster called 'amplicons', you must use this [bash command](https://github.com/farhadm1990/Microbiome_analysis/blob/main/scripts/import.sh). The input path must be defined by a [manifest file](https://github.com/farhadm1990/Microbiome_analysis/blob/main/manifestArranged.tsv), which includes the name of the sample and a path to each sample sequence for both forward (in one column) and reverse reads (in another column). Sequence data are paired end in the format of FASTA with quality score (Fastaq); therefore, in qiime2 the type will be "SampleData[PairedEndSequencesWithQuality]" and their imput format asigned as PairedEndFastqManifestPhred33V2.  
Here is the command:
```python
#!/bin/bash
#SBATCH -p ghpc_v3     #here is the name of the cluster node assigned for my work
#SBATCH -N 1           #number of the cores of processor on the cluster
#SBATCH --mem=32G      #required processor capacity
#SBATCH -t 24:00:00    #resrvation time

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path ~/data/dss/amplicons/manifestArranged.tsv \ #link to the folder in which my manifest file is located
  --output-path ~/data/dss/demuxed-dss.qza                 # link to the path where I want my demultiplexed data to be exported in

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
```

You can submit this bash script, which has a '.sh' format to the cluster by the following command:

```python
sbatch import.sh
```

This might take a while before you get the results. The output of this, is a `.qza` file that you have already specified it in the command, in this example `demuxed-dss.qza`. You can then create a visualized file from this artifact, with the following command:

```python
qiime demux summarize \
  --i-data ./demuxed-dss.qza \
  --o-visualization ./demuxed-dss.qzv

```
This `deuxed-dss.qzv`is a visualized format of our `demuxed-dss.qza`, which you can view it on [qiime2 viewer](https://view.qiime2.org/). Once you are there you can either drag-and-drop the [downloaded](https://github.com/farhadm1990/Microbiome_analysis/blob/main/artifacts/demuxed-dss.qzv?raw=true) artifact into the designated area or simpley copy the link to the artifact from [this repository](https://github.com/farhadm1990/Microbiome_analysis/blob/main/artifacts/demuxed-dss.qzv) and paste it in the box *file from the web*. Once there, you must come across the following picture:

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Demux.PNG)

On this `overview` page you can see counts of demultiplexed sequences for the entire samples for both forward and reverse reads, with min, median, mean and max and total counts. 

In the `Interactive Quality Plot` at the top left, you can see the quality (Phred, Q) score for the sequence bases in forward and reverse reads. 

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/demux-interactive.PNG)

Understanding this plot is very important for your denoising step, i.e. you must decide a truncation length for each forward and reverse reads in a way to keep at least 50% of the reads equal or above Q = 30. You can see this changes by hovering over the interactive box plots. You can see that in forward reads for example, the quality of reads starts falling in sequence baces of above 260 nt and already in 220 nt for reverse reads.

## 2. Filtering, dereplication, sample inference, chimera identification, and merging of paired-end reads by DADA2 package in qiime2.
In this step we use [DADA2](https://www.nature.com/articles/nmeth.3869) package to denoise the read and create the ASV tables and representative sequences (repseqs). 
Denoising is refered to sequence inference by resolving the sequencing erroros for each sequence as small as 1 nucleotide, which creates Amplicon Sequence Variants (ASVs). Denoising is superior over Operationals Taxonomic Units (OTU), which is clustering sequence reads into [OTUs](https://peerj.com/articles/5364/) at a defined cut-off rate for similarity (≥ 97), in a way that OTU clusters two different sequnces as the same OTU even if they are different in 3 nt. While, ASV can be as acurate as 1 nt difference. 
In denoising, we use 
## 3. Training a primer region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)
In this step RESCRIPr will be used for creating more region specific, more sensitive based on our primerset.

## 4. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method
Using ASV table and repseqs we create a phylogenetic tree using SEPP package in Qiime2

# Steps in **R**

## 5. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
Using qiime2R package, we can bring all generated artifacts from qiime2 into R and integrate them into one phyloseq object by qiime_to_phyloseq()

## 6. R-based analysis of microbiome data
When we imported all the artifacts from **qiime2** into **R**, we can use different packages and costume functions to render different *preprocessing* and *analitycal* steps.
