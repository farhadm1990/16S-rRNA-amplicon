# Microbiome Analysis

## Analysis of (gut) microbiome in Qiime2 and R.

In this module, I will walk you through the necessary steps involved in the analysis of 16S rRNA microbiome amplicone data from raw sequences to publication quality visualisations and statistical analysis. 
Note: all this workflow has been done on Jupyter notebook on a cluster node with 120 GB processer from Aarhus University, Denmark. In order to multitask in different nodes, tasks on Qiime2 have been submited to the cluster by seperate bash scripts.

This module includes the following steps:

## [1. Importing raw data into Qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Import.sh)
Here we import the raw sequences into a qiime artifact. Sequence data are paired end in the format of FASTA with quality score (Fastaq); therefore, in qiime2 the type will be "SampleData[PairedEndSequencesWithQuality]" and their imput format asigned as PairedEndFastqManifestPhred33V2. The input path must be defined by a manifest file which includes the name of the sample and a path to each sample sequence. 

## 2. Quality control and merging paired-end sequences (in Qiime2)
In this step we use DADA2 package to create the ASV tables and representative sequences (repseqs).

## 3. Training a primer region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)
In this step RESCRIPr will be used for creating more region specific, more sensitive based on our primerset.

## 4. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method
Using ASV table and repseqs we create a phylogenetic tree using SEPP package in Qiime2

## 5. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
Using qiime2R package, we can bring all generated artifacts from qiime2 into R and integrate them into one phyloseq object by qiime_to_phyloseq()

```R

ex = function(x){
x * 2
}
```
