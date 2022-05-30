# Microbiome Analysis

## Analysis of (gut) microbiome in **Qiime2** and **R**.

In this module, I will walk you through the necessary steps involved in the analysis of 16S rRNA microbiome amplicone data from raw sequences to publication quality visualisations and statistical analysis. 
Note: all this workflow has been done on Jupyter notebook on a cluster node with 120 GB processer from Aarhus University, Denmark. In order to multitask in different nodes, tasks on Qiime2 have been submited to the cluster by seperate bash scripts.

This module includes the following steps

# [Steps in **Qiime**](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md)
### [1. Importing raw data into Qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#1-importing-raw-data-into-qiime2)
### [2. Filtering, dereplication, sample inference, chimera identification, and merging of paired-end reads by DADA2 package in qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#2-filtering-dereplication-sample-inference-chimera-identification-and-merging-of-paired-end-reads-by-dada2-package-in-qiime2)
### [3. Training a primer-based region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#3-training-a-primer-region-specific-classifier-for-taxonomic-classification-by-na%C3%AFve-bayes-method-in-qiime2)
### [4. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#4-creating-a-phylogenetic-tree-using-sate-enabled-phyhlogenetic-placement-sepp-method)


# [Steps in **R**](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Steps_R.md)

### 1. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
### 2. R-based analysis of microbiome data
