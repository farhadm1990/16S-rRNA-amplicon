# ![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Cover.jpg)

#  Microbiome Analysis 


## Analysis of (gut) microbiome in **Qiime2** and **R**.

In this module, I will walk you through the necessary steps involved in the analysis of 16S rRNA microbiome amplicone data from raw sequences to publication quality visualisations and statistical analysis. Non-cultured 16S rRNA metagnomics is a promissing method for understanding the ecology of an environemnt in regards with the number and the structure of the microbiome in association with the environmental factors, e.g. host-microbiome interactions. In prokaryots there is a ubiquitous gene compartment integrated in the ribosome, so-called 16S rRNA genes, which are highly conserved among prokaryots and at the same time having hypervariable regions (HVRs) V1 to V9, which are good targets for evolutionary and ecological studies on prokaryotes [Jünemann et. al (2017)](https://www.sciencedirect.com/science/article/pii/S0168165617315985). This module is mainly focused on 16S rRNA gene data, but I can carefully say that you can apply most of the techniques explained here to genome data and count multivariate datasets.  
Note: all this workflow has been done on Jupyter notebook on a cluster node with 120 GB processer from Aarhus University, Denmark. In order to multitask in different nodes, tasks on Qiime2 have been submited to the cluster by seperate bash scripts.

This module includes the following steps

# [Steps in **Qiime**](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md)
### [1. Importing raw data into Qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#1-importing-raw-data-into-qiime2)
#
### [2. Filtering, dereplication, sample inference, chimera identification, and merging of paired-end reads by DADA2 package in qiime2](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#2-filtering-dereplication-sample-inference-chimera-identification-and-merging-of-paired-end-reads-by-dada2-package-in-qiime2)
#
### [3. Training a primer-based region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#3-training-a-primer-region-specific-classifier-for-taxonomic-classification-by-na%C3%AFve-bayes-method-in-qiime2)
#
### [4. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Qiime_Steps.md#4-creating-a-phylogenetic-tree-using-sate-enabled-phyhlogenetic-placement-sepp-method)
#

# [Steps in **R**](https://github.com/farhadm1990/Microbiome_analysis/blob/main/R_steps.md)

### [1. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R.](https://github.com/farhadm1990/Microbiome_analysis/blob/main/R_steps.md#1-importing-unrooted-tree-asv-table-and-repseqs-with-the-metadata-to-r-using-qiime2r-package-into-r) 
#
### [2. Preprocessing and cleaning up the sequences](https://github.com/farhadm1990/Microbiome_analysis/blob/main/R_steps.md#2-preprocessing-and-cleaning-up-the-dataset)
#
### [3. Alpha diversity](https://github.com/farhadm1990/Microbiome_analysis/blob/main/R_steps.md#3-alpha-diversity)
#
### [4. Statistical analysis on alpha diversity metrics.](https://github.com/farhadm1990/Microbiome_analysis/blob/main/R_steps.md#4-statistical-analysis-on-alpha-diversity-metrics-generalized-linear-mixed-effect-model-glmem)
#
### [5. Beta diversity.](https://github.com/farhadm1990/Microbiome_analysis/blob/main/R_steps.md#5-beta-diversity-diversity-between-samples)
#
### 6. Statistical analysis on beta diversity metrics.
