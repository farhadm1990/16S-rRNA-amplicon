# Microbiome_Analysis
=======
## Analysis of (gut) microbiome in Qiime2 and R.

In this module, I will walk you through the necessary steps involved in analysis of microbiome sequence data from raw sequences to publication quality graphs and visualizations. The module includes the following steps:

##1. Quality control and merging paired-end sequences (in Qiime2)
In this step we use DADA2 package to create the ASV tables and representative sequences (repseqs).

##2. Training a primer region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)
In this step RESCRIPr will be used for creating more region specific, more sensitive based on our primerset.

##3. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method
Using ASV table and repseqs we create a phylogenetic tree using SEPP package in Qiime2

##4. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
Using qiime2R package, we can bring all generated artifacts from qiime2 into R and integrate them into one phyloseq object by qiime_to_phyloseq()

