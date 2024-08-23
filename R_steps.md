After we created all artifacts in qiime2 such as ASV table, repseqs, rooted tree, taxonomy table and metadata, now we can bring them into R by usign qiime2_R package in R. But, we must first install and load the required packages:

```R
github.pkg <- c("jfukuyama/phyloseqGraphTest", "jbisanz/qiime2R") 
bioc.pkg <- c("phyloseq", "DESeq2", "decontam", "MicrobiotaProcess")
cran.pkg <- c("tidyverse", "pacman", "glue", "vegan", "devtools", "ggrepel", "reshape2", 
              "ggnetwork", "DT", "intergraph","VennDiagram", "lsmeans", "pheatmap", "phyloseqGraphTest")
inst.pkg <- cran.pkg %in% installed.packages()

 

if (any(!inst.pkg)){ install.packages(cran.pkg[!inst.pkg],repos = "http://cran.rstudio.com/") } 
 
inst.pkg <- github.pkg %in% installed.packages() 
if (any(!inst.pkg)){ devtools::install_github(github.pkg[!inst.pkg], force = TRUE) } 


inst.pkg <- bioc.pkg %in% installed.packages() 
if(any(!inst.pkg)){ BiocManager::install(bioc.pkg[!inst.pkg])
}
```


## 
## 1. Importing rooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
Using qiime2R package, we can bring all generated artifacts from qiime2 into R and integrate them into one phyloseq object by qiime_to_phyloseq()

```R
library("qiime2R")
library("phyloseq")

pst <- qza_to_phyloseq(features = "./tableNoFilt.qza", tree = "./treeNoFilt.qza", 
                       taxonomy = "./Taxonomy/taxonomy-dssNoFilt.qza", metadata = "./metadataNoFilt.tsv")
 
#we need to merge the refseqs like this since in above fun, there is no argument for that
repseqs <- read_qza("./repseqsNoFilt.qza")$data

pst <- merge_phyloseq(pst, repseqs)

  
```

Now `pst` is a `phyloseq` object containing all the artifacts and the metadata to our samples.
And if you want to see what does pst object contain, you can simply type `pst` and hit `ctrl/option(on mac)` + `enter` 

 ```R
pst
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 4540 taxa and 108 samples ]
sample_data() Sample Data:       [ 108 samples by 11 sample variables ]
tax_table()   Taxonomy Table:    [ 4540 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 4540 tips and 4520 internal nodes ]
refseq()      DNAStringSet:      [ 4540 reference sequences ]
```

If you want to revoke any artifacts from `pst` object, you can use the designated calling functons in phyloseq:

```R
otu_table(pst)    #To call ASV count table
tax_table(pst)    #To call taxonomy table 
phy_tree(pst)     #To call the phylogenetic tree
refseq(pst)       #To call the representative sequences for each ASV
sample_data(pst)  #To call the metadata
```

```R
otu_table(pst)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/ASV%20table%20in%20R.PNG)
> Figure 1. ASV count table. Rows are ASV ids and columns are the Sample IDs.
> 


```R
tax_table(pst)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/taxa%20table%20in%20R.PNG)
> Figure 2. Taxonomy table. Rows are ASV IDs and columns are taxonomic levels.
> 

##

```R
sample_data(pst)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Metadata%20in%20R.PNG)
> Figure 3. Metadata table. Rows are Sample IDs and columns are the variables.
> 


##
```R
refseq(pst)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Refseqs%20in%20R.PNG)
> Figure 4. Reference/representative sequence table. Width seq is for the neucleotide counts and names are ASV IDs.
> 

You can extract the sequences into your local directory as `fasta` format.

```R

#Using writeXStringSet
writeXStringSet(refseq(pst), "./dss.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
```


You can also convert the class of each variable in your metadata by the following code:

```R
#Converting all non-numeric/logical charatcters in the metadata into factors including pig_no, treatment, dss and gb

for(i in seq_len(ncol(sample_data(pst)))) {
 if(!is.numeric(sample_data(pst)[[i]]) && !is.logical(sample_data(pst)[[i]])) {
    sample_data(pst)[[i]] = as.factor(sample_data(pst)[[i]])  } else {
     sample_data(pst)[[i]]
 } 
}

#Optional: renaming the variable levels
#Treatment
iflese(sample_data(pst)$treatment == "ct", "CT", ifelse(sample_data(pst)$treatment == "gb", "GB",
ifelse(sample_data(pst)$treatment == "dss", "DSS", ifelse(sample_data(pst)$treatment == "gbdss", "GBDSS", "NegCtrl")))) #NegCtrl marks the negative samples.

#Segment
iflese(sample_data(pst)$sample_type == "digesta_25", "Proximal", ifelse(sample_data(pst)$sample_type == "digesta_50", "Midcolon",
ifelse(sample_data(pst)$sample_type == "digesta_75", "Distal", "Feces")))

# Releveling the treatment factors

sample_data(pst)$treatment<-factor(sample_data(pst)$treatment, levels = c('CT','GB','DSS', 'GBDSS','NegCtrl'))
```



**You can follow up more advanced bioinformatics and statistics on microbiome data on [Microloop](https://farhad1990.shinyapps.io/microloop)**
