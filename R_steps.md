After we created all artifacts in qiime2 such as ASV table, repseqs, unrooted tree, taxonomy table and metadata, now we can bring them into R by usign qiime2_R package in R. 

## 1. Importing unrooted tree, ASV table and repseqs with the metadata to R using qiime2R package into R. 
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



```R
sample_data(pst)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Metadata%20in%20R.PNG)
> Figure 3. Metadata table. Rows are Sample IDs and columns are the variables.
> 



```R
refseq(pst)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Refseqs%20in%20R.PNG)
> Figure 4. Reference/representative sequence table. Width seq is for the neucleotide counts and names are ASV IDs.
> 


You can also convert the class of each variable in your metadata by the following code:

```R
#Converting all non-numeric/logical charatcters in the metadata into factors including pig_no, treatment, dss and gb

for(i in seq_len(ncol(sample_data(pst)))) {
 if(!is.numeric(sample_data(pst)[[i]]) && !is.logical(sample_data(pst)[[i]])) {
    sample_data(pst)[[i]] = as.factor(sample_data(pst)[[i]])  } else {
     sample_data(pst)[[i]]
 } 
}

#Optional: renaming the treatment levels
iflese(sample_data(pst)$treatment == "ct", "CT", ifelse(sample_data(pst)$treatment == "gb", "GB",
ifelse(sample_data(pst)$treatment == "dss", "DSS", ifelse(sample_data(pst)$treatment == "gbdss", "GBDSS", "NegCtrl"))
#NegCtrl marks the negative samples.

# Releveling the treatment factrs

sample_data(pst)$treatment<-factor(sample_data(pst)$treatment, levels = c('CT','GB','DSS', 'GBDSS','NegCtrl'))
```

## 2. Preprocessing and cleaning up the dataset

### Decontamination of the reads by `Decontam` package: unsupervised
It is based on the reads in the negative samples (here we have 5 neg samples) according to their prevalenc (not frequency), which is presence/absence of ASVs across samples.
The default threshold for a contaminant is that it reaches a probability of 0.1 in the statistical test being performed. In another word, an ASV should be present in at least 10% of the samples to pass the filter.  We can also try it with 0.5.

```R
library("decontam")

#First make a vector containing the name of the negative samples. 
neg_cont<-c("KDN1NK","KDN1N","KDN2N","KDNK","KDN3N")

contamPrev01 = isContaminant(pst, method = "prevalence", neg = "is.neg", threshold = 0.1)
contamPrev05 = isContaminant(pst, method = "prevalence", neg = "is.neg", threshold = 0.5)
```

Now you can extract the ID of the contaminant ASVs and remove the reads correlated to them.

```R 
contamASV = contamPrev01[contamPrev01$contaminant==TRUE, ] %>% rownames
#contamPrev05[contamPrev05$contaminant==TRUE, ]#this one is a bit strict!

#note that not all ASVs are contams (e.g. Escherichia, Pseudomonas, Corynebacterium) so we need to check them with the taxa
tax_table(pst)[rownames(tax_table(pst))%in%contamASV,]

#removing the contaminant ASVs from our phyloseq object
pst <- subset_taxa(pst, taxa_names(pst)!=contamASV)

#Now that we are done with the negative controls, we can remove them 
pst<- subset_samples(pst, !rownames(sample_data(pst))%in%neg_cont)
```

### Removing unasigned or NA reads from Phylum taxonomic level:

```R
pst <- subset_taxa(pst, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "unassigned"))

```


### Taxonomic filtering based on prevalence: supervsed

```R
#monitoring the number of the samples in which the prevalence of a taxon is at least one
prevdf <- apply(otu_table(pst),ifelse(taxa_are_rows(pst), 1, 2), function(x){sum(x>0)})
prevdf <-data.frame(ASVprev = prevdf,           #Number of the samples containing a given ASV; this shows how many samples each ASV occured.
                   TaxaAbund = taxa_sums(pst),  #The abundance of each ASV across all samples. 
                   tax_table(pst))
head(prevdf)
#Find out the phyla that are of mostly low-prevalence features by computing the total and average prev of features in each Phylum
plyr::ddply(prevdf, "Phylum", function(df){cbind(means = round(mean(df$ASVprev), 2), sums = round(sum(df$ASVprev),
                                         2))}) %>% mutate(sanity = ifelse(means == sums, "TRUE", "FALSE"))
#in the phylum level, Deinococcota and Myxococcota appeared only in one percent of samples and therefore, they will be filtered from the dataset
junkphyl = c("Myxococcota", "Deinococcota")
pst = subset_taxa(pst, !Phylum %in% junkphyl)

```

You could also do another level of filtering based on the prevalence/apearance of ASVs in a certain number of samples, e.g. 1 sample.
Warning: this might influence your alpha diversity measures, i.e. might cause underestimation of the changes as for instance richness is calculated based on the number of singletones in each sample. 

```R
#Filtering ASVs based on their prevalence threshold of 5 samples accross the samples. This means that each ASV should have appeared at least in 5 samples to be kept.
prevdf <- apply(otu_table(pst),ifelse(taxa_are_rows(pst), 1, 2), function(x){sum(x>0)}) #this shows the prevalence (if an ASV occured or not) of an ASV across the samples, therefore the range is from 1-108

table_count <- apply(otu_table(pst), 2, function(x) ifelse(x>0, 1, 0)) #this also gives the same results but with the count tables.
suspected_ASV = table_count[which((rowSums(table_count)/ncol(table_count))*100 < 5.4),] %>% rownames     #5.4 % of the samples equals 5 sample, since we have 108 samples here.

            
# you can derive the exact same results by the build-in phyloseq function 
# Prevalence: in more than 5 samples

TaxaTokeep <- genefilter_sample(pst.no.single, function(x) x>0, 5)                     

#let's remove these below-5 sample ASVs from the dataset
pst.prev = subset_taxa(pst, !taxa_names(pst)%in%suspected_ASV)
#pst.prev = prune_taxa(TaxaTokeep, pst.no.single)  # you could also use this code for filtering     
pst.prev   
pst.ancom=pst.no.single     
```


Since in this study we are only interested in the changes of bacteria depending on our environmental factors, we remove all non-bacterial ASVs in the Kingdom level.

```R
#now, we only keep the bacterial kingdom and remove the rest.
pst = subset_taxa(pst, Kingdom == "d__Bacteria")
pst
```


## 2. R-based analysis of microbiome data
When we imported all the artifacts from **qiime2** into **R**, we can use different packages and costume functions to render different *preprocessing* and *analitycal* steps.
