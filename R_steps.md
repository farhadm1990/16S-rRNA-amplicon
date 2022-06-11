After we created all artifacts in qiime2 such as ASV table, repseqs, unrooted tree, taxonomy table and metadata, now we can bring them into R by usign qiime2_R package in R. 
##
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
ifelse(sample_data(pst)$treatment == "dss", "DSS", ifelse(sample_data(pst)$treatment == "gbdss", "GBDSS", "NegCtrl")) #NegCtrl marks the negative samples.

#Segment
iflese(sample_data(pst)$sample_type == "digesta_25", "Proximal", ifelse(sample_data(pst)$sample_type == "digesta_50", "Midcolon",
ifelse(sample_data(pst)$sample_type == "digesta_75", "Distal", "Feces")))


# Releveling the treatment factors

sample_data(pst)$treatment<-factor(sample_data(pst)$treatment, levels = c('CT','GB','DSS', 'GBDSS','NegCtrl'))
```
# 

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

table_count <- apply(otu_table(pst), 2, function(x) ifelse(x>0, 1, 0)) #this also gives the same results.
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
#
### Removing singletones based on abundance: Supervised
#
#### Mximum threshold abundance for ASVs in the count table

Singletones are those ASVs that occure only once across all samples. Therefore, their relative abundance `[ASV/sum(ASVs)]` equals to `1`. Since we have 4 different gourps of treatments, an ASV with appearance in only one sample cannot be due to biological differences, otherwise we would expect it to appear in at least one group of treatments. If you have done the prevoius 'prevalence filteing' step, there must be no singletones left in the dataset. However, removing singletones are not recommended for estimation of beta diversity metrics, especially richness as it is based on singletones and doubletones. While, you can remove them for differentially abundance analysis to get a clear signal from the dataset.

First make a relative abundance from the ASV count table and see which ASVs occured only once with the following function. 

```R
#A function to find singletones You need to be careful about this step!
out.ASV = function (phyloseq, threshold =1, binwidth = 0.01) {
       
#This function requires phyloseq, tidyverse and glue packages to be loaded. 
    
                      rel_abund = t(apply(otu_table(phyloseq), ifelse(taxa_are_rows(phyloseq), 1,2), 
                      function(x) x/sum(x))) #making the relative abundance table
                      names.single = t(apply(rel_abund, 1, function(x){ifelse(x == threshold, TRUE, ifelse(x == sum(x),
                      TRUE, FALSE))})) %>% melt %>% filter(value == TRUE) %>% select(1) %>% pull  %>% unfactor()
                      single.ASV = rel_abund[rownames(rel_abund)%in%as.vector(names.single),]
                      single.ASV[single.ASV == 0] <- NA # A seperate dataset for annotation of singletones on the barplot
                        
                        if (length(names.single) == 0 ) {
                        print(glue("{length(names.single)} singletones detected in this dataset"))
                        qplot.noSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                        color = rel_abund, main = "Frequency count of relative abundance") +
                        xlab ("Relative abundance in samples") + ylab("Frequency")
                            
                        
                         return(structure(list(qplot.noSing, names.single))   )
                            
                        } else {
                                         
                       qplot.withSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, color = rel_abund, 
                       main = "Frequency count of relative abundance with singletones") +
                       geom_bar(aes(single.ASV),  color = "red", width = binwidth)+
                       xlab ("Relative abundance in samples") + ylab("Frequency") + 
                       geom_label_repel(aes(x = 0.95, y =length(rel_abund)/10), 
                       label.padding =  unit(0.55, "lines"), label = glue("{length(names.single)}\n Singletones"), color = "black")
                            
                       qplot.rmSing = qplot(rel_abund[!rownames(rel_abund) %in% names.single, ], geom = "histogram",
                       binwidth = binwidth, main = "Frequency count of relative abundance without singletones") +
                       xlab ("Relative abundance in samples") + ylab("Frequency")
                            
                       print(glue('{length(names.single)} singletones detected in the dataset'))
                       return(structure(list(qplot.withSing, qplot.rmSing, unlist(names.single))) )
                    
                        }                        
    
                             
        }
                        
single.test = out.ASV(phyloseq = pst, threshold = 1, binwidth = 0.1)
singletones = single.test[[3]] #here you can extract the names of the singletones and remvoe them from the dataset.

#Now you can remove the singletones from your pst file as follows:
pst.no.single = subset_taxa(pst, !taxa_names(pst)%in% singletones)
pst = pst.no.single
```
The outcome of this function is a barplot showing the relative abundance of taxa on `x` axis and their counts across the count table on the `y` axis. If there is not singletones, the barplot should show no counts on relative abundance of `1` on `x` axis. 


```R
single.test[[1]]
```

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Singletone_plot.PNG)
> Figure 5. Barplot of ASV relative abundance with their frequency across ASV table. The red bar represents the count of singletones.
> 

#

#### Taxonomic filtering based on minimum threshold abundance 
Unlike the prevoius step, which we removed the ones with outmost abundance `1` or `100%`, here, we look for those with a minimum threshold abundacne, i.e. an ASV should occure in > 0.01% across all samples to pass the filter.

```R
# Abumdance: ASV > 0.01% overall abundance across all samples
total.depth <- sum(otu_table(pst.prev))
totAbuThreshold <- 1e-4 * total.depth
pst <- prune_taxa(taxa_sums(pst.prev)>totAbuThreshold, pst.prev)

```

Since in this study we are only interested in the changes of bacteria depending on our environmental factors, we remove all non-bacterial ASVs in the Kingdom level.

```R
#now, we only keep the bacterial kingdom and remove the rest.
pst = subset_taxa(pst, Kingdom == "d__Bacteria")
pst
```

#

### Rarefaction to an equal sampling depth

Due to the technical reasons, Next Generation Sequencing (NGS) machines like Illumina do not usually run the squencing equally for all samples. This results in an unequal sampling/read depth among samples. Therefore, it is advisable to nurmalize the data to an equal sampling depth. Rarefaction is resampling without replacement, which has been widely used for this type of normalization. However, depending on different pipelines for denoising or clustering, e.g. be it DADA2 or Deblur, new concerns have risen in regard with rarefaction as it might partially affect our estimate of richness and consequently the biological interpertation of microbiota ([Bardenhorst et. al, 2022](https://www.sciencedirect.com/science/article/pii/S2001037021005456)).

You can do rarefaction using `rarefy_even_depth` funciton from `rarefy_even_depth` package. But, in order to get an understanding on the sampling depth and consequent changes in the reads, you can create a rarefaction curve and choose the depth from which on you don't see an increase in the richness or diversity by increasing the reading depth. 

```R
library(MicrobiotaProcess)

ps_rar_curve <- ggrarecurve(obj = pst, 
                     indexNames = c("Observe", "Shannon"),
                     chunks = 400, 
                     theme(legend.spacing.y = unit(0.02, "cm"),
                           legend.text = element_text(size = 6)), show.legend=F)

ps_rar_curve
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/rarefaction%20curve.PNG)
> Figure 6. Rarefaction curve for richness and Shannon indices. You can see that around 17000 sampling depth the curves reach platue.


With the following code you can perform rarefaction:

```R
pst_rar_repF <- phyloseq::rarefy_even_depth(pst, sample.size = 17000, replace = F) #without replacement and 17000 sampling depth

```
# 

### Converting reads from relative abundance to absolute abundance

For each sample we extracted total load of 16S rRNA copy-number genes by perfomring a qPCR and we use them to convert the relative abundnce counts into aboslute abundance. Since we only did that for samples from proximal, distal and feces, we prune our samples to these types. Then we convert the counts into relative abundance and multiply them to the total load of 16S rRNA genes for each correspondent samples. 

```R
#Prunning the samples

pst.qPCR = prune_samples(sample_data(pst.rar.clean)$sample_type %in% c("Proximal", "Distal", "Feces"), pst_rar_repF) 

# Converting the observed counts to relative abundant.
pst.qPCR.relab <- transform_sample_counts(physeq = pst.qPCR, fun = function(x) x/sum(x))

 # Normalizing the counts into their absolute abundance
pst.abs.rel = t(otu_table(pst.qPCR.relab)) * sample_data(pst.qPCR)$copy_number
pst.abs.rel = t(pst.abs.rel)                                          
pst.abs.rel = apply(pst.abs.rel, 2, function(x) round(x))

otu_table(pst.qPCR) <- otu_table(pst.abs.rel, taxa_are_rows = TRUE) 

# You can also make a copy of your data as log-transformed

pst.qPCR.log = transform_sample_counts(pst.qPCR, function(x) {log(x+1)})
```
# 

## 3. Alpha diversity

Estimating alpha diversity, diversity within sample, using `estimate_richness` function of `phyloseq` package.

```R

#Calculating the alpha diversity indexes for qPCR data:
#Richness is the number of observed ASVs
 
Chao1 =estimate_richness(pst, split = TRUE, measures = "Chao1") #for richness, we don't use rarefied table

#Diversity: takes both observed and evenness into account in a  way that it shows that if the community has similar abundances. Therefore, a community with a high diversity has a lot of species with similar abundances. A community with many species but only a single dominant one has a lower diversity. Shannon is an index for diversity measurement and it can be 4 or hihger in value. 
Shannon = estimate_richness(pst.qPCR, split = TRUE, measures = "Shannon")

#Evenness: it is the probability of if two bacteria belong to the same speices. In other word, how even the abundance of bacteria are. since it is probabilility, the value for evenness index, e.g. Pielou, is always between 0-1.
normcount = apply(otu_table(pst.qPCR), 2, function(x) x/sum(x))#we need a relative abundance count to measure Pielou
ASVcount = colSums(normcount != 0)                  
Pielou = Shannon / log(ASVcount)     #pielou corrects the shannon index for speices number             
Simpson = estimate_richness(pst.qPCR, split = TRUE, measures = "Simpson")                  

#Phylogenetic distance                   
library(picante)
FaithPD = pd(t(otu_table(pst.qPCR)), tree = phy_tree(pst.qPCR), include.root = F)$PD
                  
#Adding the indexes to the metadatas              
sample_data(pst.qPCR) <- data.frame(sample_data(pst.qPCR), Chao1=Chao1[[1]],
                                    Shannon = Shannon$Shannon,  FaithPD = FaithPD)   
sample_data(pst.qPCR)   


```

#

## Creating a stacked barplot for visualizating the compositon of microbiota in different taxonomic level.

Since we can create stacked barplot for different taxon levels and the ASV IDs are rather complicated and long to read, you can use the following function to agglomerate your taxa into a certain level by giving them the level name as the row names.

```R
#A function to create unique names for each ASV. It removes any NA in Order level then attempts to use the name of one level higher taxa for those 
#who have similar names, e.g. uncultured_bacterium

gloomer = function(ps = data, taxa_level = taxa_level, NArm = "TRUE"){
    rank.names = c('Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    

#====================Sometimes in genus level, we might have multiple uncultured organisms, which if we want to make unique out of them for the species level it won't work====
    #since adding uncultured to uncultered is sill duplication. therefore if the taxa_level is set to species we first make a unique genus and then we go further to the speices===#

#Removing unculured Family
ps = subset_taxa(ps, !Family %in% c("uncultured", "NA"))
    
if(taxa_level == "Species") {
ps = subset_taxa(ps, !Genus %in% NA)#we remove genus tagged NA
tax_table(ps)[, taxa_level] <- ifelse(tax_table(ps)[, taxa_level] %in% NA, paste0("unknown"), paste(tax_table(ps)[, taxa_level]))#convert NA in species into unknown
    
  physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)
    
#first take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] == "uncultured", 
       paste0(taxdat[ , length(rank.names[1:which(rank.names=="Genus")])-1], "_", taxdat[,6]), paste(taxdat[,6]))

spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% melt %>% filter(value == "TRUE") 

if(dim(duplis)[[1]] > 0) {
duplis = uni %>% melt %>% filter(value == "TRUE") %>% dplyr::select(1) %>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis, 
                    paste0(taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "_", taxdat[,taxa_level]), taxdat[,taxa_level]))

    taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as.matrix(taxdat)   
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
    
} else {
    
taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
    
}
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))

    
    
#==========================================# 
} else if (taxa_level == "Genus") {
    
    physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)
    
# take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] == "uncultured", 
       paste(taxdat[ , length(rank.names[1:which(rank.names=="Genus")])-1], "_", taxdat[,6]), paste(taxdat[,6]))

  
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[taxdat[,taxa_level] %in% rownames(otudat), taxa_level]
taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))
 
    
    
} else {
    
    
physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = TRUE)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
otudat = otu_table(physeq)
    
spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE")

if(dim(duplis)[[1]] > 0){#if there is not duplications, we can simply use the taxa names as the row name
    
    duplis = uni %>% reshape2::melt %>% filter(value == "TRUE") %>% dplyr::select(1)%>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis, 
                    paste(taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "_", taxdat[,taxa_level]), taxdat[,taxa_level]))
    
taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as.matrix(taxdat)   
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
} else {

taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
}
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))
 

}
return(physeq) 
    }
    
 ```
 
 Now we can use the `gloomer` function to agglomerate our reads to `Phylum` level and create the stacked barplot for the absolute counts.
 
```R
#Making a costumized color list for the phylum
phylcol=c('coral4', "darkorange",'antiquewhite4','gold', 'cornflowerblue', 'plum4',
          'darkgoldenrod3','aquamarine4', 'cadetblue2', 'darkgreen', 'mediumseagreen', 'red','Gray',
        'steelblue2','darkmagenta')

#Agglomeration
glom.phyl <- gloomer(pst.qPCR, "Phylum", NArm = TRUE)

#Merging counts for treatments
trans.ps1 <- merge_samples(glom.phyl, "treatment") 

                                     
#Barplot                                    
plot_bar(trans.ps1, fill="Phylum", x = "treatment2") + scale_fill_manual(values = phylcol) + 
    xlab("Treatments") + ylab(" 16S rRNA gene copies per g") + 
scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+#to print the scale in scientific notation
            theme_bw() + 
            theme(text = element_text(size =15, face = "bold"))
                   
ggsave("./barplot_phylum_absolute.jpeg", device = "jpeg", dpi = 300, width = 7) 
```
 
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/barplot_phylum_absolute.jpeg) 
 > Figure 7. Stacked barplot for different phyla in 4 treatment groups.

#

## Creating a venn diagram for shared taxa between different treatments.
#
### Species level

```R
library(eulerr)
library(gplots)
library(VennDiagram)

spec.pst <- gloomer(pst.qPCR, taxa_level = "Species", NArm = TRUE)

#making a subset of dataset for different treatments
ct.ps <- subset_samples(spec.pst, treatment == "CT")
gb.ps <- subset_samples(spec.pst, treatment == "GB")
dss.ps <- subset_samples(spec.pst, treatment == "DSS")
gbdss.ps <- subset_samples(spec.pst, treatment == "GBDSS")

#filtering data for counts bigger than 0
asv.ct <- round(otu_table(ct.ps)[apply(otu_table(ct.ps), 1, function(x) any(x > 0)),],0) %>% rownames()
                                      
asv.gb <- round(otu_table(gb.ps)[apply(otu_table(gb.ps), 1, function(x) any(x > 0)),],0)%>% rownames
                                       
asv.dss <- round(otu_table(dss.ps)[apply(otu_table(dss.ps), 1, function(x) any(x > 0)),],0) %>% rownames()
                                        
asv.gbdss <- round(otu_table(gbdss.ps)[apply(otu_table(gbdss.ps), 1, function(x) any(x > 0)),],0)  %>% rownames()
                                             
#First making a simple venn diagram to see the actual counts and unions of ASVs for each treatment
gplots::venn(data =list(asv.ct, asv.gb, asv.dss, asv.gbdss)) 

#Then we use the numbers to make a prettier venn diagram 
ven.diag = draw.quad.venn(col = "white", alpha = 0.75, area1 = 0+23+3+2+3+167+0+17, area2 = 23+1+3+1+167+16+17+1, area3 = 2+3+167+3+16+1+1+14, area4 = 3+0+17+167+1+16+14+1,
                          fontface = "bold", 
               n12 = 23+3+167+17, n13 = 2+3+3+167, n14 = 3+0+167+17, n1234 = 167, n123 = 3+167, n124 = 167+17, n134 =3+167,
              n234 = 167+16, n23=3+1+167+16, n24=167+16+17+1, n34 = 3+167+16+14, category = c("CT", "GB", "DSS", "GBDSS"), 
              fill = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4"))
ggsave(plot = ven.diag, "./venn.taxa.species.jpeg", device = "jpeg", dpi = 500)

```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/venn.taxa.species.jpeg)
> Figure 8. Venn diagram of shared species between treatments.

#
### Phylum level

```R
#Phylum level
phyl.pst <- gloomer(pst.qPCR, taxa_level = "Phylum", NArm = TRUE)

#making a dataset for different treatments
ct.ps <- subset_samples(phyl.pst, treatment == "CT")
gb.ps <- subset_samples(phyl.pst, treatment == "GB")
dss.ps <- subset_samples(phyl.pst, treatment == "DSS")
gbdss.ps <- subset_samples(phyl.pst, treatment == "GBDSS")

#filtering data for counts bigger than 0
asv.ct <- round(otu_table(ct.ps)[apply(otu_table(ct.ps), 1, function(x) any(x > 0)),],0) %>% rownames()
                                      
asv.gb <- round(otu_table(gb.ps)[apply(otu_table(gb.ps), 1, function(x) any(x > 0)),],0)%>% rownames
                                       
asv.dss <- round(otu_table(dss.ps)[apply(otu_table(dss.ps), 1, function(x) any(x > 0)),],0) %>% rownames()
                                        
asv.gbdss <- round(otu_table(gbdss.ps)[apply(otu_table(gbdss.ps), 1, function(x) any(x > 0)),],0)  %>% rownames()
                                             
#First making a simple venn diagram to see the actual counts and unions of ASVs for each treatment
gplots::venn(data =list(asv.ct, asv.gb, asv.dss, asv.gbdss)) 

#Then we use the numbers to make a prettier venn diagram 
ven.diag = draw.quad.venn(col = "white", alpha = 0.75, area1 = 0+1+0+0+10+0+0+2, area2 = 1+1+0+0+10+1+2+0, area3 = 0+0+0+10+0+1+0+0, fontface = "bold",
                area4 = 0+0+10+2+1+0+0+0, 
               n12 = 1+0+10+2, n13 = 0+0+0+10, n14 = 0+0+2+10, n1234 = 10, n123 = 0+10, n124 = 10+2, n134 =10+0,
              n234 = 10+1, n23=0+0+10+1, n24=10+1+0+2, n34 = 0+10+1+0, category = c("CT", "GB", "DSS", "GBDSS"), 
              fill = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4"))
ggsave(plot = ven.diag, "./venn.taxa.phylum.jpeg", device = "jpeg", dpi = 500)                                            
```

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/venn.taxa.phylum.jpeg)
> Figure 9. Venn diagram of shared phylum between treatments.


# 

## 4. Statistical analysis on alpha diversity metrics: Generalized Linear Mixed Effect Model (GLMEM)

This experiment was done in a 3 complete blocks with a `2 x 2` factorial design for the treatment arrangement on 23 6-week old piglets. In each block/round, we have used 2 litters (sow from which we got the piglets). After `19` days of experiemnt, piglets were sacrificed and digesta samples were taken from their proximal and distal colon and feces. In each round of experiment (n = 3), there were `4` pens in which 2 pigs were housed.

```R

#Loading required packages
library(nlme)
library(lme4)
library(lmerTest)
library(car)
library(lsmeans)
library(postHoc)
library(multcomp)

# Creating a data frame from the metadata.
alpha.qpcr = sample_data(pst.qPCR) %>% data.frame

# Inspecting the distribution of alpha metrics

alpha.qpcr$Chao1 %>% hist(main = "Chao")
alpha.qpcr$Shannon  %>% hist(main = "Shannon")
alpha.qpcr$FaithPD %>% hist(main = "FaithPD")

# Checking the experiment desing and rank efficiency of variables
xtabs(~ gb + dss + litter +  Pen, exp.test)

# Fitting the data to a GLMEM by glmer function with Gamma distribution and log link. We can use the optimizer control = glmerControl(optimizer = c("bobyqa", "bobyqa")) to fix convergance issue, if there was any.
faith.lm = glmer(FaithPD ~ gb * dss + sample_type + blok  + (1|litter) + (1|pig_no), 
            data = alpha.qpcr, family = Gamma(link = "log")) # Full model
            
faith.red = glmer(FaithPD ~ gb * dss   + (1|litter) + (1|pig_no), 
            data = alpha.qpcr, family = Gamma(link = "log")) # Reduced model
anova(faith.lm, faith.red)                                   # Checking for differences in the models after reduction

car::Anova(faith.red, type = 2) #A qi-square test. if the interaction is not significant use the default type 2 and otherwise type 3 (reports the intercept). 

faith.pairs = pairs(emmeans(faith.lm, ~ gb + dss), type = "response", adjust = "BH") %>% cld(Letters = letters) # Extracting responses and their pairwise comparison. Type="response" can give you a back transformed pair-wise comparision of EMMs. 

#contrast(emmeans(chao.red, ~ gb * dss), adjust = "BH", type = "response", interaction =  T)

cld(emmeans(faith.lm, ~ gb + dss, adjust = "BH",  type = "response"),Letters = letters) 
#summary(emmeans(chao.red, ~ gb + dss, adjust = "BH",  type = "response"))

```

Now you can extract the parwise comparisons from the model to be used for figure annotations.

```R
#extracted pairwise comparisons from the model

chao.pairs = chao.pairs %>% as.data.frame %>% rownames_to_column("order") %>% arrange(order)
shannon.pairs = shan.pairs %>% as.data.frame %>% rownames_to_column("order") %>% arrange(order)
faith.pairs = faith.pairs %>% as.data.frame %>% rownames_to_column("order") %>% arrange(order)


#creating df for the manual pvalue
chao.stat = chao.pairs %>% mutate(group1 = factor(c("CT", "CT", "CT", "GB", "GB", "DSS"), levels = c("CT", "GB", "DSS") ), 
                       group2 = factor( c("GB", "DSS", "GBDSS", "DSS", "GBDSS", "GBDSS"), levels = c("GB", "DSS", "GBDSS" )),
                       p =ifelse(chao.pairs$p.value < 0.06 & chao.pairs$p.value > 0.01 , "*", 
                           ifelse(chao.pairs$p.value < 0.01 &  chao.pairs$p.value > 0.001, "**",
                           ifelse(chao.pairs$p.value < 0.001 , "***", "ns")))) %>% as_tibble

shannon.stat =  shannon.pairs %>% mutate(group1 = factor(c("CT", "CT", "CT", "GB", "GB", "DSS"), levels = c("CT", "GB", "DSS") ), 
                       group2 = factor( c("GB", "DSS", "GBDSS", "DSS", "GBDSS", "GBDSS"), levels = c("GB", "DSS", "GBDSS" )),
                       p =ifelse(shannon.pairs$p.value < 0.06 & shannon.pairs$p.value > 0.01 , "*", 
                           ifelse(shannon.pairs$p.value < 0.01 &  shannon.pairs$p.value > 0.001, "**",
                           ifelse(shannon.pairs$p.value < 0.001 , "***", "ns")))) %>% as_tibble

faith.stat =  faith.pairs %>% mutate(group1 = factor(c("CT", "CT", "CT", "GB", "GB", "DSS"), levels = c("CT", "GB", "DSS") ), 
                       group2 = factor( c("GB", "DSS", "GBDSS", "DSS", "GBDSS", "GBDSS"), levels = c("GB", "DSS", "GBDSS" )),
                       p =ifelse(faith.pairs$p.value < 0.06 & faith.pairs$p.value > 0.01 , "*", 
                           ifelse(faith.pairs$p.value < 0.01 &  faith.pairs$p.value > 0.001, "**",
                           ifelse(faith.pairs$p.value < 0.001 , "***", "ns")))) %>% as_tibble
```

Making a violin plot for FaithPD, as an example, with pairwise comparisons between treatments. 

 

```R
library(ggpubr)

#Creating a long dataframe
long_mtdat = long_data[long_data$variable %in% c("FaithPD"),]

ggboxplot(long_mtdat, x = "treatment", y = "value", width = 0, facet.by = "variable") +
geom_violin( aes(fill = treatment), alpha=0.7, trim = FALSE)+ 
geom_boxplot(width =0.15, outlier.color = NA)+ stat_pvalue_manual( data = faith.stat, y.position = 35.8,
xmin = "group1", xmax = "group2", step.increase =0.1, label = "p")  +
geom_jitter(aes(color = treatment), alpha = 0.4, size = 2, show.legend = F) +
scale_fill_manual(values = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4"))+ 
scale_color_manual(values = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4"))+
labs(fill= "Treatment", title = "Violin plot of FaithPD based on absolute counts", y = "Alpha diversity",
x = "Treatment") + theme_bw() + 
theme(text = element_text(size =11, face = "bold"))

ggsave(filename = "./alpha_dss_abs_faith.jpeg", device = "jpeg", dpi = 300, width = 5)
```
![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/alpha_dss_abs_faith.jpeg)
> Figure 10. Violin plot of FaithPD metric for four treatments. Pairwise comparisons were adjusted by BH at P < 0.5.

#

## 5. Beta diversity: diversity between samples

Before measuring beta diveristy, it is worth to know that unlike alpha diversity which deals with differences of diversity within samples, beta diversity accounts for diversity and differences in community between groups of samples. In this chapter we `phyloseq` and `vegan` packages. There are different beta diversity metrics such as Jaccard index, Bray-Curtis dissimilarity or UniFrac distance.  

#### The Jaccard distance: 
To calculate the Jaccard index, we take the intersection (the number of species present in both samples) and divide this 
by the total number of species, or their union. 

#### Bray-Curtis dissimilarity:
Bray-courtis, unlike the Jaccard (only taking presence absence into account), accounts for the number of features in each groups (evenness). Therefore, it is to say that Bray-curtis is the beta equivalent of Shannon in alpha. To measure Breay-Curtis dissimilarity matrix, we first look at the abundances, or the total counts, of species shared between the two samples. Next, we take the smallest value for each of the shared species. We add all these smallest values and multiply it by two. Next, we divide this by the sum of the values in both samples. This dissimilarity matrix has range of 0 (similar) to 1 (dissimilar).

#### UniFrac Phylogenetic Distance
The UniFrac distance does not use species directly. Instead, it works with a phylogenetic tree. The phylogenetic tree is made from differences between genes. Species on branches close together in the tree have more similar genes. Instead of counting species present in both samples, we sum the lengths of branches that are not shared and divide them by the sum of all branch lengths (unweighted UniFrac) and if you multiply the branch lenghts to the abundance of taxa it will be called weighted UniFrac.


To estimate beta diversity, it is advisable to `log` transform your data to account for the zero inflation. Here we literaly add a psodocount to the zero counts as log of zero is undefined.

```R
pst.qPCR.log <- transform_sample_counts(pst.qPCR, function(x) log(1+x))

```

If you want to know different distance metrics, you can do as follows:
```R
#Here is a list of all distance metrics
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
```

You can create a metric distance scaling (MDS) also termed as Principal Coordinate Analysis (PCoA) or non-metric distance scaling (NMDS) plots to visualize beta diversity of your bacterial data on reduced dimentional plots. In this example we only do that for Weighted UniFrac Phylogenetic Distance by `ordinate` function from phyloseq package.

```R
#Weighted UniFrac: based on log-transformed data 

#Wunifrac PCoA
wunifrac.pcoa = ordinate(pst.qPCR.log, method = "PCoA", distance = "wunifrac")
evals<-wunifrac.pcoa$values$Eigenvalues

wunifrac.pcoa.log = plot_ordination(pst.qPCR.log, wunifrac.pcoa, color="treatment", shape= "sample_type", 
   title = "PCoA plot of log Weighted UniFrac")+
   labs(col="Treatment", shape = "Segment")+
   coord_fixed(sqrt(evals[2]/evals[1]))+  #+stat_ellipse() you could also add ellipse to your clusters.
   labs(x = sprintf("PCoA1 [%s%%]", round(evals/sum(evals)*100,1)[1]),
   y = sprintf("PCoA2 [%s%%]", round(evals/sum(evals)*100, 2)[2])) +
   scale_color_manual(values = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4")) + 
   theme_classic() + geom_vline(xintercept = 0, lty = 2, alpha=0.5) + 
   geom_hline(yintercept = 0, lty = 2, alpha =0.5)  + geom_point(size = 4) +coord_fixed()+
   scale_y_continuous(limits = c(-0.12, 0.15)) + scale_x_continuous(limits = c(-0.15, 0.21)) +  theme_bw() +
   theme(text = element_text(size = 15, face = "bold"))
wunifrac.pcoa.log
```
![(wunifrac.pcoa.log)](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/wunifrac.pcoa.log.jpeg)
> Figure 11. PCoA plot for weighted UniFrac distance. Each point represents one sample with different colors correspondent to treatments and the segment/sample type as the shape for the points. 

```R
#Wunifrac NMDS
wunifrac.nmds.log=ordinate(pst.qPCR.log, method="NMDS", distance = "wunifrac")
wunifrac.stress = stressplot(wunifrac.nmds.log)#in this dataset, R squared shows that the nMDS is perfectly able to capture variation in the data.
# plot(wunifrac.nmds.log)#the red points are Transformed taxon-wise dissimilarities and the circles are Transformed sample-wise dissimilarities.

wunifrac.nmds.log = plot_ordination(pst.qPCR.log, wunifrac.nmds.log, color="treatment", shape= "sample_type", 
  title = "Weighted UniFrac distance NMDS plot of log data")+
  labs(col="Treatment", shape = "Segment")+
  scale_color_manual(values = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4")) +
  theme_bw() + geom_vline(xintercept = 0, lty = 2, alpha=0.5) + 
  geom_hline(yintercept = 0, lty = 2, alpha =0.5)+scale_x_continuous(limits = c(-0.15, 0.25)) + 
  scale_y_continuous(limits = c(-0.1, 0.1)) + theme(text = element_text(size = 15, face = "bold")) + geom_point(size =4)+
  coord_fixed()
wunifrac.nmds.log 
```
![stress plot](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/stress%20plot%20of%20Wnweighted%20UniFrac%20distance%20for%20NMDS%20plot%20with%20log%20transformation.jpeg)
> Figure 12. Stress plot for NMDS of weighted UniFrac. R squares shows how well the NMDS captures variation in the data, hence the higher the better. 


![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/wunifrac.nmds.logtransformed.jpeg)
> Figure 13. NMDS plot for weighted UniFrac distance. Each point represents one sample with different colors correspondent to treatments and the segment/sample type as the shape for the points. 

### Network-based beta diversity

You can also create a #graph or #network for understaning the similairty/dissimilarity of bacterial distributions between treatments based on distance metrics, e.g. Bray-Curtis dissimilairty matrix. The samples with the similar distributions form solid edges, and with close relations they form mixed edges and with different distributions they form no edges with one another. Then you can perform a permutational test to verify the graph's validity.

```R
library("phyloseqGraphTest")
library("igraph")
library("ggnetwork")

# Creating the based on Bray-Curtis dissimaliry matrix

gt = graph_perm_test(pst.qPCR, "treatment",
                    distance = "bray", type = "mst")
gt$pval

# Graph 

plotNet = plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text("Network based analysis of treatments based on\n Bray-Curtis, P = 0.002", size = 9)) + 
labs(col = "Treatment") + geom_nodes(aes(color = sampletype), inherit.aes = T, size = 5) + 
scale_color_manual(values = c("deeppink1", "darkorange", "deepskyblue",   "springgreen4"))

#Permutaion test histogram
plotPerm1 = plot_permutations(gt) + theme_bw() + geom_text(label = glue("{gt$pval}"), aes(x = 55, y = 10), color = "red")

#Patching two plots together
net.treat = grid.arrange(ncol = 2, plotNet, plotPerm1) 

ggsave(filename =  "./Network/treatment.net.bray_pval0.002.jpeg", net.treat, dpi = 300, device = "jpeg", width = 12, height = 10)

```

![graph_bray](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/treatment.net.bray_pval0.002.jpeg)
> Figure 14. Graph-based plot based on Bray-Curtis dissimilarity with minimum-spanning tree (MST) test statistic.


#

## 6. Statistical analysis on beta diversity: a Distance-based Redundancy Analysis (dbRDA)

In the prevoius section, we showed how to visualized differences between samples belonging to different treatments in diversity of microbiota.
First we create a distance matrix by either `vdist` from `vegan` package or by `distance` from `phyloseq` pakcage. 
```R
##WUniFrac log
wunifrac.dist.qpcr.log = phyloseq::distance(pst.qPCR.log, method = "wunifrac")#weighted unifrac distance on log-transformed qpcr data dataset
```

Then we test to see if the disperssion of the variance around the treatment centroids are homogenous (variance homogeniety test) by `betadisper` function from `vegan` package. If the p-value of test statistic is significant, it means that there is a significant difference in variance for any of the tested levels and the results from your `dbrda` model will not be reliable.

```R
set.seed(10)

# Making a permutational iteration to be passed on later
h <- with(data = data.frame(sample_data(pst.qPCR.log)), how(blocks = litter, nperm = 9999))
wunifrac.disp <- vegan::betadisper(wunifrac.dist.qpcr.log, group = sample_data(pst.qPCR.log)$treatment, 
                                 type = "centroid")
# Anova test
anova(wunifrac.disp, permutations = h)
#permutest(wunifrac.disp, permutation =h)
```
![beta.disper](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/beta.disper.PNG)
> Figure 15. Test statistics for variance homogeniety test. Betadisper test statistics are not significant, meaning that we can accept the null hypothesis that our groups have the same dispersions of variance. This means we can be confident that our adonis/dbrda result is maninly due to biological differences due to the treatment effect rather than differences in variance dispersions.


You can plot and save the PCoA of beta dispersion test:
```R
#to remove all open graphic devices if jpeg function doesn't save the picture
for(i in dev.list()[1]:dev.list()[length(dev.list())]){
   dev.off()
    }

jpeg( "./Beta diversity/dispersion of variance_wunifrac.jpeg", quality = 100)

plot(wunifrac.disp, col = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4"), bty = "n",
  las = 1, main = "Dispersion of variance around the centroids, WUniFrac", sub=NULL,
  xlab = sprintf("PCo1 [%s%%]", round(eig.vals/sum(eig.vals)*100,0)[1]),
 ylab = sprintf("PCo2 [%s%%]", round(eig.vals/sum(eig.vals)*100,1)[2])); text("P = 0.09", 
                                                               x = 0.15, y = -0.12, cex = 1.5)
```

![beta.disper.pcoa](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/dispersion%20of%20variance_wunifrac.jpeg)
> Figure 16. PCoA of variance homogeniety test around the centroids for different treatments. In significant test statistic shows indicates a homogenous variance dispersion for all centroids.

#
### Distance-based Redundancy Analysis (dbRDA)
Distance-based Redundancy Analysis (dbRDA) is a Redundancy Analysis on the eigenvalue-scaled result of Principal Coordinates Analysis. This is a method for carrying out constrained ordinations on data using non-Euclidean distance measures, with the assumption of liniear relationship between response and environmental variables ([this is not the case especially in ecological data](https://sites.ualberta.ca/~ahamann/teaching/graphics/LabRDA.pdf)). The usual methods for constrained ordinations (CCA, RDA) use Euclidean distance, but this does not work for all data, such as community count data, e.g. wunweighted UniFrac distance. 

```R
set.seed(1990)
h <- with(data = data.frame(sample_data(pst.qPCR.log)), how(blocks = litter, nperm = 9999))

db.rda.wunifrac = vegan::dbrda(wunifrac.dist.qpcr.log ~  gb * dss + sample_type  + Condition(litter), 
data = sample_data(pst.qPCR.log)%>%data.frame)#dbRDA does not do any permutation test

#but it is anova that does the permutations in order to generate the psudo-F statistics. 

#inertia is the squared wunifrac distance here

db.rda.wunifrac 
```
![dbrda.model](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/dbrda.model.PNG)
> Figure 17. The dbrda model summary. Constrained is related to  the variance defiend by our treatments and the model component, and unconstrained is unexplained variance and is shows it for different MDS axis. Conditional term, here is our litter, is what we limitted the permutation within. Proportion column shows the proportion of variance explained by different components of our model.


```R
anova(db.rda.wunifrac, permutation =h) #first do the omnibus test to give you an overal test on the model
anova(db.rda.wunifrac, permutation =h, by = "margin")
anova(db.rda.wunifrac, permutation =h, by = "term")
anova(db.rda.wunifrac, permutation =h, by = "axis")
```
![dbrda.omnibus](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/dbrda.omnibus.PNG)

> Figure 18. dbRDA test statistics and psudo-F for the whole model. F is the psudo-F, which is the ratio between the variance of the tested variables and the residual variance.

![dbrda.margin](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/dbrda.margin.PNG)

Figure 19. dbRDA test statistics and psudo-F for the marginal effects.

![dbrda.term](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/dbrda.term.PNG)

Figure 20. dbRDA test statistics and psudo-F for different independent variables passed on to the model.

![dbrda.axis](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/dbrda.axis.PNG)

Figure 21. dbRDA test statistics and psudo-F for different dbRSA axis.


Now you can extract the `site` and `centroid` scores from your `db.rda.wunifrac` artifact and make an ordination plot out of the model, which is the variation explained only by the treatments. 

```R
# Extreacting the goodies :)
score.site = vegan::scores(db.rda.wunifrac, display = "sites") %>% as.data.frame
score.centroid = vegan::scores(db.rda.wunifrac, display = "cn") %>% as.data.frame
rownames(score.centroid) <- levels(sample_data(pst.qPCR)$treatment)

eig.vals = db.rda.wunifrac$CCA$eig
inertia.total = db.rda.wunifrac$tot.chi #Total variation (inertia) explained. This number should be used as the denominator for measuring the amount of variance out of totoal variance wxplained by each dbrda axis.

#Ordination Plot
score.site %>% ggplot(aes(dbRDA1, dbRDA2, shape = sample_data(pst.qPCR)$sample_type, 
                          color = sample_data(pst.qPCR)$treatment)) +
geom_point(size =5 ) + geom_hline(yintercept = 0, lty = 2, alpha =0.5) + geom_vline(xintercept = 0, lty = 2,
                                                                                  alpha = 0.5)+
coord_fixed() + scale_color_manual(values = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4")) + 
theme_bw() + scale_y_continuous(na.value = c(-2, 3), n.breaks = 10) +
scale_x_continuous(na.value = c(-1, 1), n.breaks = 10) + labs(col="Treatment", shape = "Segment") + 
xlab(label = paste("dbRDA1 [", round(eig.vals[[1]]/sum(eig.vals)*100, 1), 
                    "% of fitted and", round(eig.vals[[1]]/inertia.total*100, 1), "% of total variation]")) + 
ylab(label = paste("dbRDA2 [", round(eig.vals[[2]]/sum(eig.vals)*100, 1), "% of fitted and", round(eig.vals[[2]]/inertia.total*100, 1), "of total variation]")) + 
theme(axis.title = element_text(size = 15), text = element_text(size = 13, face = "bold"),
    axis.text.x =element_text(size =10, face = "bold"),
axis.text.y =element_text(size =10, face = "bold")) + ggtitle(label = "dbRDA plot of WUF scores for site")  # + 
#stat_ellipse(data = score.site, aes(dbRDA1, dbRDA2), type = "euclid", linetype = 2, inherit.aes = F, geom = "polygon" )

ggsave("./Beta diversity/wunifrac.dbRDA.scores.site.jpeg", height = 8, width = 9, dpi =300)
```

![dbrda.plot](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/wunifrac.dbRDA.scores.site.jpeg)
> Figure 22. Ordination plot of dbRDA site or sample scores. Each point reprents one sample. On different axis you can see the variance explained by our constrained varialbes, i.e. treatments after fitting the model and the other variance explained by treatment out of total variance. 

You can also fit environmental factors and taxa on this ordination plot. In this instance, I will try to see which chemicals, e.g. SCFA had associations with any treatments on our dbRDA plot.

In order to do so I load chemical data first:
```R

# Loading chemical data 
chem.dat = read.table("./qPCR_metadat.tsv")
sample_data(pst.qPCR) <- chem.dat

#changing the variable classes to factor if applicable

for(i in seq_len(ncol(sample_data(pst.qPCR)))) {
 if(!is.numeric(sample_data(pst.qPCR)[[i]]) && !is.logical(sample_data(pst.qPCR)[[i]])) {
    sample_data(pst.qPCR)[[i]] = as.factor(sample_data(pst.qPCR)[[i]])  } else {
     sample_data(pst.qPCR)[[i]]
 } 
}

```

Then we can fit chemical data to our dbrda model and make the ordination data frames for ploting.

```R
pst.spec <- gloomer(ps = pst.qPCR.log, taxa_level = "Species", NArm = TRUE)

#Envfit for chemical data
env.chem = sample_data(pst.spec)[,15:24] %>% as.matrix
env.chem = apply(env.chem, 2, function(x) log(1 + x)) #log transformation of chemical data
pst.spec.chem <-pst.spec
env.chem = env.chem[complete.cases(env.chem),]        # Removing NA

otu_table(pst.spec.chem) <- otu_table(pst.spec.chem)[, colnames(otu_table(pst.spec.chem)) %in% rownames(env.chem) ]

wunifrac.dist.qpcr = phyloseq::distance(pst.spec.chem, method = "wunifrac") # Making a new distance

set.see
dbrda.wunifrac.chem = dbrda(wunifrac.dist.qpcr ~  gb * dss  + Condition(litter), data = sample_data(pst.spec.chem)%>%data.frame)
h <- with(data = data.frame(sample_data(pst.spec.chem)), how(blocks = litter, nperm = 9999))

#Fitting data with envfit
rda.env  = envfit(dbrda.wunifrac.chem, env =  env.chem, perm = h, choice = c(1, 2), na.rm = TRUE)
                 
envfit = data.frame(rda.env$vectors[1]$arrows, p.vals = rda.env$vectors[4] ) 
                 
sig.envfit = envfit[envfit$pvals  < 0.05, ]
                 
env.chem = env.chem[, colnames(env.chem) %in% rownames(sig.envfit)]
                 
x = sig.envfit[,1]*0.15
y = sig.envfit[,2]*0.2
x0 = rep(0, length(x))
y0 = rep(0, length(y))

#plotting species
Y= t(otu_table(pst.spec.chem))
n = nrow(Y)

#standardizing the first two eigenvectors
ev.stand  = scale(dbrda.wunifrac.chem$CA$u[,c(1,2)])

#Then the cov function calculates covariance between chemical concentration and the standardized eigenvectors.

n = nrow(env.chem)

S = cov(env.chem, ev.stand)

U = S %*% diag((dbrda.wunifrac.chem$CCA$eig[c(1,2)]/(n-1))^(-0.19))

arrow.df = data.frame(PCoA1 = U[,1], PCoA2 = U[,2], taxon = rownames(U), x0 = rep(0, nrow(U), y0 = rep(0, nrow(U))))
arrow.df$y0 = rep(0, nrow(U))
                 
arrow.df2 = data.frame(PCoA1 = U[,1], PCoA2 = U[,2], taxon = rownames(U), x0 = rep(0, nrow(U), y0 = rep(0, nrow(U))))
arrow.df$y0 = rep(0, nrow(U))

# Making the labeled dbrda plot

score.site %>% ggplot(aes(dbRDA1, dbRDA2))  +
geom_point(size =4, aes( color = sample_data(pst.qPCR.log)$treatment, shape = sample_data(pst.qPCR.log)$sample_type) ) + geom_hline(yintercept = 0, lty = 2, alpha =0.5) + geom_vline(xintercept = 0, lty = 2,
                                                                                  alpha = 0.5)+
coord_fixed() + scale_color_manual(values = c("deeppink1", "deepskyblue", "darkorange",  "springgreen4")) + 
theme_bw() + scale_y_continuous(na.value = c(-2, 3), n.breaks = 10) +
scale_x_continuous(na.value = c(-1, 1), n.breaks = 10) + labs(col="Treatment", shape = "Segment") + 
xlab(label = paste("dbRDA1 [", round(eig.vals[[1]]/sum(eig.vals)*100, 1), 
                    "% of fitted and", round(eig.vals[[1]]/inertia.total*100, 1), "% of total variation]")) + 
ylab(label = paste("dbRDA2 [", round(eig.vals[[2]]/sum(eig.vals)*100, 1), "% of fitted and", round(eig.vals[[2]]/inertia.total*100, 1), "of total variation]")) + 
theme(axis.title = element_text(size = 15), text = element_text(size = 13, face = "bold"),
    axis.text.x =element_text(size =10, face = "bold"),
axis.text.y =element_text(size =10, face = "bold")) + ggtitle(label = "dbRDA plot of WUF scores for site") + 
geom_segment(data = arrow.df, aes(x0, y0, xend = PCoA1, yend = PCoA2), color = "purple", size = 0.5, 
             arrow = arrow(type = "closed", angle = 15, length = unit(0.75,"line"))) + 
geom_text(data = arrow.df, aes(PCoA1, PCoA2, label = taxon), position = position_nudge(0.1), size = 4, color = "black") 
#
ggsave("./ordinations/wunifrac.dbRDA.scores.site.labeled.jpeg", height = 8, width = 10, dpi =300)
```
![dbrda.ord.label](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/wunifrac.dbRDA.scores.site.labeled.jpeg)
> Figure 23. Ordination plot of dbRDA site or sample scores and labeled with chemical data. The direction of arrows indicates association between biogenic amines and red meat consumption for instance.

## 7. Differential abundance analysis of taxa: DESeq2
 
Now that we found out the main effects of red meat consumption and DSS treatment and their interaction were significant on beta diversity, we can move on to identify those bacteria/ASVs that were differentially abundant depending on different treatments in a certain taxonomic level. First, we gloom the taxa to `phylum` and `species` levels. Second, make a model design formula. Third, you can estimate the geometric means for the estimation of size factors and then you can pas it on to the `geomMeans` argument of `estimateSizeFactors` function and fit it by `DESeq` function with `Wald` test in a `parametric` fit type.

```R
#glooming the taxa 
#phyl.qpcr =  gloomer(ps = pst.qPCR, taxa_level = "Phylum", NArm = TRUE)
#spec.qpcr = gloomer(pst.qPCR, taxa_level =  "Species", NArm = TRUE)

#Creating the experimental design
coldat = sample_data(phyl.qpcr)[,1:8] %>% data.frame
coldat$Pen <- factor(coldat$Pen, levels = c("pn1", "pn2", "pn3", "pn4"))
coldat$blok <- relevel(coldat$blok, 1)
#Here we can also do the factorial model.
design = model.matrix(~ gb * dss + blok + sample_type, data = coldat)

#since different litters have been used in different blocks, there is a rank defficiency for the model when it comes to the interaciton. 
#Therefore, we remove those all-zero columns from our desing and use the updated desing for fitting the model
zeros = apply(design, 2, function(x) all(x==0))
zeros = which(zeros)              
#design = design[,-zeros]

# Making a deseq artifact with the facrotial design             
phyl_dds <- phyloseq_to_deseq2(phyl.qpcr, design = design )

# Calculating geometric means prior to estimate size factors
gm.mean = function(x, na.rm= TRUE) {
    exp(sum(log(x[x>0]), na.rm=na.rm)/length(x))
}

geo.mean = apply(counts(phyl_dds), 1, gm.mean)

phyl_dds = estimateSizeFactors(phyl_dds, geoMeans = geo.mean)

phyl_dds<-DESeq(phyl_dds, test = "Wald", fitType = "parametric")

# You can do the same for the species

spec_dds<-phyloseq_to_deseq2(spec.qpcr, design  = design) #for phylum the deseq object is, phyl_dds
geo.mean = apply(counts(spec_dds), 1, gm.mean)
spec_dds = estimateSizeFactors(spec_dds, geoMeans = geo.mean)
spec_dds<-DESeq(spec_dds, test = "Wald", fitType = "parametric")


```

You can now see the results anames by typing `resultsNames(phyl_dds)`, and based on that you can extract the results from the model.

```R
results(phyl_dds, name = "treatmentGB") %>% data.frame %>% filter(padj <= 0.05)
results(phyl_dds, name = "treatmentDSS") %>% data.frame %>% filter(padj <= 0.05)
results(phyl_dds, name = "treatmentGBDSS") %>% data.frame %>% filter(padj <= 0.01)
```


It is always good to test if there are outliers in your result table.

```R
#Test for outliers
phyl.res<-results(phyl_dds, cooksCutoff = FALSE)#to get the results of gb and gbdss and blok vs. ct and block 1
alpha=0.05
sigtab=phyl.res[which(phyl.res$padj<=alpha),]
sigtab.phyl=cbind(as(sigtab, "data.frame"), as(tax_table(phyl.qpcr)[rownames(sigtab),], "matrix"))

spec.res<-results(spec_dds, cooksCutoff = FALSE)#to get the results of gb and gbdss and blok vs. ct and block 1
alpha=0.05
sigtab=spec.res[which(spec.res$padj<=alpha),]
sigtab.spec=cbind(as(sigtab, "data.frame"), as(tax_table(spec.qpcr)[rownames(sigtab),], "matrix"))

all(rowMeans(counts(phyl_dds, normalized = TRUE, replaced = TRUE))==phyl.res$baseMean) #if there were outliers, it should have returned FALSE

all(rowMeans(counts(spec_dds, normalized = TRUE, replaced = TRUE))==spec.res$baseMean) #if there were outliers, it should have returned FALSE
```

You can also use `contrast` to extract the responds of different covariates in the model:

```R
# Phylum tables

gb.ct = results(phyl_dds, lfcThreshold = 0, contrast = list("gbYes")) %>% data.frame %>% filter(padj <= 0.05)
dss.ct = results(phyl_dds, lfcThreshold = 0, contrast = list("dssYes")) %>% data.frame %>% filter(padj <= 0.05)
gbdss.ct = results(phyl_dds, lfcThreshold = 0, contrast = list("gbYes.dssYes")) %>% data.frame %>% filter(padj <= 0.05)

# Species tables

gb.ct.spec = results(spec_dds, lfcThreshold = 2, contrast = list("gbYes")) %>% data.frame %>% filter(padj <= 0.01)
dss.ct.spec = results(spec_dds, lfcThreshold = 2, contrast = list("dssYes")) %>% data.frame %>% filter(padj <= 0.01)
gbdss.ct.spec = results(spec_dds, lfcThreshold = 2, contrast = list("gbYes.dssYes")) %>% data.frame %>% filter(padj <= 0.01)

# Downloading the phylum table results

write.table(gb.ct, "./gb_vs_ct.deseq.tsv", sep = ";")
write.table(dss.ct, "./dss_vs_ct.deseq.tsv", sep = ";")
write.table(gbdss.ct, "./gbdss_vs_ct.deseq.tsv", sep = ";")

# Downloading the species table results

write.table(gb.ct.spec, "./gb_vs_ct.deseq.tsv", sep = ";")
write.table(dss.ct.spec, "./dss_vs_ct.deseq.tsv", sep = ";")
write.table(gbdss.ct.spec, "./gbdss_vs_ct.deseq.tsv", sep = ";")
```

If you want to see the number of differentially abundant taxa between two treatments or their magnitude of change, you can do this:

```R
 #comparison between treatments in species levels
dss.n <- dss.ct.spec %>% filter(abs(log2FoldChange) > 2, padj <=0.01 ) %>% rownames()
gb.n <- gb.ct.spec %>% filter(abs(log2FoldChange) > 2, padj <=0.01 ) %>% rownames()
gbdss.n <- gbdss.ct.spec %>% filter(abs(log2FoldChange) > 2, padj <=0.01 ) %>% rownames()

# DSS vs GB comparison
compare.dss.gb <- matrix(NA, nrow = length(dss.n), ncol = length(gb.n), dimnames = list(dss.n, gb.n))

for(i in 1:length(dss.n)){
    for(j in 1:length(gb.n)){
        compare.dss.gb[i, j] = ifelse(dss.n[i] == gb.n[j], 1, 0)
    }
}

#GBDSS vs GB comparison
compare.gbdss.gb <-  matrix(NA, nrow = length(gbdss.n), ncol = length(gb.n), dimnames = list(gbdss.n, gb.n))

for(i in 1:length(gbdss.n)){
    for(j in 1:length(gb.n)){
        compare.gbdss.gb[i, j] = ifelse(gbdss.n[i] == gb.n[j], 1, 0)
    }
}

compare.dss.gb = compare.dss.gb %>% melt %>% filter(value > 0) %>% select (Var1) %>% pull %>% as.vector
compare.gbdss.gb = compare.gbdss.gb %>% melt %>% filter(value > 0) %>% select (Var1) %>% pull %>% as.vector

#Number of in-common species
compare.dss.gb %>% length
compare.gbdss.gb %>% length

#Different species between gbdss and dss
compare.gbdss.gb[!compare.gbdss.gb %in% compare.dss.gb]
```

Now it is time to make some figures out of our results. You can do Volcano and waterfall plots or make heatmap for the `Log2FoldChange` and/or `adjusted p-vlalue`.

```R
#Volcano plot for species

spec.gb.ct = results(spec_dds, contrast = list("gbYes.dssYes")) %>%data.frame %>% filter(!padj%in%NA)
spec.gb.ct$Significant = ifelse(spec.gb.ct$padj < 0.01, "FDR < 0.05", "Not Sig")
spec.gb.ct = rownames_to_column(spec.gb.ct, "Species")

spec.gb.ct %>% ggplot(aes(x = log2FoldChange, -log10(pvalue))) + 
geom_point( aes(color = Significant), size = 2) + 
scale_color_manual(values = c("forestgreen", alpha("deeppink", 0.2))) + 
theme_bw(base_size = 12) + theme(legend.position= "bottom") + 
               geom_text_repel(data= top_n(spec.gb.ct[spec.gb.ct$log2FoldChange < -25 & -log10(spec.gb.ct$pvalue) > 9,], 
               -2, log2FoldChange), #- means the lowest p.value (hihger significant)
               aes(label = Species), 
               size = 3, 
               box.padding = unit(0.6, "lines"),
               point.padding = unit(0, "lines"), max.overlaps = 3)+
geom_text_repel(data= top_n(spec.gb.ct[spec.gb.ct$log2FoldChange > 25 & -log10(spec.gb.ct$pvalue) > 15,], 2, log2FoldChange),   #- means the lowest p.value (hihger significant)
               aes(label = Species),  
               size = 3, 
               box.padding = unit(0.6, "lines"),
               point.padding = unit(0, "lines"), max.overlaps = 3) +
ggtitle (label = "Volcano Plot of the top 12 most significant log2FoldChange\n genus for GBDSS vs. CT") #+ scale_y_continuous(limits = c(-0.5, 40))

ggsave("./Deseq_species/volc_deseq_gbdss_vs_ct.jpeg", device = "jpeg", dpi = 300)
```

![volcano](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/volc_deseq_gbdss_vs_ct.jpeg) 
> Figure 24. Volcano plot of differential abundance species. X axis represents Log2FoldChange and Y axis is -log10 of adjusted p-value.


```R
#plotting for the phylum alone
theme_set(theme_bw())
sigtabphyl = dss.ct %>% rownames_to_column("Phylum")

#filtering results below 2 logFC
#sigtabphyl = sigtabphyl[abs(sigtabphyl$log2FoldChange)>2,]

#a costumized color scheme
sigtabphyl = sigtabphyl %>% mutate(phyl.col = ifelse(log2FoldChange > 0, "increased", "decreased"))


#filtering results above 0.01 padjust
alpha = 0.05
sigtabphyl = sigtabphyl[sigtabphyl$padj <=alpha,]
# Phylum order
x = tapply(sigtabphyl$log2FoldChange, sigtabphyl$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabphyl$Phylum = factor(as.character(sigtabphyl$Phylum), levels=names(x))

ggplot(sigtabphyl, aes(y=Phylum, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
geom_col( aes(fill = phyl.col), width = 1, show.legend = F, color = "white") +
  theme(title = element_text(size = 15, color = "black", face = "bold"),text = element_text(size =15, face = "bold"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 12, face = "bold"),
       axis.text.y = element_text(size = 13, face = "italic" )) + scale_x_continuous(limits = c(-30, 25),n.breaks = 10) + 
           ggtitle("Log2FoldChange of Phylum, 
DSS")  +  geom_text (mapping = aes(x=-23, y = 0.65, label = "FDR < 0.05, |LFC| > 0"), 
 color = "red", size = 4) + scale_fill_manual(values = c(  "deepskyblue", "darkgray"))  +
           geom_text(aes( label = round(log2FoldChange),2))
            #geom_text(aes(x=log2FoldChange, y = Phylum, label = -log10(padj) %>% round(1)), nudge_y = 0.16, color = "black", size = 4) 
ggsave("./Deseq_phylum/difabund_dss.jpeg", device = "jpeg", dpi = 300)
```
![waterfall](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/difabund_dss.jpeg)
> Figure 25. Barplot of differentially abundance for phylum. 


```R
#plotting for the Species and Phylum

theme_set(theme_bw())

spec.taxa = tax_table(spec.qpcr)
 

sigtabspec = cbind(as(dss.ct.spec, "data.frame"), as(spec.taxa[rownames(spec.taxa) %in% rownames(dss.ct.spec),], "matrix"))




#a costumized color scheme
phylcol=c('coral4', "darkorange",'gold','antiquewhite4', 'cornflowerblue', 'plum4',
          'darkgoldenrod3','aquamarine4', 'cadetblue2', 'red', 'darkblue', 'Maroon', 'Gray',
        'steelblue2','darkmagenta', 'tomato1', 'cyan4', 'darkgreen')

colindex = data.frame(color = phylcol[1:length(unique(tax_table(spec.qpcr)[,2]))], phylum = sort(unique(tax_table(spec.qpcr)[,2])))


phyla = unique(data.frame(tax_table(spec.qpcr)[,2])) %>% pull
colors = c()
for(i in phyla){
   colors[i] = colindex[colindex$Phylum == i,1]
}

#filtering out the taxa below 2 LFC
sigtabspec = sigtabspec[abs(sigtabspec$log2FoldChange)>2,]

#filtering results above 0.01 padjust
alpha = 0.01
sigtabspec = sigtabspec[sigtabspec$padj <=alpha,]
# Phylum order
x = tapply(sigtabspec$log2FoldChange, sigtabspec$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabspec$Phylum = factor(as.character(sigtabspec$Phylum), levels=names(x))


           
#Species reorder
x = tapply(sigtabspec$log2FoldChange, sigtabspec$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabspec$Species = factor(as.character(sigtabspec$Species), levels=names(x))
           
ggplot(sigtabspec, aes(y=Species, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
  geom_point(size=3)+
  theme(legend.text = element_text(face = "bold"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 12, face = "bold"),
       axis.text.y = element_text(size = 9, face = "italic"), title = element_text(face = "bold") ) +
           ggtitle("Log2FoldChange of Species 
for DSS") + geom_text (mapping = aes(x=-25,
       y =0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") + scale_x_continuous(limits = c(-31, 31),n.breaks = 10) + 
           scale_y_discrete(expand = c(0.00005,0.8)) + scale_color_manual(values = colors) + ylab("Species")
           #+ geom_text(aes( label = -log10(padj) %>% round(1)),inherit.aes = TRUE, nudge_y = 0.4, color = "black")
ggsave("./Deseq_species/difabund_spec_dss.jpeg", device = "jpeg", dpi  = 300, height = 15, width = 12)


```
![diff.abund.dss](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/difabund_spec_dss.jpeg)
> Figure 26. Waterfal plot of differentally abundant species for groups treated with DSS.
