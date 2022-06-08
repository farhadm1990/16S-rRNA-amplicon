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

Loading chemical data

```R

# I have added FaithPD to the dataset 
chem.dat = read.table("./qPCR_metadat.tsv")
sample_data(pst.qPCR) <- chem.dat



for(i in seq_len(ncol(sample_data(pst.qPCR)))) {
 if(!is.numeric(sample_data(pst.qPCR)[[i]]) && !is.logical(sample_data(pst.qPCR)[[i]])) {
    sample_data(pst.qPCR)[[i]] = as.factor(sample_data(pst.qPCR)[[i]])  } else {
     sample_data(pst.qPCR)[[i]]
 } 
}

```
