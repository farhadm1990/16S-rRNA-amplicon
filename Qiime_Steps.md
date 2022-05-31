Before we continue on this jurney, I highly recommend you to check out the doc vignette of the [qiime tutorials](https://docs.qiime2.org/2022.2/tutorials/) related to the [installation](https://docs.qiime2.org/2022.2/install/native/) and case studies, such as [Parkinson's Mouse Tutorial](https://docs.qiime2.org/2022.2/tutorials/pd-mice/) for instance. Therefore, I declare that all qiime-related codes are basically copied and costumized based on [qiime tutorials](https://docs.qiime2.org/2022.2/tutorials/) and at any applicable point, I tried to provided the reader with the relevant references. If you come across something without reference, please don't go mad at me and just let me know about it ;)

## 1. Importing raw data into Qiime2
Here we import the raw sequences into a qiime artifact. But first you need to download the amplicon data [here](https://www.dropbox1.com/scl/fo/o2j5uwiaynh6owsom9kf0/h?dl=0&rlkey=4qvl191j9zfx4332tfm4pul1k_notyet) and save it on your local drive. In order to bring the sequences from the local path (here they are saved on a folder in my cluster called 'amplicons', you must use this [bash command](https://github.com/farhadm1990/Microbiome_analysis/blob/main/scripts/import.sh). The input path must be defined by a [manifest file](https://github.com/farhadm1990/Microbiome_analysis/blob/main/manifestArranged.tsv), which includes the name of the sample and a path to each sample sequence for both forward (in one column) and reverse reads (in another column). Sequence data are paired end in the format of FASTA with quality score (Fastaq); therefore, in qiime2 the type will be "SampleData[PairedEndSequencesWithQuality]" and their imput format asigned as PairedEndFastqManifestPhred33V2.  
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

You can submit this bash script, which has a `.sh` format to the cluster by the following command:

```python
sbatch import.sh
```

This might take a while before you get the results. The output of this, is a `.qza` file that you have already specified it in the command, in this example `demuxed-dss.qza`. You can then create a visualized file from this artifact, with the following command:

```python
qiime demux summarize \
  --i-data ./demuxed-dss.qza \
  --o-visualization ./demuxed-dss.qzv

```
This `demuxed-dss.qzv`is a visualized format of our `demuxed-dss.qza`, which you can view it on [qiime2 viewer](https://view.qiime2.org/). Once you are there you can either drag-and-drop the [downloaded](https://github.com/farhadm1990/Microbiome_analysis/blob/main/artifacts/demuxed-dss.qzv?raw=true) artifact into the designated area or simpley copy the link to the artifact from [this repository](https://github.com/farhadm1990/Microbiome_analysis/blob/main/artifacts/demuxed-dss.qzv) and paste it in the box *file from the web*. Once there, you must come across the following picture:

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Demux.PNG)
> **Figure 1**. Demultiplexed pairedEnd reads.

On this `overview` page you can see counts of demultiplexed sequences for the entire samples for both forward and reverse reads, with min, median, mean and max and total counts. 

In the `Interactive Quality Plot` at the top left, you can see the quality (Phred, Q) score for the sequence bases in forward and reverse reads. 

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/demux-interactive.PNG)
> Figure 2. Interacvive plot for demultiplexed pairedEnd reads. 

Understanding this plot is very important for your denoising step, i.e. you must decide a truncation length for each forward and reverse reads in a way to keep at least 50% of the reads equal or above Q = 30. You can see this changes by hovering over the interactive box plots. You can see that in forward reads for example, the quality of reads starts falling in sequence baces of above 260 nt and already in 220 nt for reverse reads.

## 2. Filtering, dereplication, sample inference, chimera identification, and merging of paired-end reads by DADA2 package in qiime2.
In this step we use [DADA2](https://www.nature.com/articles/nmeth.3869) package to denoise the reads and create the ASV tables and representative sequences (repseqs). 
Denoising is refered to sequence inference by resolving the sequencing erroros for each sequence as small as 1 nucleotide, which creates Amplicon Sequence Variants (ASVs). Denoising is superior over Operationals Taxonomic Units (OTU), which is clustering sequence reads into [OTUs](https://peerj.com/articles/5364/) at a defined cut-off rate for similarity (≥ 97), in a way that OTU clusters two different sequnces as the same OTU even if they are different in 3 nt. While, ASV can be as acurate as 1 nt difference. 
For denoising in qiime2, you can use the following code chunk in a bash file and submit it to your local cluster by `sbatch denoise.sh`:)

```python
#!/bin/bash
#SBATCH -p ghpc_v3 #name of the cluster node
#SBATCH -N 1
#SBATCH -n 10     #number of the processor cores, you will chose this more than 1, only if your qiime2 package/command can support it, e.g. it has e --p-n-threads 
#SBATCH --mem=64G #capacity of the total processors reserved for the job
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ~/data/dss/demuxed-dss.qza \         #the input file for denoising is our demultiplexed pairedEnd reads that was generated in the previous step.
--p-trim-left-f 17 \                                        #length of my forward primer (17 nt) 
--p-trim-left-r 21 \                                        #length of my reverse primer (21 nt)
--p-trunc-len-f 260 \                                       #truncation length for forward reads. I.e. we truncate reads over 260 base since their quality started dropping from this point onwards. 
--p-trunc-len-r 220 \                                       #truncation length for reverse reads. I.e. we truncate reads over 220 base since their quality started dropping from this point onwards. 
--o-table ~/data/dss/tableNoFilt.qza \                      #the ASV count table
--o-representative-sequences ~/data/dss/repseqsNoFilt.qza \ #the representative sequences for in each read
--o-denoising-stats ~/data/dss/denoising-statsNoFilt.qza \  #the status of the denoising process in a full catalogue
--p-n-threads 10                                            #number of the processors assinged for this task.

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID

```
Now you can visualize the ASV table to create a `.qzv` by the following command:
```python
qiime feature-table summarize --i-table ./tableNoFilt.qza \ #The ASV table as the input
--m-sample-metadata-file ./metadataArranged.tsv \           #The metadata as the input
--o-visualization ./tableNoFilt.qzv                         #The qzv format of our ASV table
```
As you can see we have a [metadata](https://docs.qiime2.org/2022.2/tutorials/metadata/?highlight=metadata) field, which is required for our visualizations and downstream analysis. A metadata is basically a tab-delimited, usually `.tsv` or `.csv`, file which links each sample to its original parents, i.e. the SampleID, the design of the experiemnt, the type of treatment etc.
Bellow you can see our [metadata file](https://github.com/farhadm1990/Microbiome_analysis/blob/main/metadataArranged.tsv):

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Metadata.PNG)
> Figure 3. A small metadata including the sample IDs, treatment, block etc.

You can also do the visualization for the `repseqsNoFilt.qza` file as follows:

```python
qiime feature-table tabulate-seqs \
--i-data repseqsNoFilt.qza \
--o-visualization repseqsNoFilt.qzv
```
And visualization for the `denoising-statsNoFilt.qza` file:

```python
qiime metadata tabulate \
--m-input-file ./denoising-statsNoFilt.qza \
--o-visualization ./denoising-statsNoFilt.qzv
```
If you drag and drop the `tableNoFilt.qzv` file in [qiime2 view](https://view.qiime2.org/), you can see three main menues; `Overview`, `Interactive Sample Detail` and `Feature Detail`. If you click on `Interactive Sample Detail` you can see a slider to the left of the picture which could be changed, based which you can arbiterarily decide, to which depth of reading you can do your rarefaction. Nonetheless, I am not going to do any rarefaction or preprocessings in qiime, but rather I will continue to create other artifacts in qiime2 and further transfer them to R. For now, you can take a look at the `tableNoFilt.qzv` file in qiime viewer. 

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/ASV%20table%20in%20qiime.PNG)
> Figure 4. ASV table indicating number of reads per sample.

You can move the left right slider to see how many features you would keep in how many samples. If you want to keep doing the downstream analysis, you can use this indicator as a premise to decide which reading depth you choose for rarefaction. 

## 3. Training a primer-based region-specific classifier for taxonomic classification by Naïve-Bayes method (in Qiime2)
For taxonomic classifications, you need to have a classifier to which you blast your sequences against to find out which taxonomic groups each sequence belongs to. This is also called reference phylogeny, which is a cruitial step in identifying the marker genes (in this case 16S rRNA) taken from different environmental a in saco samples. In order to do so, there are different 16S rRNA databases, of which [Greengens](https://www.nature.com/articles/ismej2011139) and [SILVA](https://www.arb-silva.de/) are well-known databases for the full length of 16S rRNA genes. You can always download the pre-trained classifiers at the [Data Resources](https://docs.qiime2.org/2022.2/data-resources/) of qiime2 website. However, it is always safe to train your classfier based on your own primersets and based on a Naive-Bayesian method. In order to do so, I have followed [this toturial](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494#heading--sixth-header) by [Mike Robeson](https://forum.qiime2.org/u/SoilRotifer) and I have used SILVA 138 dataset for training my classifier with the follwoing command:

```python
#!/bin/bash
#SBATCH --mem-per-cpu 64G
#SBATCH -c 1
#SBATCH -t 24:00:00

source activate qiime2
TMPDIR=/scratch/$SLURM_JOB_ID 
qiime feature-classifier fit-classifier-naive-bayes \   # here you can use train your classifier based on Naive-Bayes method
--i-reference-reads ~/derepseqs-uniq-341f-805r.qza \    # Dereplicated sequences based on the primer set as input
--i-reference-taxonomy ~/dereptaxa-uniq-341f-805r.qza\  # Dereplicated taxonomic annotations based on the primer set as input
 --o-classifier ~/silva-classifier-primered4.qza        # The final classifier as the output
```
This `bash` task took around 12h on a cluster with 64G memory capacity. After you got `silva-classifier-primered4.qza` classifier file, you can use it for your taxonomic classifications as follows:

```python
#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 10  
#SBATCH --mem=64G
#SBATCH -t 2:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime feature-classifier classify-sklearn \                               # Sklearn package for classification
--i-reads ~/data/dss/repseqsNoFilt.qza \                                  # Representative sequences as the (input)
--i-classifier ~/data/dss/classifier/silva138-classifier-341f-805r.qza \  # Our costumized SILVA 138 classifier (input)
--o-classification ~/data/dss/Taxonomy/taxonomyNoFilt.qza --p-n-jobs 10   # The taxonomic annotation linked to the repseqs (output)

cd $SLURM_SUBMIT_DIR 
rm -rf /scratch/$USER/$SLURM_JOBID
```
And then you can visualize your taxa table with the following code:

```python
qiime metadata tabulate \
--m-input-file ~/data/dss/Taxonomy/taxonomyNoFilt.qza \
--o-visualization ~/data/dss/Taxonomy/taxonomyNoFilt.qzv
```
It might look like the bellow figure. You can see a feature ID correspondent to each sequences, the taxonomic order, and the confidence interval for this classification. 

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Taxa%20table.PNG)
> Figure 5. A taxonomy table should containe the sequence/feature ID, taxon and maybe confidence interval. 

## 4. Creating a phylogenetic tree using SATE-enabled phyhlogenetic placement (SEPP) method for diversity analysis
There are different methods for generating the phylogenetic tree, e.g. *de novo* clustering and sequence insertion. According to the [authors](https://www.worldscientific.com/doi/abs/10.1142/9789814366496_0024) on this method, Phylogenetic placement is refered to inserting short molecular sequences (called query sequences) into an existing phylogenetic tree and aligning it on full-length sequences for the same gene. This can provide information beyond pure “species identification” but also about the evolutionary relationships between these query sequences and to known species. Phylogenetic placement operates in two steps; first, an alignment is estimated for each query sequence to the alignment of the full-length sequences, and then that alignment is used to find the optimal location in the phylogenetic tree for the query sequence ([Tandy Warnow, 2015](https://link.springer.com/referenceworkentry/10.1007/978-1-4899-7478-5_711)). You can download the `SEPP` dataset [here](https://data.qiime2.org/2022.2/common/sepp-refs-gg-13-8.qza).  
In this mdoule, we can create our phylogenetic tree using ASV table and repseqs as input using `fragment-insertion sepp` package in Qiime2 in the following bash script:

```python
#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=120G
#SBATCH -t 24:00:00

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPTDIR
mkdir -p $TMPDIR

source activate qiime2.8
qiime fragment-insertion sepp \                             # The package for generating the phylogenetic tree
--i-representative-sequences ~/data/dss/repseqsNoFilt.qza \ # Our representative sequences (input)
--i-reference-database ~/data/dss/sepp-ref-gg-13-8.qza \    # The SEPP sequence dataset (input)
--o-tree ~/data/dss/treeNoFilt.qza \                        # The unrooted tree (output)
--o-placements ~/data/dss/tree-placementsNoFilt.qza \
--p-threads 10                                              # Number of processor cores

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
```

You can visualize your tree on [iTol](https://itol.embl.de/upload.cgi) which might look like this:

![alt text](https://github.com/farhadm1990/Microbiome_analysis/blob/main/Pix/Tree.PNG)
> Figure 6. The unrooted phylogenetic tree generated by SEPP method displayed on iTol.



