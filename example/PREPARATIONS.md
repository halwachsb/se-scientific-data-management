# gnomAD allele frequency example
Start with creating an account on [Posit Cloud](https://posit.cloud/). And create a new workspace for your new R project. Alternatively you can open your local R installation. 

**Note:** You should append all the commands consecutively in one single R script; the parts have been split only for presentation purposes and to prevent getting lost in a disorganized page on GitHub

**Note:** Any line that begins with a hash symbol "#"  is a comment and will not be executed. Read the comments. They describe what is happening and are helpful.

We recommend you add a header to the top of every coding script. This header should include basic details, similar to the example below. Note that you are free to choose the style and information included in your header.
```
## ---------------------------
##
## Script name: gnomadDataProcessing.R
##
## Purpose of script: worked-through example for the data management
##                    in R PhD seminar
##
## Author: Bettina Halwachs-Wenzl
##
## Date Created: 2025-07-10
##
## Copyright (c) Bettina Halwachs-Wenzl, 2025
## Email: bettina.halwachs-wenzl@uni-graz.at
##
## ---------------------------
##
## Notes: Minimal example header
##        feel free to adapt
##   
##
## ---------------------------
```


To make use of Bioconductor packages you need the BiocManger. Think of it as an AppStore. It has to be installed only once. Once it is available you can us it to install any Bioconductor package you want to use.

```
#Install Bioconductor packages if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
```

As we are prepared for our analysis, we know which packages have to be installed. This can be done either one by one, or as shown here with a comma separated list of packages at once.

Example of installing every needed packes one be one:
```
BiocManager::install("VariantAnnotation")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("AnnotationDbi")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
```


Example using a comma separated vector - c(element1,element2,element3,element4,...)
```
BiocManager::install(c("VariantAnnotation", "GenomicRanges","GenomicFeatures","AnnotationDbi","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db"))
```
After all the packages have been installed, they must be activated or loaded using the _library()_ command. Think of it as telling R which packages you intend to use.

```
# Load libraries - tell R with which function you wanna work
library(VariantAnnotation) # exploration and annotation of genetic variants
library(GenomicRanges) # representing genomic locations
library(GenomicFeatures) # For TxDb and mapping
library(AnnotationDbi)   # For the 'select' function
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # R interface to a pre-fabricated database 
                                           # of gene models for the human genome, 
                                           # based on the UCSC hg38 assembly and the knownGene track.
library(org.Hs.eg.db) # the organism database package

library(data.table) # Extension of 'data.frame'. Fast aggregation of large data
library(RCurl) # Functions to fetch URIs
library(ggplot2) # The Chuck Norris of data visualization
library(dplyr) # facilitate manipulation of data frames
library(tidyr) # For data cleaning
```
Congratulations you have now successfully installed and loaded all necessary packages which we need in our example.

# Continue with [Retrieve Data from gnomAD](GETVCF.md).


