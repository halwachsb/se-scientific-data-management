Since the installation of the necessary packages will take some time, we will start this process earlier during the seminar to be ready for the afternoon practical session.

```
#Install Bioconductor packages if you haven't already
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# Install all the needed Bioconductor packages
BiocManager::install(c("VariantAnnotation", "GenomicRanges",
                       "GenomicFeatures","AnnotationDbi","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db"))

# install all the needed CRAN packages
install.packages(c("data.table","ggplot2","dplyr","rmarkdown","RCcurl","tidyr"))

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

#### Be patient. It will take a while. Good luck!
