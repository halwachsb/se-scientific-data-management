### [< Back to previous topic](/example/PREPARATIONS.md)

# Get VCF data about ornithine decarboxylase 1 - ODC1

If you go to gnomAD and [search for ODC1](https://gnomad.broadinstitute.org/gene/ENSG00000115758?dataset=gnomad_r4) you can have a look at the available data in the gnomAD browser online and extract information manually. BUT image you are interesed in more than ONE gene.

All the data can be downloaded/accessed via https://gnomad.broadinstitute.org/downloads#v4-variants.
This example shows how this can be done programmatically for one gene of interest, but of course, it can also be scaled up to hundreds.

Since the provided files are several gigabytes in size, we will use Bioconductor packages, which allow us to stream the data for processing without loading it all into memory.

# 1. Define the gnomAD VCF files 

```
# load the location/path/link of the vcf information to a variable
gnomad_exome_chr2_vcf <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr2.vcf.bgz"
```

To stream the VCF information we use the tbi Tabix index file.
Tabix index file for a compressed Variant Call Format (VCF) file, which allows for rapid, random access to specific regions of a large, compressed VCF file without needing to decompress the entire file

```
# create a variable which holds the name/path of the tbi file
gnomad_exome_chr2_tbi <- paste0(gnomad_exome_chr2_vcf,".tbi")

```

# 2. Define the region for the ODC1 gene on GRCh38
This region is used to query the gnomAD v4 VCF files to stream for variants. Note: This region changes with every gene of interest.

**Question** Who knows the region for our gene of interest by hard? ;)
No worries. We can look it up somewhere OR we can query it from the annotated genome db here in R.

We have installed the annotated reference genome hg38 before. Now we "load" it into a variable _txdb_ so we can work with the information.

```
#load the annotated reference genome into a variable
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
```

Define the gene of interest, ODC1

```
#assign ODC1 as character to variable gene_symbol
gene_symbol <- "ODC1"
```
As the gene symbol sometimes is ambigues we will use geneID/ entrezID. Again: You can look it up OR you select it from the Bioconductor ressource AnnotationDbi.

## 2.A translate gene symbol into geneID/entrezID
```
# get geneID for ODC1 using AnnotationDbi
gene_id_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_symbol,
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

# print the result of the mapping to the R console
print(gene_id_mapping)

```

The manually assesed geneID is 4953. Furtunately it matches the result of the AnnotationDbi query. Let's store it in a variable for later use.
```
# assign the value of column ENTREZID of the gene_id_mapping object into a new variable called, gene_id
gene_id<-gene_id_mapping$ENTREZID
```

Now we are still missing the Range of the gene, which is important for screening for variants. Again, we can look it up OR we select it again from the AnnotationDbi.
In addition we will narrow the range down to cds exons.

## 2.B Get the CDS regions for the ODC1 gene (by gene ID)
We will again make use of the information provided by a Bioconductur package [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
```
# 'cdsBy' returns a GRangesList, where each element is a list of CDS exons for one transcript
cds_by_gene <- GenomicFeatures::cdsBy(txdb, by = "gene")
```
Now extract the information for our gene of interest ODC1
```
# Extract the CDS for our specific gene ID
# The result is a GRangesList for all transcripts of ODC1
odc1_cds_list <- cds_by_gene[[gene_id]]
```
Combine all coding exons which are stored in the GRangesList object.
```
#single GRanges object containing the non-overlapping coding regions.
odc1_cds_reduced <- reduce(odc1_cds_list)
```
Hang in there. We are almost there
Process the results to get the overall genomic range for the Tabix query.
Find the chromosome, strand, and overall start/end of the combined coding sequence.

```
chrom_cds <- as.character(seqnames(odc1_cds_reduced))[1]
strand_cds <- as.character(strand(odc1_cds_reduced))[1]
```

Determine the overall Start and End points of the *entire* coding region
```
# Determine the overall Start and End points of the *entire* coding region by min and max of the combined coding regions
min_cds_start <- min(start(odc1_cds_reduced))
max_cds_end <- max(end(odc1_cds_reduced))
```


**Create the final GRanges object for the full CDS region**

```
target_range <- GRanges(
  seqnames = chrom_cds,
  ranges = IRanges(start = min_cds_start, end = max_cds_end),
  strand = strand_cds
)

# print the result with a message/comment to the R console
print("New target range (encompassing CDS region):")
print(target_range)
```

Congratulations! You have created an GRanges object for your gene of interest providing the information about the region which should be screend for variants.

### Continue now with analyzing the VCF data for you gene of interest [Analyze VCF stream](PROCESSVCF.md)
