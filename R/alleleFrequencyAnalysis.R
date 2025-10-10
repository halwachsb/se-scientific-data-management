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

# IMPORTANT: Any line that begins with a hash symbol (#) is a comment and will not be executed.
# Read the comments. They describe what is happening and are helpful.


#Install Bioconductor packages if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install all the needed packages from Brioconductor
BiocManager::install(c("VariantAnnotation", "GenomicRanges",
                       "GenomicFeatures","AnnotationDbi","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db"))

# Install all the needed packages from the CRAN repository
install.packages(c("data.table","RCurl","ggplot2","dplyr","dplyr"))

# Once installed, the packages only need to be loaded for subsequent use.

# Load libraries - tell R with which function you wanna work
# Loading is the same for each package. No difference between Bioconductor or CRAN
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

# load the location/path/link of the vcf information to a variable
gnomad_exome_chr2_vcf <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr2.vcf.bgz"

# create a variable which holds the name/path of the tbi file
gnomad_exome_chr2_tbi <- paste0(gnomad_exome_chr2_vcf,".tbi")

# create a variable which holds the name/path of the tbi file
gnomad_exome_chr2_tbi <- paste0(gnomad_exome_chr2_vcf,".tbi")

#load the annotated reference genome into a variable
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#assign ODC1 as character to variable gene_symbol
gene_symbol <- "ODC1"

# get geneID for ODC1 using AnnotationDbi
gene_id_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_symbol,
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

# print the result of the mapping to the R console
print(gene_id_mapping)

# assign the value of column ENTREZID of the gene_id_mapping object into a new variable called, gene_id
gene_id<-gene_id_mapping$ENTREZID

# 'cdsBy' returns a GRangesList, where each element is a list of CDS exons for one transcript
cds_by_gene <- GenomicFeatures::cdsBy(txdb, by = "gene")

# Extract the CDS for our specific gene ID
# The result is a GRangesList for all transcripts of ODC1
odc1_cds_list <- cds_by_gene[[gene_id]]

#single GRanges object containing the non-overlapping coding regions.
odc1_cds_reduced <- reduce(odc1_cds_list)

chrom_cds <- as.character(seqnames(odc1_cds_reduced))[1]
strand_cds <- as.character(strand(odc1_cds_reduced))[1]

# Determine the overall Start and End points of the *entire* coding region by min and max of the combined coding regions
min_cds_start <- min(start(odc1_cds_reduced))
max_cds_end <- max(end(odc1_cds_reduced))

target_range <- GRanges(
  seqnames = chrom_cds,
  ranges = IRanges(start = min_cds_start, end = max_cds_end),
  strand = strand_cds
)

# print the result with a message/comment to the R console
print("New target range (encompassing CDS region):")
print(target_range)


#############################################################################################

# Check to ensure the VCF and its index are available
if(!url.exists(gnomad_exome_chr2_vcf) || !url.exists(gnomad_exome_chr2_tbi)) {
  stop("gnomAD VCF or TBI index file not found at the specified URL. Please check the current link.")
}

# tabix_file object contains now the gigabites of information of the VCF files
tabix_file <- TabixFile(gnomad_exome_chr2_vcf, index=gnomad_exome_chr2_tbi)

#Stream the VCF data for the ODC1 using the target_range defined before
# use cat() for displaying status updates, progress reports, or simple textual information 
cat("Streaming data for ODC1 from gnomAD exomes VCF...\n")

# stream the VCF information stored in the tabix_file
vcf_odc1 <- readVcf(tabix_file, genome = "GRCh38", param = target_range)
cat(paste("Successfully imported", nrow(vcf_odc1), "variants in the ODC1 region.\n"))

#Extract the INFO field, which contains the allele frequency (AF) data
info_df <- as.data.frame(info(vcf_odc1))

#Create a data frame for results with the parameters we are interesed in
results_odc1 <- data.frame(
  CHROM = as.character(seqnames(vcf_odc1)),
  POS = start(vcf_odc1),
  REF = as.character(ref(vcf_odc1)),
  ALT = unlist(lapply(alt(vcf_odc1), paste, collapse=",")), # Handle multi-allelic ALTs
  ID = names(vcf_odc1)
)

#The AF field often contains one frequency per ALT allele.
af_list <- info_df$AF
af_df <- as.data.frame(do.call(rbind, af_list))
colnames(af_df) <- paste0("AF_ALT", 1:ncol(af_df))

ac_list <- info_df$AC
ac_df <- as.data.frame(do.call(rbind, ac_list))
colnames(ac_df) <- paste0("AC_ALT", 1:ncol(ac_df))

an_vector <- info_df$AN

results_odc1 <- cbind(results_odc1, af_df, ac_df, AN=an_vector)

#Calculate AC / AN for sanity check for the first allele
results_odc1$Calculated_AF_ALT1 <- results_odc1$AC_ALT1 / results_odc1$AN
results_odc1$AF_Difference <- abs(results_odc1$AF_ALT1 - results_odc1$Calculated_AF_ALT1)


# Format for cleaner output
results_odc1$AF_ALT1 <- round(results_odc1$AF_ALT1, 5)
results_odc1$Calculated_AF_ALT1 <- round(results_odc1$Calculated_AF_ALT1, 5)
results_odc1$AF_Difference <- round(results_odc1$AF_Difference, 6)

# Select and display the most relevant columns
final_output <- results_odc1[, c("CHROM", "POS", "ID", "REF", "ALT", "AC_ALT1", "AN", "AF_ALT1")]

# Select and display the most relevant columns
final_output <- results_odc1[, c("CHROM", "POS", "ID", "REF", "ALT", "AC_ALT1", "AN", "AF_ALT1")]

cat("\n")
cat("## Filtered Allele Frequencies for ODC1 (Top 30 variants) ##\n")
print(head(final_output[order(final_output$POS), ], 30))




######################################################################################################

# Extract and convert INFO data to a tibble
info_data <- info(vcf_odc1)
info_df <- as_tibble(info_data)

# Extract fixed columns (REF, ALT) 
# here REF and ALT correspond to the REFerence allele and the ALTernate allele
fixed_data <- fixed(vcf_odc1)

# Convert DNAString/DNAStringSet columns (REF/ALT) to character vectors
fixed_df <- tibble(
  REF = as.character(fixed_data$REF),
  ALT = lapply(fixed_data$ALT, as.character) 
)

variant_ranges <- rowRanges(vcf_odc1)


variant_df <- tibble(
  seqnames = as.character(seqnames(variant_ranges)),
  start = start(variant_ranges),
  end = end(variant_ranges),
  # Extract the problematic rs IDs into a new column
  rsID = names(variant_ranges) 
)

# Combine ALL data sources (Coordinates + Fixed + INFO)
all_data <- bind_cols(variant_df, fixed_df, info_df) %>%
  
  # Unnest both ALT (now a list of characters) and AF to handle multi-allelic sites
  unnest(cols = c(ALT, AF)) %>% 
  
  mutate(
    # REF and ALT are now simple character columns, so the paste0 works!
    Variant_ID = paste0(seqnames, ":", start, ":", REF, ">", ALT)
  )

# Filter for the most frequent alternate alleles (AF > 0)
frequent_alleles <- all_data %>%
  filter(AF > 0) %>%
  # Select the top 30 most frequent alleles
  arrange(desc(AF)) %>%
  head(30) %>%
  # Create a cleaner label for the plot
  mutate(Label = paste0(Variant_ID, " (", rsID, ") - ", round(AF * 100, 1), "%"))


# Display the summary of the top alleles
cat("Successfully processed and filtered data:\n")
print(frequent_alleles)

# Create the Bar Plot
allele_plot <- ggplot(frequent_alleles, aes(x = reorder(Variant_ID, -AF), y = AF)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  geom_text(aes(label = paste0(round(AF * 100, 1), "%")), 
            vjust = -0.5, size = 3) +
  labs(
    title = paste("Top 30 Most Frequent Alleles in 4953 (ODC1)"),
    x = "Variant ID (CHROM:POS:REF>ALT)",
    y = "Allele Frequency (AF)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5)
  )

# the plot is stored in the _allele_plot_ variable. you either have to print or save it.
allele_plot

#  Optional: Save the plot. This is very helpful when running multiple analysis at once
# ggsave("most_frequent_alleles_ODC1.png", plot = allele_plot, width = 8, height = 6)
