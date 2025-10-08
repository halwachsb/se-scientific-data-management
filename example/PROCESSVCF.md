### [< Back to previous Topic](/example/GETVCF.md)

# 3. Check the defined genomic range

The following step is optional BUT ensures that the perviously defined URLs exist
```
# Check to ensure the VCF and its index are available
if(!url.exists(gnomad_exome_chr2_vcf) || !url.exists(gnomad_exome_chr2_tbi)) {
  stop("gnomAD VCF or TBI index file not found at the specified URL. Please check the current link.")
}
```

## 3.A. Create a TabixFile object from the predefined cloud based ressources
```
# tabix_file object contains now the gigabites of information of the VCF files
tabix_file <- TabixFile(gnomad_exome_chr2_vcf, index=gnomad_exome_chr2_tbi)
```

## 3.B.  Stream the VCF data for the ODC1 region or any other region you define
```
#Stream the VCF data for the ODC1 using the target_range defined before
# use cat() for displaying status updates, progress reports, or simple textual information 
cat("Streaming data for ODC1 from gnomAD exomes VCF...\n")

# stream the VCF information stored in the tabix_file
vcf_odc1 <- readVcf(tabix_file, genome = "GRCh38", param = target_range)
cat(paste("Successfully imported", nrow(vcf_odc1), "variants in the ODC1 region.\n"))
```


# 4. Summarize the result for our gene

```
#Extract the INFO field, which contains the allele frequency (AF) data
info_df <- as.data.frame(info(vcf_odc1))
````

# 5. Extract various parameters from the result

## 5.A. Extract the overall Allele Frequency (AF) and Allele Count (AC)

Note: The AF field can be a list if there are multi-allelic sites, so we need to process it.
We'll use the 'AF' (overall Allele Frequency) and 'AC' (Allele Count) fields.
```
#Create a data frame for results with the parameters we are interesed in
results_odc1 <- data.frame(
  CHROM = as.character(seqnames(vcf_odc1)),
  POS = start(vcf_odc1),
  REF = as.character(ref(vcf_odc1)),
  ALT = unlist(lapply(alt(vcf_odc1), paste, collapse=",")), # Handle multi-allelic ALTs
  ID = names(vcf_odc1)
)
```

## 5.B Extract Allele Frequencies (AF)
```
#The AF field often contains one frequency per ALT allele.
af_list <- info_df$AF
af_df <- as.data.frame(do.call(rbind, af_list))
colnames(af_df) <- paste0("AF_ALT", 1:ncol(af_df))
```

## 5.C Extract Allele Counts (AC) and Allele Number (AN)
```
ac_list <- info_df$AC
ac_df <- as.data.frame(do.call(rbind, ac_list))
colnames(ac_df) <- paste0("AC_ALT", 1:ncol(ac_df))

an_vector <- info_df$AN
```

## 5.D Combine results into one object
```
results_odc1 <- cbind(results_odc1, af_df, ac_df, AN=an_vector)
```

# 6. Make the output pretty

```
#Calculate AC / AN for sanity check for the first allele
results_odc1$Calculated_AF_ALT1 <- results_odc1$AC_ALT1 / results_odc1$AN
results_odc1$AF_Difference <- abs(results_odc1$AF_ALT1 - results_odc1$Calculated_AF_ALT1)
````

```
# Format for cleaner output
results_odc1$AF_ALT1 <- round(results_odc1$AF_ALT1, 5)
results_odc1$Calculated_AF_ALT1 <- round(results_odc1$Calculated_AF_ALT1, 5)
results_odc1$AF_Difference <- round(results_odc1$AF_Difference, 6)
```

```
# Select and display the most relevant columns
final_output <- results_odc1[, c("CHROM", "POS", "ID", "REF", "ALT", "AC_ALT1", "AN", "AF_ALT1")]
```

Finally, we can take a look at the top 30 results.
```
cat("\n")
cat("## Filtered Allele Frequencies for ODC1 (Top 30 variants) ##\n")
print(head(final_output[order(final_output$POS), ], 30))
```

While that is effective, a graphical visualization would offer a better presentation; therefore, let's proceed with creating one.

### Continue with [Graphical Visualization](VISUALIZATION.md)
