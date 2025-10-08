# Visualization of gnomAD allele frequencis for ODC1

A table is nice but a graphical visualization is better. Before visualization the generated frequency data has to be processed/structured so it can be used with one of the most powerfull packages of R, _ggplot2_. 

Let's start with structuring ...
```
# Extract and convert INFO data to a tibble
info_data <- info(vcf_odc1)
info_df <- as_tibble(info_data)
```
Get the needed information to summarize it at the end in one object

```
# Extract fixed columns (REF, ALT) 
# here REF and ALT correspond to the REFerence allele and the ALTernate allele
fixed_data <- fixed(vcf_odc1)

# Convert DNAString/DNAStringSet columns (REF/ALT) to character vectors
fixed_df <- tibble(
  REF = as.character(fixed_data$REF),
  ALT = lapply(fixed_data$ALT, as.character) 
)
```

Extract variant coordinates and rs IDs (as before)
```
variant_ranges <- rowRanges(vcf_odc1)
```

Create a clean tibble for the genomic coordinates 
```
variant_df <- tibble(
  seqnames = as.character(seqnames(variant_ranges)),
  start = start(variant_ranges),
  end = end(variant_ranges),
  # Extract the problematic rs IDs into a new column
  rsID = names(variant_ranges) 
)
```
Continue with data styling ...

```
# Combine ALL data sources (Coordinates + Fixed + INFO)
all_data <- bind_cols(variant_df, fixed_df, info_df) %>%
  
# Unnest both ALT (now a list of characters) and AF to handle multi-allelic sites
unnest(cols = c(ALT, AF)) %>% 
  
mutate(
    # REF and ALT are now simple character columns, so the paste0 works!
    Variant_ID = paste0(seqnames, ":", start, ":", REF, ">", ALT)
)
``` 

Congratulations, almost done! The data looks already quite pretty. Now let's select for the top 30 alleles

```
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
```
Finally, generate an allele frequency bar graph using _ggplot2_

```
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
````

Optionally you can save the created plot to disk. Helpful when you create more than one and whant to look at them afterwards.

```
# 2. Save the plot
#ggsave("most_frequent_alleles_odc1.png", plot = allele_plot, width = 8, height = 6)
```

## Congratulations, you did it!

