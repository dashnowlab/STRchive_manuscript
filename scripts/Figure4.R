
#so, lets count motifs from both alleles
total_with_row_number <- total %>%
  mutate(row_number = row_number())

STR_table_motif <- STR_table[c("gene", "benign_motif_reference_orientation",
                               "unknown_motif_reference_orientation")]


#normalize motifs
STR_table_motif$repeatunit_ref_normalized <- sapply(STR_table$reference_motif_reference_orientation, function(x) normalise_str(as.character(x)))
STR_table_motif$repeatunit_path_normalized <- sapply(STR_table$pathogenic_motif_reference_orientation, function(x) normalise_str(as.character(x)))
STR_table_motif$benign_motif_reference_orientation <- sapply(STR_table$benign_motif_reference_orientation, function(x) normalise_str(as.character(x)))
STR_table_motif$unknown_motif_reference_orientation <- sapply(STR_table$unknown_motif_reference_orientation, function(x) normalise_str(as.character(x)))

# account for CNG ambiguity
STR_table_motif <- STR_table_motif %>%
  mutate(repeatunit_path_normalized = ifelse(repeatunit_path_normalized == 'CNG', 'CAG,CCG,CGG,CTG', repeatunit_path_normalized))
STR_table_motif <- STR_table_motif %>%
  mutate(repeatunit_ref_normalized = ifelse(repeatunit_ref_normalized == 'CNG', 'CAG,CCG,CGG,CTG', repeatunit_ref_normalized))

# Convert motif_norm and motif_norm_small to separate rows, keeping only gene, row number, and the motifs
#you don't need row number but it's a good sanity check I think
motifs_by_allele <- total_with_row_number %>%
  select(gene, motif_norm, motif_norm_small, row_number) %>%
  pivot_longer(cols = c(motif_norm, motif_norm_small), names_to = "motif_type", values_to = "motif")

# Count unique motifs by grouping by motif and Id
motif_counts <- motifs_by_allele %>%
  group_by(gene, motif) %>%
  reframe(MotifCount = n()) %>%
  pivot_wider(names_from = motif, values_from = MotifCount, values_fill = 0)

# motif_counts_table <- motif_counts %>%
#   mutate(MotifCount = as.character(MotifCount),  # Convert MotifCount to character
#          MotifCount = paste(motif, MotifCount, sep = ": ")) %>%
#   group_by(Id) %>%
#   summarise(MotifCount = paste(MotifCount, collapse = ", ")) %>%
#   ungroup()

#get unique motif counts
unique_motif_counts <- motifs_by_allele %>%
  group_by(gene) %>%
  reframe(UniqueMotifCount = n_distinct(motif),
          motif = paste0("", paste(unique(motif), collapse = ","), ""))

# more than 1 is more interesting
unique_motif_counts <- unique_motif_counts %>%
  filter(UniqueMotifCount != 1)

unique_motif_counts_long <- unique_motif_counts %>%
  mutate(motif = strsplit(as.character(motif), ",")) %>%
  unnest(motif)

##merge it
unique_motif_counts_long <- merge(unique_motif_counts_long, STR_table_motif)


unique_motif_counts_long$classification <- apply(unique_motif_counts_long, 1, function(row) {
  m <- row['motif']
  if (any(m == unlist(strsplit(row['repeatunit_ref_normalized'], ',')))) {
    return("ref")
  } else if (any(m == unlist(strsplit(row['unknown_motif_reference_orientation'], ',')))) {
    return("unknown")
  } else if (any(m == unlist(strsplit(row['repeatunit_path_normalized'], ',')))) {
    return("path")
  } else if (any(m == unlist(strsplit(row['benign_motif_reference_orientation'], ',')))) {
    return("benign")
  } else {
    return("NEW_MOTIF?")
  }
})


# Reorder levels of classification factor
unique_motif_counts_long$classification <- factor(unique_motif_counts_long$classification,
                                                  levels = c("path", "unknown", "benign", "ref"))

# Plot with reordered levels
ggplot(unique_motif_counts_long, aes(x = gene, fill = classification)) +
  geom_bar() +
  labs(x = "Gene", y = "Unique Motif Count", fill = "Classification") +
  theme_minimal() +
  scale_fill_manual(values = c("benign" = "#00c7d5", "ref" = "darkblue", "unknown" = "lightgray", "path" = "#FD6853"),
                    labels = c("benign"= "Benign",
                               "ref" = "Reference",
                               "unknown" = "Unknown",
                               "path" = "Pathogenic"))
