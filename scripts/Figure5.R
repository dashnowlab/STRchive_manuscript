library(tidyverse)
library(dplyr)

#so, lets count motifs from both alleles
total_with_row_number <- total %>%
  mutate(row_number = row_number())

STR_table_motif <- STR_table[c("gene", "benign_motif_reference_orientation",
                               "unknown_motif_reference_orientation")]


# longreadmotifcounts <- tibble(
#   gene = c("YEATS2",
#            "RAPGEF2","C9ORF72","XYLT1",
#            "RFC1","FGF14","BEAN1",
#            "STARD7","SAMD12","ZFHX3"),
#   UniqueMotifCount = c(3L, 2L, 3L, 7L, 17L, 2L, 3L, 10L, 2L, 3L),
#   motif = c("TTTAT,TTATG,TGTTA",
#              "TTTAT,TTATG","GCCCCG,ACCGCA,CCCCGGGCCCGC",
#              "GCC,CGCGG,GGGCGC,GC,AGG,AGGCGGG,AGGG",
#              "GAAAG,GGGAA,AAAGG,GGAAAG,GGAA,AAAGGG,AGGAA,A,AAGAG,GAAAA,GAAGA,AAGGG,GGGAAGGAA,GGAAAA,GGAAA,AAGGA","AGC,GAA",
#              "ATAAC,ATAAA,A",
#              "AAATAAATAAAATA,AAAAT,AAATAAAATAACATA,AAAAG,ACAAA,AACAT,ATAAC,AACAC",
#              "ATAAC,AAAAT","CCGCCGCCACTGCCA,GCCACT,CCG")
# )
#
# longreadmotifcounts$motif <- sapply(longreadmotifcounts$motif, function(x) normalise_str(as.character(x)))


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
  dplyr::filter(UniqueMotifCount != 1)

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
  labs(x = "Gene", y = "Unique motif count", fill = "Classification") +
  theme_minimal() +
  scale_fill_manual(values = c("benign" = "#00c7d5", "ref" = "darkblue", "unknown" = "lightgray", "path" = "#FD6853"),
                    labels = c("benign"= "Benign",
                               "ref" = "Reference",
                               "unknown" = "Unknown",
                               "path" = "Pathogenic"))



### RFC1 analysis
# motifs_by_allele <- total_prev_with_row_number %>%
#   select(gene, motif_norm, motif_norm_small, row_number) %>%
#   pivot_longer(cols = c(motif_norm, motif_norm_small), names_to = "motif_type", values_to = "motif")
#
# RFC1_motif_counts <- subset(motifs_by_allele, motifs_by_allele$gene == 'RFC1')
#
# RFC1_motif_counts_motif_info <- merge(RFC1_motif_counts, STR_table_motif, by = "gene")
#
#
# RFC1_motif_counts_motif_info$classification <- apply(RFC1_motif_counts_motif_info, 1, function(row) {
#   m <- row['motif']
#   if (any(m == unlist(strsplit(row['repeatunit_ref_normalized'], ',')))) {
#     return("ref")
#   } else if (any(m == unlist(strsplit(row['unknown_motif_reference_orientation'], ',')))) {
#     return("unknown")
#   } else if (any(m == unlist(strsplit(row['repeatunit_path_normalized'], ',')))) {
#     return("path")
#   } else if (any(m == unlist(strsplit(row['benign_motif_reference_orientation'], ',')))) {
#     return("benign")
#   } else {
#     return("NEW_MOTIF?")
#   }
# })
#
#
#
# # Reorder levels of classification factor
# RFC1_motif_counts_motif_info$classification <- factor(RFC1_motif_counts_motif_info$classification,
#                                                       levels = c("path", "unknown", "benign", "ref"))
#
# ggplot(RFC1_motif_counts_motif_info, aes(x=motif, fill = classification)) +
#   geom_histogram(stat = "count") + scale_y_log10() +
#   labs(title = "Motif Counts for RFC1", x = "Motif", y = "Frequency") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust size as needed
#         axis.text.y = element_text(size = 14),  # Adjust size as needed
#         axis.title.y = element_text(size = 16)) +
#   scale_fill_manual(values = c("benign" = "#00c7d5", "ref" = "darkblue", "unknown" = "lightgray", "path" = "#FD6853"),
#                     labels = c("benign"= "Benign",
#                                "ref" = "Reference",
#                                "unknown" = "Unknown",
#                                "path" = "Pathogenic"))


# #### lr for long reads
# lr_unique_motif_counts_long <- longreadmotifcounts %>%
#   mutate(motif = strsplit(as.character(motif), ",")) %>%
#   unnest(motif)
#
# ##merge it
# lr_unique_motif_counts_long <- merge(lr_unique_motif_counts_long, STR_table_motif)
#
#
# lr_unique_motif_counts_long$classification <- apply(lr_unique_motif_counts_long, 1, function(row) {
#   m <- row['motif']
#   if (any(m == unlist(strsplit(row['repeatunit_ref_normalized'], ',')))) {
#     return("ref")
#   } else if (any(m == unlist(strsplit(row['unknown_motif_reference_orientation'], ',')))) {
#     return("unknown")
#   } else if (any(m == unlist(strsplit(row['repeatunit_path_normalized'], ',')))) {
#     return("path")
#   } else if (any(m == unlist(strsplit(row['benign_motif_reference_orientation'], ',')))) {
#     return("benign")
#   } else {
#     return("NEW_MOTIF?")
#   }
# })
#
#
# # Reorder levels of classification factor
# lr_unique_motif_counts_long$classification <- factor(lr_unique_motif_counts_long$classification,
#                                                   levels = c("path", "unknown", "benign", "ref"))
#
# # Plot with reordered levels
# ggplot(lr_unique_motif_counts_long, aes(x = gene, fill = classification)) +
#   geom_bar() +
#   labs(x = "Gene", y = "LR unique motif count", fill = "Classification") +
#   theme_minimal() +
#   scale_fill_manual(values = c("benign" = "#00c7d5", "ref" = "darkblue", "unknown" = "lightgray", "path" = "#FD6853"),
#                     labels = c("benign"= "Benign",
#                                "ref" = "Reference",
#                                "unknown" = "Unknown",
#                                "path" = "Pathogenic"))
#
#
#
# # Get common genes
# common_genes <- intersect(longreadmotifcounts$gene, unique_motif_counts$gene)
#
# lr_unique_motif_counts_filtered <- lr_unique_motif_counts_long %>%
#   filter(gene %in% common_genes)
#
# unique_motif_counts_filtered <- unique_motif_counts_long %>%
#   filter(gene %in% common_genes)
#
# # Add source column to lr_unique_motif_counts_filtered and unique_motif_counts_filtered
# lr_unique_motif_counts_filtered$source <- "lr_unique_motif_counts"
# unique_motif_counts_filtered$source <- "unique_motif_counts"
#
# # Add suffixes to the gene column in each data frame
# lr_unique_motif_counts_filtered$gene <- paste(lr_unique_motif_counts_filtered$gene, "_lr", sep = "")
# unique_motif_counts_filtered$gene <- paste(unique_motif_counts_filtered$gene, "_unique", sep = "")
#
# # Merge the data frames with suffixes
# merged_df <- rbind(lr_unique_motif_counts_filtered, unique_motif_counts_filtered)
#
# ggplot(merged_df, aes(x = gene, fill = classification)) +
#   geom_bar() +
#   geom_text(data = merged_df[grepl("_lr$", merged_df$gene), ],
#             aes(x = gene,
#                 y = (max(UniqueMotifCount))*1.1,
#                 label = "LR"),
#             vjust = -0.5,
#             size = 5,
#             color = "black") +
#   geom_text(data = merged_df[grepl("_unique$", merged_df$gene), ],
#             aes(x = gene,
#                 y = (max(UniqueMotifCount))*0.9,
#                 label = "gnomAD"),
#             vjust = -0.5,
#             size = 4,
#             color = "black") +
#   labs(fill = "Classification") +
#   scale_fill_manual(values = c("benign" = "#00c7d5", "ref" = "darkblue", "unknown" = "lightgray", "path" = "#FD6853"),
#                     labels = c("benign"= "Benign",
#                                "ref" = "Reference",
#                                "unknown" = "Unknown",
#                                "path" = "Pathogenic")) +
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#

