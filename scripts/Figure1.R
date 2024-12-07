library(RColorBrewer)
library(waffle)

### Figure 1A
blue_palette <- colorRampPalette(c("lightgray", "darkgray", "#ADD8E6", "#4169E1", "#00008B", "black"))(6)


# Count the number of characters in pathogenic_motif_reference_orientation
STR_table$length_characters <- str_length(STR_table$reference_motif_reference_orientation)

#group everything larger than 6 (the VNTRs)
STR_table$length_characters[STR_table$length_characters > 6] <- ">6"

#order for plot
limits <- c("1","3", "4", "5", "6", ">6")

# Mutate the type for POLG gene
STR_table <- STR_table %>%
  mutate(type = ifelse(gene == 'POLG', 'Coding', type))

# Plot
ggplot(STR_table, aes(x = type, fill = factor(as.character(length_characters), levels = limits))) +
  geom_bar() +
  labs(x = "Genomic region type",
       y = "Number of loci") +
  scale_fill_manual(values = blue_palette, name = "Motif length", limits = limits) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#Further determining by coding subtype; taken from literature review
STR_table$type2 <- ifelse(grepl("polyglutamine", STR_table$mechanism_detail, ignore.case = TRUE), "PolyQ",
                          ifelse(grepl("polyalanine", STR_table$mechanism_detail, ignore.case = TRUE), "PolyA", STR_table$type))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'NIPA1', 'PolyA', type2))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'POLG', 'PolyQ', type2))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'HTT', 'PolyQ', type2))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'PRDM12', 'PolyA', type2))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'THAP11', 'PolyQ', type2))

#subset for pie chart
STR_table_coding <- subset(STR_table, type2 == "Coding" | type2 == "PolyA" | type2 == "PolyQ")

STR_table_coding <- STR_table_coding %>%
  mutate(type2 = ifelse(gene == 'ZFHX3', 'PolyG', type2))

STR_table_coding <- STR_table_coding %>%
  mutate(type2 = ifelse(gene == 'COMP', 'PolyD', type2))

STR_table_coding <- STR_table_coding %>%
  mutate(type2 = ifelse(gene == 'JPH3', 'Other', type2))

STR_table_coding <- STR_table_coding %>%
  mutate(type2 = ifelse(gene == 'PRNP', 'Other', type2))

STR_table_coding <- STR_table_coding %>%
  mutate(type2 = ifelse(gene == 'VWA1', 'Other', type2))

#take out MUC1
STR_table_coding <- subset(STR_table_coding, gene != "MUC1")

# further detail is desired
# STR_table_coding <- STR_table_coding %>%
#   mutate(type2 = ifelse(gene == 'JPH3', 'Coding:Poly-A or PolyL', type2))
#
# STR_table_coding <- STR_table_coding %>%
#   mutate(type2 = ifelse(gene == 'PRNP', 'Coding:P, G, Q', type2))
# further detail is desired
# STR_table_coding <- STR_table_coding %>%
#   mutate(type2 = ifelse(gene == 'JPH3', 'Coding:Poly-A or PolyL', type2))
#
# STR_table_coding <- STR_table_coding %>%
#   mutate(type2 = ifelse(gene == 'PRNP', 'Coding:P, G, Q', type2))

# Generate a qualitative color palette with 5 colors
qual_palette <- brewer.pal(5, "Set2")

# Display the color palette
qual_palette

# Plot the pie chart with labels
type_counts <- sort(table(STR_table_coding$type2), decreasing = TRUE)

waffle(type_counts, rows = 5, size = 2, colors = qual_palette,
       title = "Coding Subtypes", flip = TRUE, reverse = TRUE, legend = "bottom")

### Figure 1B
STR_table_adjusted <- STR_table %>%
  mutate(
    benign_min = coalesce(benign_min, benign_max),
    pathogenic_max = coalesce(pathogenic_max, pathogenic_min)
  )

STR_table_adjusted <- STR_table_adjusted %>%
  mutate(
    benign_min = ifelse(gene == "MARCHF6", 0, benign_min),
    benign_max = ifelse(gene == "MARCHF6", 0, benign_max)
  )

# PMID: 20399836
STR_table_adjusted <- STR_table_adjusted %>%
  mutate(
    pathogenic_min = ifelse(gene == "POLG", 6, pathogenic_min),
    pathogenic_max = ifelse(gene == "POLG", 14, pathogenic_max)
  )

STR_table_adjusted$norm_min_bp = STR_table_adjusted$benign_min * STR_table_adjusted$motif_len
STR_table_adjusted$norm_max_bp = STR_table_adjusted$benign_max * STR_table_adjusted$motif_len
STR_table_adjusted$int_min_bp = STR_table_adjusted$intermediate_min * STR_table_adjusted$motif_len
STR_table_adjusted$int_max_bp = STR_table_adjusted$intermediate_max * STR_table_adjusted$motif_len
STR_table_adjusted$path_min_bp = STR_table_adjusted$pathogenic_min * STR_table_adjusted$motif_len
STR_table_adjusted$path_max_bp = STR_table_adjusted$pathogenic_max * STR_table_adjusted$motif_len


### presumably, I'm going to remove untrustworthy data

genes_to_remove <- c("DMD", "ZIC3", "TNR6CA", "YEATS2", "TBX1", "POLG", "NAXE", "RAI1")

STR_table_adjusted <- STR_table_adjusted %>%
  filter(!gene %in% genes_to_remove)

ggplot(STR_table_adjusted, aes(x = gene)) +
  theme_minimal() +
  geom_linerange(aes(ymin = norm_max_bp, ymax = path_max_bp, color = "gray"),
                 linewidth = 0.5, linetype = "dotted", alpha = 0.7) +
  geom_linerange(aes(ymin = path_min_bp, ymax = path_max_bp, color = 'Pathogenic'), linewidth = 1.5) +
  geom_point(data = subset(STR_table_adjusted, int_min_bp == int_max_bp), aes(y = int_min_bp), shape = 16, color= "#E7B800", size = 2) +
  geom_point(data = subset(STR_table_adjusted, path_min_bp == path_max_bp), aes(y = path_min_bp), shape = 16, color = "#FC4E07", size = 2) +
  geom_point(data = subset(STR_table_adjusted, norm_min_bp == norm_max_bp), aes(y = norm_min_bp), shape = 16, color= "#00AFBB", size = 2) +
  geom_linerange(aes(ymin = norm_min_bp+0.9, ymax = norm_max_bp, color = 'Normal'), linewidth = 1.5) +
  geom_linerange(aes(ymin = int_min_bp, ymax = int_max_bp, color = 'Intermediate*'), linewidth = 1.5) +
  scale_y_continuous(trans = 'log10') +
  scale_x_discrete(name = 'Gene') +
  scale_colour_manual(values = c('Normal' = '#00AFBB', 'Intermediate*' = '#E7B800', 'Pathogenic' = '#FC4E07'),
                      breaks = c('Normal', 'Intermediate*', 'Pathogenic'),
                      name = 'Allele size') +
  theme(panel.grid.major.x = element_line(color = 'lightgrey', linetype = 'longdash'),
        axis.title.x = element_text(size=20)) +
  labs(y = "Allele size in base pairs") +
  coord_flip()

#panels put into same figure in Adobe Illustrator
