library(RColorBrewer)

# Generate 5 shades of blue
blue_palette <- colorRampPalette(c("#ADD8E6", "#6495ED", "#4169E1", "#0000FF", "#00008B"))(5)

# Count the number of characters in pathogenic_motif_reference_orientation
STR_table$length_characters <- str_length(STR_table$pathogenic_motif_reference_orientation)

#group everything larger than 6 (the VNTRs)
STR_table$length_characters[STR_table$length_characters > 6] <- ">6"

#order for plot
limits <- c("3", "4", "5", "6", ">6")

# Mutate the type for POLG gene
STR_table <- STR_table %>%
  mutate(type = ifelse(gene == 'POLG', 'Coding', type))

# Plot
ggplot(STR_table, aes(x = type, fill = as.character(length_characters))) +
  geom_bar() +
  labs(title = "Genomic Regions and Motif Length",
       x = "Genomic Region Type",
       y = "Count") +
  scale_fill_manual(values = blue_palette, name = "Motif Length", limits = limits) +
  theme_minimal()

#Further determining by coding subtype; taken from literature review
STR_table$type2 <- ifelse(tolower(trimws(STR_table$Mechanism)) == "polyglutamine", "PolyQ",
                          ifelse(tolower(trimws(STR_table$Mechanism)) == "polyalanine", "PolyA", STR_table$type))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'NIPA1', 'PolyA', type2))

STR_table <- STR_table %>%
  mutate(type2 = ifelse(gene == 'POLG', 'PolyQ', type2))

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

# further detail is desired
# STR_table_coding <- STR_table_coding %>%
#   mutate(type2 = ifelse(gene == 'JPH3', 'Coding:Poly-A or PolyL', type2))
#
# STR_table_coding <- STR_table_coding %>%
#   mutate(type2 = ifelse(gene == 'PRNP', 'Coding:P, G, Q', type2))

pie_data <- table(STR_table_coding$type2)

# Generate 5 shades of blue
green_palette <- colorRampPalette(c("lightgreen", "darkgreen"))(5)

# Plot the pie chart with labels
pie(pie_data, col = green_palette, main = "Coding Subtypes", labels = paste(names(pie_data),"(",pie_data,")"))
