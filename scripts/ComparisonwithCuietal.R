library(dplyr)

normalise_str <- function(in_dna) {
  # Split the input by comma if it contains one
  if (grepl(",", in_dna)) {
    in_dna <- unlist(strsplit(in_dna, ",", fixed = TRUE))
  } else {
    in_dna <- list(in_dna)  # Convert single string to list
  }

  result <- character(length(in_dna))  # Initialize result vector

  for (i in seq_along(in_dna)) {
    dna <- in_dna[i]

    if (is.null(dna) || length(dna) == 0) {
      result[i] <- ''  # Return empty string for NULL or empty input
    } else {
      # Generate all circular permutations of input sequence
      all_possible <- sapply(0:(nchar(dna)-1), function(j)
        paste0(substr(dna, j+1, nchar(dna)), substr(dna, 1, j)))

      # Sort permutations alphabetically and return the first
      result[i] <- sort(all_possible)[1]
    }
  }

  return(paste(result, collapse = ","))  # Combine result elements with commas
}

### Uploading Data
STR_table <- read.csv('/Users/quinlan/Documents/Git/STRchive/data/STR-disease-loci.csv',
                      stringsAsFactors = FALSE)

# change to match gnomAD and recognize different loci within same gene
STR_table$gene[STR_table$gene == "ARX" & STR_table$stop_hg38 == 25013697] <- "ARX_1"
STR_table$gene[STR_table$gene == "ARX" & STR_table$stop_hg38 == 25013565] <- "ARX_2"
STR_table$gene[STR_table$gene == "HOXA13" & STR_table$stop_hg38 == 27199966] <- "HOXA13_1"
STR_table$gene[STR_table$gene == "HOXA13" & STR_table$stop_hg38 == 27199861] <- "HOXA13_2"
STR_table$gene[STR_table$gene == "HOXA13" & STR_table$stop_hg38 == 27199732] <- "HOXA13_3"


## smaller dataset with more relevant columns
STR_table_clean <-subset(STR_table, select=c("disease_id", "gene",
                                             "reference_motif_reference_orientation", "Inheritance",
                                             "type", "normal_min", "normal_max",
                                             "intermediate_min", "intermediate_max",
                                             "pathogenic_min", "pathogenic_max",
                                             "repeatunitlen", "age_onset_min","typ_age_onset_min",
                                             "age_onset_max", "typ_age_onset_max", "novel", "pathogenic_motif_reference_orientation"))

STR_table_clean$repeatunit_path_normalized <- sapply(STR_table_clean$pathogenic_motif_reference_orientation, function(x) normalise_str(as.character(x)))


STR_table_clean$repeatunit_ref_normalized <- sapply(STR_table_clean$reference_motif_reference_orientation, function(x) normalise_str(as.character(x)))

### Cui data
# total
# 338,963
# African population
# 64,786
# East Asian population
# 1,880
# European population
# 205,116
# Hispanic population
# 49,822
# South Asian population
# 4,245

European = read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/Primary_European.txt',
                          sep = '\t', stringsAsFactors = FALSE)

Hispanic = read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/Primary_Hispanic.txt',
                    sep = '\t', stringsAsFactors = FALSE)

SouthAsian = read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/Primary_South_Asian.txt',
                    sep = '\t', stringsAsFactors = FALSE)

EastAsian = read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/Primary_East_Asian.txt',
                    sep = '\t', stringsAsFactors = FALSE)

African = read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/Primary_African.txt',
                     sep = '\t', stringsAsFactors = FALSE)

total <- rbind(European, Hispanic, SouthAsian, EastAsian, African)

disease = read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/CuiDiseaseTable.tsv',
                   sep = '\t', stringsAsFactors = FALSE, header = TRUE)

colnames(disease)[6] = "gene"

disease$gene[disease$gene == "ENSG00000225885"] <- "SAMD12"
disease$gene[disease$gene == "ARX" & disease$TR_Locus == "chrX:25013649-25013697"] <- "ARX_1"
disease$gene[disease$gene == "ARX" & disease$TR_Locus == "chrX:25013529-25013565"] <- "ARX_2"
disease$gene[disease$gene == "HOXA13" & disease$TR_Locus == "chr7:27199924-27199966"] <- "HOXA13_1"
disease$gene[disease$gene == "HOXA13" & disease$TR_Locus == "chr7:27199825-27199861"] <- "HOXA13_2"
disease$gene[disease$gene == "HOXA13" & disease$TR_Locus == "chr7:27199678-27199732"] <- "HOXA13_3"

disease_clean <-subset(disease, select=c("TR_ID", "Motif", "gene"))

STRchive_and_Cui_disease <- merge(STR_table_clean,disease_clean,by="gene")

TR_ids <- unique(STRchive_and_Cui_disease$TR_ID)

total <- total %>%
  filter(TR_ID %in% TR_ids)
