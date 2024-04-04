
### Libraries for data and Figures 2, 3,4, and 5
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(broom)
library(ggridges)
library(DescTools)
library(Biostrings)
library(plotly)
library(ggbreak)
library(patchwork)

### Functions
#normalise TR motifs; can be string with commas
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

# normalise TR motifs, no strings, as in genotype call
normalise_call <- function(in_dna) {
  if (is.null(in_dna) || length(in_dna) == 0) {
    return('')
  }
  # Generate all circular permutations of input sequence
  all_possible <- sapply(0:(nchar(in_dna)-1), function(i) paste0(substr(in_dna, i+1, nchar(in_dna)), substr(in_dna, 1, i)))
  # Sort permutations alphabetically and return the first
  return(sort(all_possible)[1])
}

### Uploading Data
# gnomAD
gnomADSTRcalls = read.csv('/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__including_all_age_and_pcr_info__2023_06_28.tsv',
                          sep = '\t', stringsAsFactors = FALSE)

STR_table <- read.csv('/Users/quinlan/Documents/Git/STRchive/data/STR-disease-loci.csv',
                      stringsAsFactors = FALSE)

# change to match gnomAD and recognize different loci within same gene
STR_table$gene[STR_table$gene == "C9orf72"] <- "C9ORF72"
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

#rename for easier merge
colnames(gnomADSTRcalls)[1] = "gene"
# split into the two motifs
gnomADSTRcalls[c("Allele1Motif", "Allele2Motif")] <- as.data.frame(do.call(rbind, strsplit(gnomADSTRcalls$Motif, "/")), stringsAsFactors = FALSE)

# split confidence interval into bounds for each allele
gnomADSTRcalls[c("Allele1LowerBound", "Allele1HigherBound", "Allele2LowerBound", "Allele2HigherBound")] <- as.data.frame(do.call(rbind, strsplit(gsub("/", "-", gnomADSTRcalls$GenotypeConfidenceInterval), "-")), stringsAsFactors = FALSE)
gnomADSTRcalls[c("Allele1LowerBound", "Allele1HigherBound", "Allele2LowerBound", "Allele2HigherBound")] <- lapply(gnomADSTRcalls[c("Allele1LowerBound", "Allele1HigherBound", "Allele2LowerBound", "Allele2HigherBound")], as.integer)


total <- merge(STR_table_clean,gnomADSTRcalls,by="gene")

total$motif_norm <- sapply(total$Allele2Motif, function(x) normalise_call(as.character(x)))
total$motif_norm_small <- sapply(total$Allele1Motif, function(x) normalise_call(as.character(x)))

