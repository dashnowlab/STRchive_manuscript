library(biomaRt)
library(rentrez)
library(dplyr)
# solution for potential error based on library versions
#devtools::install_version("dbplyr", version = "2.3.4")
#library(dbplyr)
library(easyPubMed)
library(stringr)


### Data Setup
# change to STRchive directory
data <- read.csv('/Users/quinlan/Documents/Git/STRchive/data/STR-disease-loci.csv',
                 stringsAsFactors = FALSE)

# Filter out NA values, if there are any, from the 'gene' column
filtered_data <- data[!is.na(data$gene), ]

# Extract STRchive gene names into a list
gene_list <- as.character(filtered_data$gene)

### Use biomaRt to get gene synonyms
# This is necessary because gene names have evolved over time and can differ
# by publication
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                host = "https://www.ensembl.org")

# Retrieve the HGNC symbol and Gene Synonym for the genes in the list
# Check if the mart object is created
if (exists("mart")) {
  # Retrieve the HGNC symbol and Gene Synonym for the genes in the list
  gene_info <- getBM(attributes = c("hgnc_symbol", "external_synonym"),
                     filters = "hgnc_symbol",
                     values = gene_list,
                     mart = mart)

} else {
  cat("Failed to create the mart object. Check your internet connection and try again.")
}

# Replace missing values in 'external_synonym' with 'hgnc_symbol'
# This is because some gene_names have no synonyms, so it's easier to group this way
gene_info$external_synonym[is.na(gene_info$external_synonym) |
                             gene_info$external_synonym == ""] <- gene_info$hgnc_symbol[is.na(
                               gene_info$external_synonym) | gene_info$external_synonym == ""]

# curated list of synonyms to exclude based on non-specific results
excluded_synonym_list <- c("B37", "MHP", "MED", "DM", "DM1", "FA", "GAC", "SPD",
                           "PRP", "A1", "CCD", "PHP", "VCF")

#remove from gene_info
gene_info <- gene_info %>%
  filter(!grepl(paste(excluded_synonym_list, collapse = '|'), external_synonym))

# adding quotation marks to terms for query
gene_info <- gene_info %>%
  mutate(across(c(external_synonym, hgnc_symbol), ~sprintf('"%s"', .)))

# adding gene synonym missed by biomaRt
gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"FMR1\"", external_synonym = "\"FMR-1\""))

#because pubmed hates slashes and these use slashes
gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"NUTM2B-AS1\"", external_synonym = "\"LOC642361/NUTM2B-AS1\""))
gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"ATXN10\"", external_synonym = "\"SCA10\""))

# collapse to gene names
gene_names <- unique(gene_info$hgnc_symbol)

#empty list for publications
publications <- list()

# strings for queries
consolidated_strings <- gene_info %>%
  group_by(hgnc_symbol) %>%
  summarize(consolidated_strings = paste(unique(c(hgnc_symbol, external_synonym)), collapse = ' OR ')) %>%
  pull(consolidated_strings)

#substitute to avoid unrelated PMIDs
consolidated_strings <- gsub("BMD", "Becker muscular dystrophy", consolidated_strings)

#where results will be stored
base_directory <- '/Users/quinlan/Documents/Git/STRchive/data/'


