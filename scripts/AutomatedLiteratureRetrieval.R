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
data <- read.csv('/Users/quinlan/Documents/Git/STRchive_manuscript/data/STR-disease-loci.csv',
                 stringsAsFactors = FALSE)

# Filter out NA values, if there are any, from the 'gene' column
filtered_data <- data[!is.na(data$gene), ]

# Extract STRchive gene names into a list
gene_list <- as.character(unique(filtered_data$gene))

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

# strings for queries
consolidated_strings <- gene_info %>%
  group_by(hgnc_symbol) %>%
  summarize(consolidated_strings = paste(unique(c(hgnc_symbol, external_synonym)), collapse = ' OR ')) %>%
  pull(consolidated_strings)

#substitute to avoid unrelated PMIDs
consolidated_strings <- gsub("BMD", "Becker muscular dystrophy", consolidated_strings)

#where results will be stored
base_directory <- '/Users/quinlan/Documents/Git/STRchive_manuscript/data/literature/'

# function to perform the pubmed query
# Function printout includes gene name and if there are results, confirms
# that a file has been created
perform_pubmed_query <- function(gene_info) {
  file_paths <- list()  # Initialize the list to store all publications
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    cat("Processing gene:", gene_name, "\n")
    # Use str_detect to check if gene_name is present in each consolidated string
    idx <- str_detect(consolidated_strings, regex(gene_name, ignore_case = TRUE))

    # Put together the relevant consolidated strings
    or_terms <- paste(consolidated_strings[idx], collapse = ' OR ')
    # Splitting the or_terms into individual terms
    individual_terms <- unlist(strsplit(or_terms, " OR "))

    # Adding [Title/Abstract] to each term
    joined_terms <- paste0(individual_terms, "[Title/Abstract]")
    #joined_terms <- paste0('(', paste(or_terms, collapse = '[Title/Abstract] OR '), ')[Title/Abstract]')
    # Construct the query with organized or_terms
    query <- paste0('("repeat expansion"[Title/Abstract] OR "tandem repeat"[Title/Abstract] OR "repeat expansions"[Title/Abstract] OR "tandem repeats"[Title/Abstract] OR "repeat sequence"[Title/Abstract] OR "repeat sequences"[Title/Abstract] OR "repeat length"[Title/Abstract] OR "repeat lengths"[Title/Abstract] OR "expansion"[Title] OR "expansions"[Title] OR "repeats"[Title]) AND (', paste(joined_terms, collapse = " OR "),') AND "English"[Language] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "diseases"[Title/Abstract] OR "disorders"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "syndromes"[Title/Abstract] OR "patient"[Title/Abstract] OR "patients"[Title/Abstract] OR "proband"[Title/Abstract] OR "probands"[Title/Abstract]) AND ("journal article"[Publication Type] OR "letter"[Publication Type] or "Case Reports"[Publication Type]) NOT "review"[Publication Type]')


    # Clean up any unnecessary slashes from the query
    query <- gsub("  ", " ", query)  # Remove double spaces
    print(query)
    gene_name <- gsub('"', '', gene_name)
    out_file <- paste0(base_directory, gene_name)

    # Include a separator ("/") between base_directory and gene_name
    # Modify dest_file_prefix to include the full file path
    out.A <- batch_pubmed_download(pubmed_query_string = query,
                                   format = "medline",
                                   batch_size = 10000,
                                   dest_file_prefix = out_file,
                                   encoding = "ASCII")
    # the function adds 01.txt so, gotta fix that here
    out_file <- paste0(base_directory, gene_name, "01.txt")
    print(out_file)

    # Check if the file was created successfully
    cat("Full file path:", out_file, "\n")
    if (file.exists(out_file)) {
      cat("File exists.\n")
      file_paths[[gene_name]] <- out_file
    } else {
      cat("Error: File not found -", out_file, "\n")
    }
  }

  return(file_paths)
}

file_paths <- perform_pubmed_query(gene_info)

#creating list for publications
all_publications <- list()

# Assuming file_paths is a list of file paths
#let's get the files into a dataframe in R
for (gene_name in names(file_paths)) {
  # Append "01.txt" to the file path
  file_path <- paste0(file_paths[[gene_name]])

  tryCatch({
    # Read the file into a character vector
    current_publications <- readLines(file_path)

    # Append to the overall list
    all_publications[[gene_name]] <- current_publications
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    # Append an empty character vector in case of an error
    all_publications[[gene_name]] <- character(0)
  })
}

#let's extract publication information through a function
extract_pub_info <- function(medline_data_list, gene_name) {
  # Combine the list of XML strings into a single string
  medline_string <- paste(medline_data_list, collapse = "")

  # Use regular expressions to extract PMID and publication years
  # Extract PMIDs
  pmids <- str_extract_all(medline_string, "(?<=PMID- )\\d+")[[1]]

  # Extract Publication Years
  publication_years <- str_extract_all(medline_string, "(?<=EDAT- )\\d{4}")[[1]]

  # Ensure both vectors have the same length
  length_diff <- length(pmids) - length(publication_years)
  if (length_diff > 0) {
    publication_years <- c(publication_years, rep(NA, length_diff))
  } else if (length_diff < 0) {
    pmids <- c(pmids, rep(NA, -length_diff))
  }

  # Create a dataframe with gene_name, PMID, and PublicationYear
  pub_info_df <- data.frame(GeneName = rep(gene_name, length(pmids)),
                            PMID = pmids,
                            PublicationYear = publication_years,
                            stringsAsFactors = FALSE)

  return(pub_info_df)
}


# empty list for the pubication info
pub_info_list <- list()

# Loop through each gene_name in all_publications
# print to make sure each gene is processed
for (gene_name in names(all_publications)) {
  # Get the list of XML data for the current gene_name
  medline_data_list <- all_publications[[gene_name]]
  print(gene_name)
  # Extract publication information using the function
  pub_info_df <- extract_pub_info(medline_data_list, gene_name)

  # Append the results to the list
  pub_info_list[[gene_name]] <- pub_info_df
}

# Combine all the dataframes into a single dataframe
all_pub_info_df <- do.call(rbind, pub_info_list)

# Get the current date for a unique name
current_date <- format(Sys.Date(), "%Y%m%d")

# Concatenate the date to the file name
file_name <- paste0("/Users/quinlan/Documents/Git/STRchive_manuscript/data/all_pub_info_", current_date, ".tsv")

# Write the table with the updated file name
write.table(all_pub_info_df, file_name, sep = "\t", quote = FALSE, row.names = FALSE)

