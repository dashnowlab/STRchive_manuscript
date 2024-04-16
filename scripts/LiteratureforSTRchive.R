extract_citation_info <- function(medline_data_list, gene_name) {
  # Combine the list of XML strings into a single string
  medline_string <- paste(medline_data_list, collapse = "")

  # Use regular expressions to extract PMID, publication years, and titles
  # Extract PMIDs
  pmids <- str_extract_all(medline_string, "(?<=PMID- )\\d+")[[1]]
  print(pmids)
  # Extract Publication Dates
  publication_dates <- str_extract_all(medline_string, "(?<=DP  - )\\d+")[[1]]
  print(publication_dates)
  # Extract Titles
  title <- str_extract_all(medline_string, "(?<=TI  - ).+?(?=\\.|\\?)")[[1]]
  print(title)

  # Ensure all vectors have the same length
  length_diff <- length(pmids) - length(publication_dates)
  if (length_diff > 0) {
    publication_dates <- c(publication_dates, rep(NA, length_diff))
  } else if (length_diff < 0) {
    pmids <- c(pmids, rep(NA, -length_diff))
  }

  length_diff <- length(pmids) - length(title)
  if (length_diff > 0) {
    title <- c(title, rep(NA, length_diff))
  } else if (length_diff < 0) {
    pmids <- c(pmids, rep(NA, -length_diff))
  }

  # Create a dataframe with gene_name, PMID, PublicationYear, and Title
  pub_info_df <- data.frame(GeneName = rep(gene_name, length(pmids)),
                            PMID = pmids,
                            PublicationDate = publication_dates,
                            Title = title,
                            stringsAsFactors = FALSE)

  return(pub_info_df)
}


for (gene_name in names(all_publications)) {
  # Get the list of XML data for the current gene_name
  medline_data_list <- all_publications[[gene_name]]
  print(gene_name)
  # Extract publication information using the function
  pub_info_df <- extract_citation_info(medline_data_list, gene_name)

  # Append the results to the list
  pub_info_list[[gene_name]] <- pub_info_df
}

# Combine all the dataframes into a single dataframe
all_pub_info_df <- do.call(rbind, pub_info_list)
