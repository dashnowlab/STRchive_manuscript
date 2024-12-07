library(ggplot2)
library(ggrepel)
library(scales)
library(dplyr)

#can upload data from file path or use all_pub_info_df
# either can be GeneNamerated from AutomatedLiteratureRetrieval.R
# file_path <- "/Users/quinlan/Documents/Git/STRchive_manuscript/data/all_pub_info_20240430.tsv"
 #
 #all_pub_info_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


### Adding more curated data to use in visualization of PMIDs
# first publication by PMID
firstpub <- data.frame(
  stringsAsFactors = FALSE,
  GeneName = c("ABCD3","AFF2","AFF3","AR","ARX",
               "ATN1","ATXN1","ATXN10","ATXN2","ATXN3","ATXN7",
               "ATXN8OS","BEAN1","C9orf72","CACNA1A","CBL","COMP",
               "CSTB","DAB1","DIP2B","DMD","DMPK","EIF4A3","FGF14",
               "FMR1","FOXL2","FXN","GIPC1","GLS","HOXA13","HOXD13",
               "HTT","JPH3","LRP12","MARCHF6","NIPA1","NOP56",
               "NOTCH2NLC","NUTM2B-AS1","PABPN1","PHOX2B","POLG","PRNP",
               "RAPGEF2","RFC1","RILPL1","RUNX2","SAMD12",
               "STARD7","TBP","TCF4","THAP11","TNRC6A","VWA1","XYLT1",
               "YEATS2","ZFHX3","ZIC2","ZIC3","ZNF713","CNBP",
               "PPP2R2B", "SOX3", "TBX1", "PRDM12", "MUC1", "NAXE",
               "pre-MIR7-2", "RAI1", "TAF1"),
  PMID = c(39068203L, 8334699L,24763282L,1461383L,
           19587282L,7842016L,7951322L,11017075L,8896556L,
           7874163L,8908515L,10192387L,19878914L,22154785L,9371901L,
           7603564L,9887340L,9126745L,28686858L,17236128L,
           27417533L,1310900L,24360810L,36516086L,1605194L,
           12529855L,8596916L,32413282L,30970188L,15385446L,8614804L,
           8401589L,11694876L,31332380L,31664039L,31286297L,
           21683323L,31332380L,31332380L,9462747L,12640453L,
           22963882L,1683708L,29507423L,31230722L, 35148830L,
           26220009L,29507423L,31664034L,10484774L,23185296L,24677642L,
           29507423L,33559681L,30554721L,31539032L,38035881L,
           11285244L,20452998L,25196122L,11486088L,10581021L,12428212L,
           19948535L,26005867L, 23396133L,39455596L, 38714869L, 37994247L, 29229810)
)


# Merge based on 'GeneName' and 'PMID'
merged_df <- merge(firstpub, all_pub_info_df[, c("GeneName", "PMID", "PublicationYear")],
                   by = c("GeneName", "PMID"),
                   all.x = TRUE)  # 'all.x = TRUE' retains all rows from 'firstpub'

# what doesn't come up in our data bc of data types (2) or abstract not indexed (1)
merged_df$PublicationYear[merged_df$PMID == 7951322] <- 1994
# ATXN1
merged_df$PublicationYear[merged_df$PMID == 11486088] <- 2001
# CNBP
merged_df$PublicationYear[merged_df$PMID == 10581021] <- 1999
# PPP2R2B

# what doesn't come up because of key words not in abstract
merged_df$PublicationYear[merged_df$PMID == 26005867] <- 2015
# PRDM12
merged_df$PublicationYear[merged_df$PMID == 12428212] <- 2002
# SOX3
merged_df$PublicationYear[merged_df$PMID == 19948535] <- 2010
# TBX1
merged_df$PublicationYear[merged_df$PMID == 1310900] <- 1992
#DMPK
merged_df$PublicationYear[merged_df$PMID == 12529855] <- 2003
# FOXL2
merged_df$PublicationYear[merged_df$PMID == 8614804] <- 1996
# HOXD13

merged_df$PublicationYear[merged_df$PMID == 38035881] <- 2023
# fixing pubmed index bug that has wrong year for this PMID

merged_df$PublicationYear[merged_df$PMID == 38714869] <- 2024

#note: The final HOXA13 locus was characterized in the 2004 PMID here,
# but a 2000 manuscript identifed an expansion in the third tract (PMID:10839976)

merged_df <- merged_df %>%
  rename(EarliestPublicationYear = PublicationYear)

all_pub_info_df$PublicationYear <- as.numeric(all_pub_info_df$PublicationYear)

# Join 'earliest_years' with 'all_pmids' based on 'GeneName_name'
merged_df <- inner_join(merged_df, all_pub_info_df, by = "GeneName")

# Filter out rows where Publicationyear is earlier than the EarliestPublicationYear
filtered_df <- merged_df %>%
  filter(PublicationYear >= EarliestPublicationYear)

# ensure numeric for comparison
filtered_df$PMID.y <- as.numeric(filtered_df$PMID.y)

# Find the PMID.x values that are not present in PMID.y, add'em in
missing_pmids <- anti_join(filtered_df, filtered_df %>% dplyr::select(PMID.y),
                           by = c("PMID.x" = "PMID.y"))

# Create a new dataframe with the missing PMIDs and EarliestPublicationYear
# this should be (and currently is) the 6 we didn't find from earlier
new_rows <- missing_pmids %>%
  mutate(PMID.y = PMID.x, PublicationYear = EarliestPublicationYear) %>%
  distinct()

# Append the new rows to the original dataframe
filtered_df <- bind_rows(filtered_df, new_rows)

# evidence by number of independent observations
evidence <- data.frame(
  GeneName = c("ABCD3","AFF2","AFF3","AR","ARX",
               "ATN1","ATXN1","ATXN10","ATXN2","ATXN3","ATXN7",
               "ATXN8OS","BEAN1","C9orf72","CACNA1A","CBL","COMP",
               "CSTB","DAB1","DIP2B","DMD","DMPK","EIF4A3","FGF14",
               "FMR1","FOXL2","FXN","GIPC1","GLS","HOXA13","HOXD13",
               "HTT","JPH3","LRP12","MARCHF6","NIPA1","NOP56",
               "NOTCH2NLC","NUTM2B-AS1","PABPN1","PHOX2B","POLG","PRNP",
               "RAPGEF2","RFC1","RILPL1","RUNX2","SAMD12",
               "STARD7","TBP","TCF4","THAP11","TNRC6A","VWA1","XYLT1",
               "YEATS2","ZFHX3","ZIC2","ZIC3","ZNF713","CNBP",
               "PPP2R2B", "SOX3", "TBX1", "PRDM12", "MUC1", "NAXE",
               "pre-MIR7-2", "RAI1", "TAF1"),
  ind_obs = c("< 10","100 > x > 50","< 10","> 100","50 > x > 10","100 > x > 50",
              "> 100","100 > x > 50","> 100","> 100","> 100","> 100","> 100",
              "> 100","> 100","< 10","< 10","50 > x > 10","< 10","< 10","> 100","> 100",
              "50 > x > 10","> 100","50 > x > 10","< 10","< 10","< 10","> 100",
              "100 > x > 50","100 > x > 50","> 100","100 > x > 50","> 100","> 100","< 10","> 100",
              "> 100","> 100","> 100","< 10","< 10","> 100","50 > x > 10","< 10","> 100","< 10",
              "50 > x > 10","100 > x > 50","< 10","> 100","< 10","50 > x > 10","< 10","50 > x > 10",
              "< 10","< 10","> 100","> 100","50 > x > 10","> 100","50 > x > 10","< 10","> 100",
              "50 > x > 10","> 100", "< 10", "50 > x > 10","< 10","> 100"))

# summary of PMIDS and years
summary_data <- filtered_df %>%
  group_by(GeneName) %>%
  summarize(MinPublicationYear = min(PublicationYear, na.rm = TRUE), MaxPublicationYear = max(PublicationYear, na.rm = TRUE),
            TotalPMIDs = n_distinct(PMID.y))

#merge data for final plot
merged_data <- merge(evidence, summary_data, by = "GeneName", all.x = FALSE)

#create ordered labels for legend
labels <- c("< 10" = "< 10", "50 > x > 10" = "10 < x < 50",
            "100 > x > 50" = "50 < x < 100", "> 100" = "> 100")

limits <- c("< 10", "50 > x > 10", "100 > x > 50", "> 100")

#final plot
ggplot(merged_data, aes(x = MinPublicationYear, y = TotalPMIDs, label = GeneName)) +
  geom_jitter(aes(size = ind_obs, color = ind_obs), width = 0.5) +
  geom_text_repel(data = subset(summary_data, !duplicated(GeneName)), aes(label = GeneName),
                  box.padding = 0.5, segment.color = "grey50", segment.size = 0.2,
                  nudge_y = 0.1, size = 5) +
  labs(title = "Minimum Publication Year by Total Number of Unique PMIDs",
       x = "Minimum Publication Year",
       y = "Total Number of Unique PMIDs (log scale)",
       size = "Ind. Obs.") +
  scale_y_log10() +
  scale_size_manual(values = c("< 10" = 2, "50 > x > 10" = 3,
                               "100 > x > 50" = 4, "> 100" = 6),
                    limits = limits,
                    labels = labels) +
  scale_color_manual(values = c("< 10" = "gray", "50 > x > 10" = "#DDCC77",
                                "100 > x > 50" = "#88CCEE", "> 100" = "#332288"),
                     limits = limits, labels = labels,
                     name = "Ind. Obs.") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

# what doesn't come up in our data bc of data types (2/6) or abstract not indexed (1/6)
merged_df$PublicationYear[merged_df$PMID == 7951322] <- 1994
# ATXN1
merged_df$PublicationYear[merged_df$PMID == 11486088] <- 2001
# CNBP
merged_df$PublicationYear[merged_df$PMID == 10581021] <- 1999
# PPP2R2B

# what doesn't come up because of key words not in abstract
merged_df$PublicationYear[merged_df$PMID == 26005867] <- 2015
# PRDM12
merged_df$PublicationYear[merged_df$PMID == 12428212] <- 2002
# SOX3
merged_df$PublicationYear[merged_df$PMID == 19948535] <- 2010
# TBX1

merged_df <- merged_df %>%
  rename(EarliestPublicationYear = PublicationYear)

all_pub_info_df$PublicationYear <- as.numeric(all_pub_info_df$PublicationYear)

# Join 'earliest_years' with 'all_pmids' based on 'GeneName_name'
merged_df <- inner_join(merged_df, all_pub_info_df, by = "GeneName")

# Filter out rows where Publicationyear is earlier than the EarliestPublicationYear
filtered_df <- merged_df %>%
  filter(PublicationYear >= EarliestPublicationYear)

# ensure numeric for comparison
filtered_df$PMID.y <- as.numeric(filtered_df$PMID.y)

# Find the PMID.x values that are not present in PMID.y, add'em in
missing_pmids <- anti_join(filtered_df, filtered_df %>% dplyr::select(PMID.y),
                           by = c("PMID.x" = "PMID.y"))

# Create a new dataframe with the missing PMIDs and EarliestPublicationYear
# this should be (and currently is) the 6 we didn't find from earlier
new_rows <- missing_pmids %>%
  mutate(PMID.y = PMID.x, PublicationYear = EarliestPublicationYear) %>%
  distinct()

# Append the new rows to the original dataframe
filtered_df <- bind_rows(filtered_df, new_rows)

# evidence by number of independent observations

# summary of PMIDS and years
summary_data <- filtered_df %>%
  group_by(GeneName) %>%
  summarize(MinPublicationYear = min(PublicationYear, na.rm = TRUE), MaxPublicationYear = max(PublicationYear, na.rm = TRUE),
            TotalPMIDs = n_distinct(PMID.y))

#merge data for final plot
merged_data <- merge(evidence, summary_data, by = "GeneName", all.x = FALSE)

#create ordered labels for legend
labels <- c("< 10" = "< 10", "50 > x > 10" = "10 < x < 50",
            "100 > x > 50" = "50 < x < 100", "> 100" = "> 100")

limits <- c("< 10", "50 > x > 10", "100 > x > 50", "> 100")

#final plot
library(ggplot2)

ggplot(merged_data, aes(x = MinPublicationYear, y = TotalPMIDs, label = GeneName)) +
  theme_minimal() +
  geom_jitter(aes(size = ind_obs, color = ind_obs), width = 0.5) +
  geom_text_repel(data = subset(summary_data, !duplicated(GeneName)), aes(label = GeneName),
                  box.padding = 0.5, segment.color = "grey50", segment.size = 0.2,
                  nudge_y = 0.05, size = 5, fontface = "italic") +
  labs(x = "Year of earliest indexed publication",
       y = "Total number of unique PMIDs (log scale)",
       size = "Ind. Obs.") +
  scale_y_log10() +
  scale_size_manual(values = c("< 10" = 2, "50 > x > 10" = 3,
                               "100 > x > 50" = 4, "> 100" = 6),
                    limits = limits,
                    labels = labels) +
  scale_color_manual(values = c("< 10" = "gray", "50 > x > 10" = "#DDCC77",
                                "100 > x > 50" = "#88CCEE", "> 100" = "#332288"),
                     limits = limits, labels = labels,
                     name = "Ind. Obs.") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x.bottom = element_text(size=20),
        axis.title.y = element_text(size=20),
        )


