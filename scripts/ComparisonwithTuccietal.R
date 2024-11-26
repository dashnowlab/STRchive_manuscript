### Run ProcessingSTRchiveandgnomADData.R
### Run CalculatingPGsandConfidenceIntervals.R

### Tucci data found in supplementarytables20230623.xlsx from
### https://www.researchsquare.com/article/rs-3097805/v1

### Pathogenic cutoffs found in Table S3

### Results found in Table S5 and post-QC data below
### C9orf72 intermediate dashes replaced with 0
### AR males changed to XY_AR


Tucci <- data.frame(
  stringsAsFactors = FALSE,
       check.names = FALSE,
            Cohort = c("100kGP","TOPMed","100kGP",
                       "TOPMed","100kGP","TOPMed","100kGP","TOPMed","100kGP",
                       "TOPMed","100kGP","TOPMed","100kGP","TOPMed",
                       "100kGP","TOPMed","100kGP","TOPMed","100kGP","TOPMed",
                       "100kGP","TOPMed","100kGP","TOPMed","100kGP","TOPMed"),
             gene = c("XY_AR","XY_AR","ATN1",
                       "ATN1","ATXN1","ATXN1","ATXN2","ATXN2","ATXN3",
                       "ATXN3","ATXN7","ATXN7","C9ORF72","C9ORF72","CACNA1A",
                       "CACNA1A","DMPK","DMPK","FXN","FXN","HTT","HTT",
                       "JPH3","JPH3","TBP","TBP"),
            Normal = c(15449L,17996L,34186L,47964L,
                       34156L,47945L,34140L,47920L,34190L,47980L,34189L,
                       47982L,34151L,47952L,34167L,47951L,34024L,47790L,
                       34163L,47544L,34120L,47921L,34185L,47923L,34158L,
                       47944L),
      Intermediate = c(18,12,4,21,29,
                       27,34,36,0,6,0,2,0,0,15,27,
                       144,172,27,74,57,58,5,62,32,
                       37),
   Full_mutation = c(2L,10L,0L,1L,5L,14L,16L,
                       30L,0L,0L,1L,2L,39L,34L,8L,8L,22L,24L,455L,
                       368L,13L,7L,0L,1L,0L,5L)
)

### Calculate PG frequency for cohort
Tucci$percentage<- (Tucci$Full_mutation/(Tucci$Intermediate+Tucci$Normal)*100)

Tucci <- subset(Tucci, select = -c(Normal,Intermediate, Full_mutation))

### Separate by cohort
Tucci_wide <- Tucci %>%
  pivot_wider(names_from = Cohort, values_from = percentage, names_prefix = "percentage_")

combined_df <- subset(combined_df, select = c(gene, pathogenic_percent,
                                                    Pathogenic_Count_lower_ci,
                                                    Pathogenic_Count_upper_ci,
                                              carrier_percent,
                                              Carrier_Count_lower_ci,
                                              Carrier_Count_upper_ci,
                                              inheritance))

### Tucci is really calculating a carrier frequency; FXN is the only AR in this
### data set so we need to move the carrier to pathogenic for meaningful comparison
combined_df$pathogenic_percent[combined_df$gene == "FXN"] <- combined_df$carrier_percent[combined_df$gene == "FXN"]
combined_df$Pathogenic_Count_lower_ci[combined_df$gene == "FXN"] <- combined_df$Carrier_Count_lower_ci[combined_df$gene == "FXN"]
combined_df$Pathogenic_Count_upper_ci[combined_df$gene == "FXN"] <- combined_df$Carrier_Count_upper_ci[combined_df$gene == "FXN"]

### combined Tucci and gnomAD data
combined_df_Tucci <- merge(Tucci_wide, combined_df, by="gene")

### Adjust cutoff for ATXN1, C9ORF72 and ATXN7 for gnomAD analysis
### ATXN1: 44 (we used 39)
### ATXN7: 36 (we used 34)
### C9ORF72, we used 250 and they used 31!
### all other loci had pathogenic thresholds

AD_total <- subset(AD_total, gene %in% c("ATXN1", "ATXN7", "C9ORF72", "ATXN2"))
AD_total$pathogenic_min[AD_total$gene == "ATXN1"] <- 44
AD_total$pathogenic_min[AD_total$gene == "ATXN2"] <- 33
AD_total$pathogenic_min[AD_total$gene == "ATXN7"] <- 36
AD_total$pathogenic_min[AD_total$gene == "C9ORF72"] <- 31

### get new PG estimates with confidence interval
AD_total <- AD_total %>%
  mutate(
    AD_Allele2_greater_than_path_min = (Allele2LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized_list)
  )

AD_result1 <- AD_total %>%
  group_by(gene) %>%
  summarise(
    Allele2_greater_than_path_min_sum = sum(AD_Allele2_greater_than_path_min, na.rm = TRUE),
    .groups = "drop"
  )

AD_result2 <- AD_total %>%
  group_by(gene, repeatunit_path_normalized, inheritance) %>%
  summarise(
    Allele2_non_na_count = sum(!is.na(Allele2)),
    .groups = "drop"
  )

adj_gnomAD_df <- merge(AD_result1, AD_result2)
colnames(adj_gnomAD_df) <- c("gene", "Pathogenic_Count", "Motif", "inheritance", "Total_Loci")

adj_gnomAD_df$pathogenic_percent <- adj_gnomAD_df$Pathogenic_Count/adj_gnomAD_df$Total_Loci*100

adj_gnomAD_df <- calculate_ci(adj_gnomAD_df, successes_col = "Pathogenic_Count",
                              trials_col = "Total_Loci", conf_level = 0.95)

adj_gnomAD_df$gene[adj_gnomAD_df$gene == "ATXN1"] <- "ATXN1_adj"
adj_gnomAD_df$gene[adj_gnomAD_df$gene == "ATXN7"] <- "ATXN7_adj"
adj_gnomAD_df$gene[adj_gnomAD_df$gene == "ATXN2"] <- "ATXN2_adj"
adj_gnomAD_df$gene[adj_gnomAD_df$gene == "C9ORF72"] <- "C9ORF72_adj"


adj_gnomAD_df <- subset(adj_gnomAD_df, select = -c(Motif, Total_Loci,
                                               Pathogenic_Count))

combined_df <- subset(combined_df, select = -c(carrier_percent, Carrier_Count_lower_ci,
                                               Carrier_Count_upper_ci))

combined_df_adj <- do.call(rbind, list(combined_df, adj_gnomAD_df))

### Make sure our merged result keeps the adj genes
Tucci_wide_duplicates <- Tucci_wide %>%
  filter(gene %in% c("C9ORF72", "ATXN1", "ATXN7", "ATXN2")) %>%
  mutate(gene = paste0(gene, "_adj"))

Tucci_wide <- bind_rows(Tucci_wide, Tucci_wide_duplicates)


combined_df_Tucci_adj <- merge(Tucci_wide, combined_df_adj, by="gene")

### In range?
combined_df_Tucci_adj <- combined_df_Tucci_adj %>%
  mutate(kgp_inrange = if_else(percentage_100kGP >= Pathogenic_Count_lower_ci &
                                 percentage_100kGP <= Pathogenic_Count_upper_ci, "YES", "NO"),
         TOPMed_inrange = if_else(percentage_TOPMed >= Pathogenic_Count_lower_ci &
                                    percentage_TOPMed <= Pathogenic_Count_upper_ci, "YES", "NO"))

### Plot All Together
### Order by pathogenic count, with adjusted genes adjacent
combined_df_Tucci_adj$gene <- factor(combined_df_Tucci_adj$gene, levels=c("ATN1",
                "ATXN3", "C9ORF72_adj","C9ORF72", "TBP", "CACNA1A", "JPH3", "HTT",
                "ATXN7_adj", "ATXN7", "DMPK", "XY_AR", "ATXN2_adj", "ATXN2", "ATXN1_adj","ATXN1",
                "FXN"))

ggplot(data = filter(combined_df_Tucci_adj), aes(x = gene, y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = percentage_TOPMed), fill = "#FEFE62", size = 5, alpha = 0.8, shape = 24) +
  geom_point(aes(y = percentage_100kGP), fill = "#40B0A6", size = 5, alpha = 0.8, shape = 25) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percent of individuals") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  theme_minimal()


ggplot(data = filter(combined_df_Tucci_adj), aes(x = gene, y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = percentage_TOPMed), fill = "#FEFE62", size = 5, alpha = 0.8, shape = 24) +
  geom_point(aes(y = percentage_100kGP), fill = "#40B0A6", size = 5, alpha = 0.8, shape = 25) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percent of individuals") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.27) +
  coord_flip() +
  theme_minimal()


### for making the image, which is in GoogleSlide
# ggplot(data = filter(combined_df_Tucci_adj), aes(x = gene, y = pathogenic_percent)) +
#   geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
#   geom_point(aes(y = percentage_TOPMed), fill = "#FEFE62", size = 15, alpha = 0.8, shape = 24) +
#   geom_point(aes(y = percentage_100kGP), fill = "#40B0A6", size =15, alpha = 0.8, shape = 25) +
#   geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 7, alpha = 0.9, shape = 20) +
#   labs(x = "Gene", y = "Percent of individuals") +
#   theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
#         plot.title = element_text(hjust = 0.5)) +
#   coord_flip() +
#   theme_void()
