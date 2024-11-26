### Functions
# Calculating 95% confidence interval
calculate_ci <- function(df, successes_col, trials_col, conf_level, method) {

  # Calculate the estimated proportion and confidence intervals for each row
  est_props <- df[[successes_col]] / df[[trials_col]]
  cis <- BinomCI(df[[successes_col]], df[[trials_col]], conf.level = conf_level, method = method)

  # Extract the lower and upper confidence interval bounds
  lower_cis <- cis[, "lwr.ci"]
  upper_cis <- cis[, "upr.ci"]

  # Add the confidence intervals to the original data frame with prefixed names
  df[[paste0(successes_col, "_lower_ci")]] <- lower_cis*100
  df[[paste0(successes_col, "_upper_ci")]] <- upper_cis*100

  # Return the updated data frame
  return(df)
}

# Checking lowerbound estimates against allele genotypes
# table(as.numeric(total$Allele2LowerBound) == as.numeric(total$Allele2), useNA = "ifany")
#
# table(as.numeric(total$Allele1LowerBound) == as.numeric(total$Allele1), useNA = "ifany")
#
# 60201+965298
#
# 965298/1025499
#
# 33098+1076598
#
# 1076598/1109696
#
# ### Allele 1
# differences <- abs(total$Allele1LowerBound - total$Allele1)
# average_difference <- mean(differences)
# differences <- total[total$Allele1LowerBound != total$Allele1, c("Allele1LowerBound", "Allele1")]
# differences <- abs(differences$Allele1LowerBound - differences$Allele1)
# average_difference <- mean(differences)
# range(differences)
# median(differences)
#
# ### Allele 2
# differences <- abs(total$Allele2LowerBound - total$Allele2)
# average_difference <- mean(differences, na.rm = TRUE)
#
# differences <- total[total$Allele2LowerBound != total$Allele2, c("Allele2LowerBound", "Allele2")]
#
# differences <- abs(differences$Allele2LowerBound - differences$Allele2)
# average_difference <- mean(differences, na.rm = TRUE)
# range(differences, na.rm = TRUE)
# median(differences, na.rm = TRUE)

# excluding CNG
total_prev <- subset(total, total$motif_norm != "CNG")

#including CNG
#total_prev <- total

#while both inheritances are possible... they're more AD, so, we use that here
total_prev <- total_prev %>%
  mutate(inheritance = ifelse(gene == 'ATXN2', 'AD', inheritance))
total_prev <- total_prev %>%
  mutate(inheritance = ifelse(gene == 'PABPN1', 'AD', inheritance))

# FOXL2 is AD with one allele >23, and AR if both alleles > 19
# since all allele1_lowerbound is 14, we'll use AD with a 24 path min
total_prev <- total_prev %>%
  mutate(inheritance = ifelse(gene == 'FOXL2', 'AD', inheritance))

total_prev <- total_prev %>%
  mutate(pathogenic_min = ifelse(gene == 'FOXL2',24, pathogenic_min))


total_prev$repeatunit_path_normalized_list <- strsplit(as.character(total_prev$repeatunit_path_normalized), ",")

AD_total <- subset(total_prev, inheritance == "AD")
# Step 1: Calculate AD_Allele2_greater_than_path_min
AD_total <- AD_total %>%
  mutate(
    AD_Allele2_greater_than_path_min = (Allele2LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized_list)
  )

# Step 2: Aggregate AD_Allele2_greater_than_path_min by gene
AD_result1 <- AD_total %>%
  group_by(gene) %>%
  summarise(
    Allele2_greater_than_path_min_sum = sum(AD_Allele2_greater_than_path_min, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Aggregate Allele2 by gene, repeatunit_path_normalized
AD_result2 <- AD_total %>%
  group_by(gene, repeatunit_path_normalized, inheritance) %>%
  summarise(
    Allele2_non_na_count = sum(!is.na(Allele2)),
    .groups = "drop"
  )

# Step 4: Merge the two results
AD_result <- AD_result1 %>%
  left_join(AD_result2, by = "gene")


AR_total <- subset(total_prev, inheritance == "AR")
# path result if we're NOT doing contraction for VWA1
# Step 1: Update pathogenic_min for "VWA1"
AR_total <- AR_total %>%
  mutate(pathogenic_min = ifelse(gene == "VWA1", 3, pathogenic_min))

# Step 2: Calculate logical conditions as new columns
AR_total <- AR_total %>%
  mutate(
    AR_Allele21_greater_than_path_min = (Allele2LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized_list) &
      (motif_norm_small %in% repeatunit_path_normalized_list) &
      (Allele1LowerBound >= pathogenic_min),

    AR_Allele2_greater_than_path_min = (Allele2LowerBound >= pathogenic_min) &
      (Allele1LowerBound < pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized)
  )

# Step 3: Aggregate results for AR_Allele21_greater_than_path_min and AR_Allele2_greater_than_path_min
AR_result1 <- AR_total %>%
  group_by(gene, inheritance) %>%
  summarise(
    AR_Allele21_sum = sum(AR_Allele21_greater_than_path_min, na.rm = TRUE),
    AR_Allele2_sum = sum(AR_Allele2_greater_than_path_min, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Aggregate Allele2 counts by gene, repeatunit_path_normalized, and inheritance
AR_result3 <- AR_total %>%
  group_by(gene, repeatunit_path_normalized) %>%
  summarise(
    Allele2_count = sum(!is.na(Allele2)),
    .groups = "drop"
  )

# Step 5: Merge the results
AR_result <- AR_result1 %>%
  left_join(AR_result3, by = "gene")

# XR_XX_total_sex
XR_XX_result <- total_prev %>%
  filter(inheritance == "XR", Sex == "XX") %>%
  mutate(
    Allele21_greater_than_path_min = (Allele1LowerBound >= pathogenic_min) &
      (Allele2LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized) &
      (motif_norm_small %in% repeatunit_path_normalized),

    Allele2_greater_than_path_min = (Allele1LowerBound < pathogenic_min) &
      (Allele2LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized)
  ) %>%
  group_by(gene, inheritance, repeatunit_path_normalized_list) %>%
  summarise(
    Allele21_sum = sum(Allele21_greater_than_path_min, na.rm = TRUE),
    Allele2_sum = sum(Allele2_greater_than_path_min, na.rm = TRUE),
    Allele1_count = sum(!is.na(Allele1)),
    .groups = "drop"
  )

# XR_XY_total_sex
XR_XY_result <- total_prev %>%
  filter(inheritance == "XR", Sex == "XY") %>%
  mutate(
    Allele1_greater_than_path_min = (Allele1LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized_list),

    Allele2_greater_than_path_min = (Allele2LowerBound >= pathogenic_min) &
      (motif_norm %in% repeatunit_path_normalized_list),

    Allele_greater_than_path_min = Allele1_greater_than_path_min | Allele2_greater_than_path_min
  ) %>%
  group_by(gene, inheritance, repeatunit_path_normalized_list) %>%
  summarise(
    Allele_greater_sum = sum(Allele_greater_than_path_min, na.rm = TRUE),
    Allele1_count = sum(!is.na(Allele1)),
    .groups = "drop"
  )

# XD_total_sex
XD_result <- total_prev %>%
  filter(inheritance == "XD") %>%
  mutate(
    Allele_greater_than_path_min = (
      ((Allele1LowerBound >= pathogenic_min) &
         (motif_norm_small %in% repeatunit_path_normalized_list)) |
        ((Allele2LowerBound >= pathogenic_min) &
           (motif_norm %in% repeatunit_path_normalized_list))
    )
  ) %>%
  group_by(gene, inheritance, repeatunit_path_normalized) %>%
  summarise(
    Allele_greater_sum = sum(Allele_greater_than_path_min, na.rm = TRUE),
    Allele1_count = sum(!is.na(Allele1)),
    .groups = "drop"
  )

colnames(AD_result) <- c("gene", "Pathogenic_Count", "Motif", "inheritance", "Total_Loci")
colnames(AR_result) <- c("gene", "inheritance", "Pathogenic_Count", "Carrier_Count", "Motif",  "Total_Loci")
colnames(XR_XX_result) <- c("gene", "inheritance", "Motif", "Pathogenic_Count", "Carrier_Count",  "Total_Loci")
colnames(XR_XY_result) <- c("gene", "inheritance", "Motif", "Pathogenic_Count", "Total_Loci")
colnames(XD_result) <- c("gene", "inheritance", "Motif", "Pathogenic_Count",  "Total_Loci")

XR_XY_result$gene <- paste0("XY_", XR_XY_result$gene)
XR_XX_result$gene <- paste0("XX_", XR_XX_result$gene)

AD_result <- transform(AD_result, Carrier_Count = NA)
XR_XY_result <- transform(XR_XY_result, Carrier_Count = NA)
XD_result <- transform(XD_result, Carrier_Count = NA)

# combine data frames
combined_df <- do.call(rbind, list(AD_result, AR_result, XR_XX_result, XR_XY_result, XD_result))

combined_df$pathogenic_percent <- (combined_df$Pathogenic_Count*100)/(combined_df$Total_Loci)

combined_df$carrier_percent <- combined_df$Carrier_Count/combined_df$Total_Loci*100

combined_df <- calculate_ci(combined_df, successes_col = "Pathogenic_Count", trials_col = "Total_Loci", conf_level = 0.95)
combined_df <- calculate_ci(combined_df, successes_col = "Carrier_Count", trials_col = "Total_Loci", conf_level = 0.95)

combined_df <- left_join(combined_df, STR_table %>% dplyr::select(gene, prevalence), by = "gene")


#where we don't know or have sex-specific prevalence
combined_df["prevalence"][combined_df["prevalence"] == ''] <- NA
combined_df$prevalence[combined_df$gene == "XY_AFF2"] <- 2/50000
combined_df$prevalence[combined_df$gene == "XY_AR"] <-  1/30000
combined_df$prevalence[combined_df$gene == "XY_DMD"] <- 4.8/100000

combined_df$prevalence_dec <- sapply(combined_df$prevalence, function(fraction) eval(parse(text = fraction)))

# Remove the loci that are suspicious
manual_suspects <- c("XY_AFF2", "ARX_1", "NOTCH2NLC", "TBP", "ZNF713", "ARX_2", "FOXL2", "HOXA13_1", "HOXA13_2", "HOXA13_3", "NIPA1", "PABPN1", "PHOX2B", "RUNX2", "SOX3", "TBX1", "ZIC2")

combined_df <- combined_df %>%
  filter(!gene %in% manual_suspects)
