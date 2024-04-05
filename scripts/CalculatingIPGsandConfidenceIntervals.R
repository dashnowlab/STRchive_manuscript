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

# excluding CNG
total_prev <- subset(total, total$motif_norm != "CNG")

#while both inheritances are possible... they're more AD, so, we use that here
total_prev <- total_prev %>%
  mutate(Inheritance = ifelse(gene == 'ATXN2', 'AD', Inheritance))
total_prev <- total_prev %>%
  mutate(Inheritance = ifelse(gene == 'PABPN1', 'AD', Inheritance))

total_prev <- total_prev %>%
  mutate(repeatunit_path_normalized = ifelse(repeatunit_path_normalized == 'CNG', 'CAG,CCG,CGG,CTG', repeatunit_path_normalized))
total_prev$repeatunit_path_normalized_list <- strsplit(as.character(total_prev$repeatunit_path_normalized), ",")

AD_total <- subset(total_prev, Inheritance == "AD")
AD_Allele2_greater_than_path_min <- (AD_total$Allele2LowerBound >= AD_total$pathogenic_min) &
  (AD_total$motif_norm %in% AD_total$repeatunit_path_normalized_list)
AD_result1 <- aggregate(cbind(AD_Allele2_greater_than_path_min), by = list(AD_total$gene), FUN = sum, na.rm=TRUE)
AD_result2 <- aggregate(cbind(AD_total$Allele2), by = list(AD_total$gene, AD_total$repeatunit_path_normalized, AD_total$Inheritance), FUN = function(x) sum(!is.na(x)))
AD_result <- merge(AD_result1, AD_result2)

AR_total <- subset(total_prev, Inheritance == "AR")
# path result if we're NOT doing contraction for VWA1
AR_total <- AR_total %>%
  mutate(pathogenic_min = ifelse(gene == "VWA1", 3, pathogenic_min))
AR_Allele21_greater_than_path_min <- (AR_total$Allele2LowerBound >= AR_total$pathogenic_min) &
  (AR_total$motif_norm %in% AR_total$repeatunit_path_normalized_list) &
  (AR_total$motif_norm_small %in% AR_total$repeatunit_path_normalized_list) &
  (AR_total$Allele1LowerBound >= AR_total$pathogenic_min)
AR_Allele2_greater_than_path_min <- (AR_total$Allele2LowerBound >= AR_total$pathogenic_min) &
  (AR_total$Allele1LowerBound < AR_total$pathogenic_min) &
  (AR_total$motif_norm %in% AR_total$repeatunit_path_normalized)
AR_result1 <- aggregate(cbind(AR_Allele21_greater_than_path_min), by = list(AR_total$gene), FUN = sum, na.rm=TRUE)
### pathogenic result for VWA1, which is deviation, not contraction
#AR_result1$AR_Allele21_greater_than_path_min[AR_result1$Group.1 == "VWA1"] <- sum(subset(AR_total, gene == "VWA1" & Allele1 != 2 & Allele2 != 2 & (AR_total$motif_norm == AR_total$repeatunit_path_normalized))$gene == "VWA1")
AR_result2 <- aggregate(cbind(AR_Allele2_greater_than_path_min), by = list(AR_total$gene), FUN = sum, na.rm=TRUE)
### carrier result for VWA1, which is deviation, not contraction
#AR_result2$AR_Allele2_greater_than_path_min[AR_result2$Group.1 == "VWA1"] <- sum(subset(AR_total, gene == "VWA1" & (Allele1 == 2 & Allele2 != 2) | (Allele1 != 2 & Allele2 == 2) & (AR_total$motif_norm == AR_total$repeatunit_path_normalized))$gene == "VWA1")
AR_result3 <- aggregate(cbind(AR_total$Allele2), by = list(AR_total$gene, AR_total$repeatunit_path_normalized, AR_total$Inheritance), FUN = function(x) sum(!is.na(x)))
AR_result <- merge(AR_result1, AR_result2)
AR_result <- merge(AR_result, AR_result3)

XR_XX_total_sex <- subset(total_prev, Inheritance == "XR" & Sex == "XX")
XR_XX_Allele21_greater_than_path_min <- (XR_XX_total_sex$Allele1LowerBound >= XR_XX_total_sex$pathogenic_min) &
  (XR_XX_total_sex$Allele2LowerBound >= XR_XX_total_sex$pathogenic_min) &
  (XR_XX_total_sex$motif_norm %in% XR_XX_total_sex$repeatunit_path_normalized) &
  (XR_XX_total_sex$motif_norm_small %in% XR_XX_total_sex$repeatunit_path_normalized)
XR_XX_Allele2_greater_than_path_min <- (XR_XX_total_sex$Allele1LowerBound < XR_XX_total_sex$pathogenic_min) &
  (XR_XX_total_sex$Allele2LowerBound >= XR_XX_total_sex$pathogenic_min) &
  (XR_XX_total_sex$motif_norm %in% XR_XX_total_sex$repeatunit_path_normalized)
XR_XX_result1 <- aggregate(cbind(XR_XX_Allele21_greater_than_path_min), by = list(XR_XX_total_sex$gene), FUN = sum, na.rm=TRUE)
XR_XX_result2 <- aggregate(cbind(XR_XX_Allele2_greater_than_path_min), by = list(XR_XX_total_sex$gene), FUN = sum, na.rm=TRUE)
XR_XX_result3 <- aggregate(cbind(XR_XX_total_sex$Allele1), by = list(XR_XX_total_sex$gene, XR_XX_total_sex$repeatunit_path_normalized, XR_XX_total_sex$Inheritance), FUN = function(x) sum(!is.na(x)))
XR_XX_result <- merge(XR_XX_result1, XR_XX_result2)
XR_XX_result <- merge(XR_XX_result, XR_XX_result3)

XR_XY_total_sex <- subset(total_prev, Inheritance == "XR" & Sex == "XY")
XR_XY_Allele1_greater_than_path_min <- (XR_XY_total_sex$Allele1LowerBound >= XR_XY_total_sex$pathogenic_min) &
  (XR_XY_total_sex$motif_norm %in% XR_XY_total_sex$repeatunit_path_normalized_list)
XR_XY_Allele2_greater_than_path_min <- XR_XY_total_sex$Allele2LowerBound >= XR_XY_total_sex$pathogenic_min &
  (XR_XY_total_sex$motif_norm %in% XR_XY_total_sex$repeatunit_path_normalized_list)
XR_XY_Allele_greater_than_path_min <- XR_XY_Allele1_greater_than_path_min | XR_XY_Allele2_greater_than_path_min
XR_XY_result1 <- aggregate(cbind(XR_XY_Allele_greater_than_path_min), by = list(XR_XY_total_sex$gene), FUN = sum, na.rm=TRUE)
XR_XY_result2 <- aggregate(cbind(XR_XY_total_sex$Allele1), by = list(XR_XY_total_sex$gene, XR_XY_total_sex$repeatunit_path_normalized, XR_XY_total_sex$Inheritance), FUN = function(x) sum(!is.na(x)))
XR_XY_result <- merge(XR_XY_result1, XR_XY_result2)

XD_total_sex <- subset(total_prev, Inheritance == "XD")
XD_Allele_greater_than_path_min <- ((XD_total_sex$Allele1LowerBound >= XD_total_sex$pathogenic_min) &
                                      (XD_total_sex$motif_norm_small %in% XD_total_sex$repeatunit_path_normalized_list)) |
  ((XD_total_sex$Allele2LowerBound >= XD_total_sex$pathogenic_min)) &
  (XD_total_sex$motif_norm %in% XD_total_sex$repeatunit_path_normalized_list)
XD_result1 <- aggregate(cbind(XD_Allele_greater_than_path_min), by = list(XD_total_sex$gene), FUN = sum, na.rm=TRUE)
XD_result2 <- aggregate(cbind(XD_total_sex$Allele1), by = list(XD_total_sex$gene, XD_total_sex$repeatunit_path_normalized, XD_total_sex$Inheritance), FUN = function(x) sum(!is.na(x)))
XD_result <- merge(XD_result1, XD_result2)


colnames(AD_result) <- c("gene", "Pathogenic_Count", "Motif", "Inheritance", "Total_Loci")
colnames(AR_result) <- c("gene", "Pathogenic_Count", "Carrier_Count", "Motif", "Inheritance", "Total_Loci")
colnames(XR_XX_result) <- c("gene", "Pathogenic_Count", "Carrier_Count", "Motif", "Inheritance", "Total_Loci")
colnames(XR_XY_result) <- c("gene", "Pathogenic_Count", "Motif", "Inheritance", "Total_Loci")
colnames(XD_result) <- c("gene", "Pathogenic_Count", "Motif", "Inheritance", "Total_Loci")

XR_XY_result$gene <- paste0("XY_", XR_XY_result$gene)
XR_XX_result$gene <- paste0("XX_", XR_XX_result$gene)

AD_result <- transform(AD_result, Carrier_Count = NA)
XR_XY_result <- transform(XR_XY_result, Carrier_Count = NA)
XD_result <- transform(XD_result, Carrier_Count = NA)

# combine data frames
combined_df <- do.call(rbind, list(AD_result, AR_result, XR_XX_result, XR_XY_result, XD_result))

combined_df$pathogenic_percent <- combined_df$Pathogenic_Count/combined_df$Total_Loci*100

combined_df$carrier_percent <- combined_df$Carrier_Count/combined_df$Total_Loci*100

combined_df <- calculate_ci(combined_df, successes_col = "Pathogenic_Count", trials_col = "Total_Loci", conf_level = 0.95)
combined_df <- calculate_ci(combined_df, successes_col = "Carrier_Count", trials_col = "Total_Loci", conf_level = 0.95)

combined_df <- left_join(combined_df, STR_table %>% dplyr::select(gene, Prevalence), by = "gene")


#where we don't know or have sex-specific prevalence
combined_df["Prevalence"][combined_df["Prevalence"] == ''] <- NA
combined_df$Prevalence[combined_df$gene == "XY_AFF2"] <- 2/50000
combined_df$Prevalence[combined_df$gene == "XY_AR"] <-  1/30000
combined_df$Prevalence[combined_df$gene == "XY_DMD"] <- 4.8/100000

combined_df$prevalence_dec <- sapply(combined_df$Prevalence, function(fraction) eval(parse(text = fraction)))
