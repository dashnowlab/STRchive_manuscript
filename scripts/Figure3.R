evidence <- data.frame(
  gene = c("AFF2", "AFF3", "AR", "ARX_1", "ARX_2", "ATN1", "ATXN1", "ATXN10", "ATXN2", "ATXN3",
           "ATXN7", "ATXN8OS", "BEAN1", "C9ORF72", "CACNA1A", "CBL", "COMP", "DAB1", "DIP2B",
           "DMD", "DMPK", "FMR1", "FOXL2", "FXN", "GIPC1", "GLS", "HOXA13_1", "HOXA13_2",
           "HOXA13_3", "HOXD13", "HTT", "JPH3", "LRP12", "MARCHF6", "NOP56", "NIPA1", "NOTCH2NLC",
           "NUTM2B-AS1", "PABPN1", "PHOX2B", "POLG", "PPP2R2B", "PRDM12", "RAPGEF2", "RFC1",
           "RILPL1", "RUNX2", "SAMD12", "SOX3", "STARD7", "TBP", "TBX1", "TCF4", "TNRC6A",
           "XYLT1", "YEATS2", "ZIC2", "ZIC3", "ZNF713", "CNBP", "CSTB", "EIF4A3", "PRNP",
           "VWA1", "ABCD3", "FGF14", "ZFHX3", "THAP11"
  ),
  ind_obs = c("100 > x > 50","< 10","> 100","50 > x > 10","< 10","100 > x > 50",
              "> 100","100 > x > 50","> 100","> 100","> 100","> 100","> 100",
              "> 100","> 100","< 10","< 10","50 > x > 10","< 10","< 10","> 100","> 100",
              "50 > x > 10","> 100","50 > x > 10","< 10","< 10","< 10","< 10","< 10","> 100",
              "100 > x > 50","100 > x > 50","> 100","100 > x > 50","> 100","> 100","< 10","> 100",
              "> 100","> 100","> 100","< 10","< 10","> 100","50 > x > 10","< 10","> 100","< 10",
              "50 > x > 10","100 > x > 50","< 10","> 100","< 10","50 > x > 10","< 10","50 > x > 10",
              "< 10","< 10","> 100","> 100","50 > x > 10","> 100","50 > x > 10","< 10","> 100",
              "50 > x > 10","< 10"))

STR_table_evidence <- merge(STR_table_clean, evidence,
                            by = "gene", all.x = TRUE)

labels <- c("< 10" = "< 10", "50 > x > 10" = "10 < x < 50",
            "100 > x > 50" = "50 < x < 100", "> 100" = "> 100")

limits <- c("< 10", "50 > x > 10", "100 > x > 50", "> 100")

ggplot(subset(STR_table_evidence, !is.na(age_onset_min) & Inheritance != ''),
       aes(color = ind_obs, x = reorder(gene, -age_onset_max))) +
  scale_x_discrete(name = 'Gene') +
  geom_linerange(aes(ymin = age_onset_min, ymax = age_onset_max), size = 3, alpha = 0.5) +
  geom_linerange(aes(ymin = typ_age_onset_min, ymax = typ_age_onset_max), size = 3) +
  geom_point(data = subset(STR_table_evidence, age_onset_min == age_onset_max), aes(y = age_onset_min), size = 2.5,
             alpha = 1, shape = 17) +
  geom_segment(aes(x = 1, y = 18, xend = 67, yend = 18), linetype = 'longdash', color = 'lightgrey') +
  coord_flip() + scale_y_continuous(breaks=c(0,18, 25, 50, 75, 100)) +
  labs(title = "Ranges in Age of Onset With Evidence Levels", y = "Age of Onset (years)") +
  scale_color_manual(values = c("< 10" = "darkgray", "50 > x > 10" = "lightblue",
                                "100 > x > 50" = "#3182BD", "> 100" = "darkblue"),
                     limits = limits, labels = labels,
                     name = "Ind. Obs.") +  # Rename legend for color (ind_obs)
  theme_minimal()
