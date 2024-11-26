# To plot all variation even without prevalence
ggplot(data = combined_df, aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "#E66100", size = 3, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 5, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percentage") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  scale_y_continuous(sec.axis = dup_axis()) + theme_minimal()

# for CNG only
ggplot(data = filter(combined_df, combined_df$Motif == "CNG"), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "#E66100", size = 5, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percentage") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  scale_y_continuous(sec.axis = dup_axis()) + theme_minimal()

ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "#E66100", size = 5, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percent of individuals") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  theme_minimal()

ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "#E66100", size = 5, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percent of individuals") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  ylim(0,0.64) +
  theme_minimal()


### plots overlaid with legend in PowerPoint File available in Data

### to get large icons for legend
ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)),
       aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.5, size = 3, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "#E66100", size = 16, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "#5D3A9B", size = 27, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percentage") +
  coord_flip() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.64) +
  theme_void()


