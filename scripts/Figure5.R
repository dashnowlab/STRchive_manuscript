color_labels <- c('black', 'black', 'black', 'black', 'black',
                  'black', 'black', 'black', 'blue','black', 'black', 'black', 'black', 'black', 'black',
                  'purple', 'lightpink', 'orange')

ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 3, alpha = 0.9, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "red", size = 3, alpha = 0.9, shape = 20) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.2, color = "black") +
  labs(x = "Gene", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =8),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  theme(axis.text.x = element_text(colour = color_labels))

ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)),
       aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 3, alpha = 0.9, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "red", size = 3, alpha = 0.9, shape = 20) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.2, color = "black") +
  labs(x = "Gene", y = "Percentage") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =8),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.64) +
  theme(axis.text.x = element_text(colour = color_labels))
