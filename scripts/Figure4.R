# To plot all variation even without prevalence
# ggplot(data = combined_df, aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
#   geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
#   geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 3, alpha = 0.8, shape = 23) +
#   geom_point(aes(y = pathogenic_percent), color = "red", size = 5, alpha = 0.9, shape = 20) +
#   labs(x = "Gene", y = "Percentage") +
#   theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
#         plot.title = element_text(hjust = 0.5)) +
#   coord_flip() +
#   scale_y_continuous(sec.axis = dup_axis())



### If we exclude CNG
color_labels <- c('black', 'black', 'black', 'black', 'black',
                  'black', 'black', 'black', 'black','black', 'black', 'black', 'black', 'black', 'black',
                  'black', 'black', 'black')

### If we include CNG
# color_labels <- c('black', 'black', 'black', 'black', 'black',
#                   'black', 'black', 'black', 'blue','black', 'black', 'black', 'black',
#                   'black', 'black', 'black', 'black',
#                   'purple', 'lightpink', 'orange')
#
#
# ggplot(data = filter(combined_df, combined_df$Motif == "CNG"), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
#   geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
#   geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 5, alpha = 0.8, shape = 23) +
#   geom_point(aes(y = pathogenic_percent), color = "red", size = 7, alpha = 0.9, shape = 20) +
#   labs(x = "Gene", y = "Percentage") +
#   theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
#         plot.title = element_text(hjust = 0.5)) +
#   coord_flip() +
#   scale_y_continuous(sec.axis = dup_axis())

ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 5, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "red", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percentage") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  theme(axis.text.y = element_text(colour = color_labels)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  theme_minimal()



ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)), aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
  geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.3, size = 2, color = "black") +
  geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 5, alpha = 0.8, shape = 23) +
  geom_point(aes(y = pathogenic_percent), color = "red", size = 7, alpha = 0.9, shape = 20) +
  labs(x = "Gene", y = "Percentage") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
        plot.title = element_text(hjust = 0.5)) +
  coord_flip() +
  theme(axis.text.y = element_text(colour = color_labels)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  ylim(0,0.64) +
  theme(axis.text.y = element_text(colour = color_labels)) +
  theme_minimal()


### plots overlaid with legend in PowerPoint File available in Data

### to get large icons for legend
# ggplot(data = filter(combined_df, !is.na(combined_df$prevalence)),
#        aes(x = reorder(gene, pathogenic_percent), y = pathogenic_percent)) +
#   geom_errorbar(aes(ymin = Pathogenic_Count_lower_ci, ymax = Pathogenic_Count_upper_ci), width = 0.5, size = 3, color = "black") +
#   geom_point(aes(y = prevalence_dec*100), fill = "blue", size = 15, alpha = 0.8, shape = 23) +
#   geom_point(aes(y = pathogenic_percent), color = "red", size = 30, alpha = 0.9, shape = 20) +
#   labs(x = "Gene", y = "Percentage") +
#   coord_flip() +
#   theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size =15),
#         plot.title = element_text(hjust = 0.5)) +
#   ylim(0,0.64) +
#   theme(axis.text.y = element_text(colour = color_labels))
#

