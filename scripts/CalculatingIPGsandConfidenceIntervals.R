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
