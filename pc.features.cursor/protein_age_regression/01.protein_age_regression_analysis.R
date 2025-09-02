#!/usr/bin/env Rscript

# Protein-Age Regression Analysis
# This script performs linear regression of age on each protein's NPQ expression
# in the Methuselah dataset, adjusts p-values with FDR, and creates a volcano plot

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(ggrepel)

# Set random seed for reproducibility
set.seed(42)

cat("Starting Protein-Age Regression Analysis...\n")

# Load Methuselah data
cat("Loading Methuselah data...\n")

# Load inflammation data
alpha.data.inflam <- read_csv("data/methuselah.alpha/inflammation/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
  filter(!grepl("SC", SampleType)) %>%
  mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))

# Load CNS data
alpha.data.cns <- read_csv("data/methuselah.alpha/CNS/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
  filter(!grepl("SC", SampleType)) %>%
  mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))

# Combine inflammation and CNS data
org.alpha.data = rbind(alpha.data.inflam, alpha.data.cns) %>% 
  group_by(SampleName, subject.id, Target, SampleType) %>%
  summarise(NPQ = mean(NPQ), .groups = 'drop')

# Load metadata
metadata = read_csv('data/methuselah.alpha/sample_metadata.csv')

# Merge data with metadata
full.data = org.alpha.data %>% 
  left_join(metadata, by = c('SampleName', 'SampleType'))

# Create wide format data
full.data.wide = full.data %>% 
  select(SampleName, Sex, Age2, Target, NPQ) %>%
  pivot_wider(names_from = 'Target', values_from = 'NPQ') %>%
  mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$")) %>% 
  mutate(SampleName = subject.id) %>%
  select(-subject.id)

# Convert Sex to numeric (M=1, F=0)
full.data.wide$Sex = ifelse(full.data.wide$Sex == 'M', 1, 0)

cat("Data loaded successfully.\n")
cat("Number of samples:", nrow(full.data.wide), "\n")
cat("Number of proteins:", ncol(full.data.wide) - 3, "\n") # -3 for SampleName, Sex, Age2

# Function to perform regression for a single protein
perform_protein_regression <- function(protein_data, age, sex) {
  # Create data frame for regression
  reg_data <- data.frame(
    age = age,
    sex = sex,
    protein = protein_data
  )
  
  # Remove rows with missing values
  reg_data <- reg_data[complete.cases(reg_data), ]
  
  if(nrow(reg_data) < 10) {
    return(list(
      estimate = NA,
      p_value = NA,
      r_squared = NA,
      n_samples = nrow(reg_data)
    ))
  }
  
  # Perform linear regression: age ~ protein + sex
  tryCatch({
    model <- lm(age ~ protein + sex, data = reg_data)
    summary_model <- summary(model)
    
    # Extract protein coefficient (not intercept or sex)
    protein_coef <- coef(model)["protein"]
    protein_p <- summary_model$coefficients["protein", "Pr(>|t|)"]
    r_sq <- summary_model$r.squared
    
    return(list(
      estimate = protein_coef,
      p_value = protein_p,
      r_squared = r_sq,
      n_samples = nrow(reg_data)
    ))
  }, error = function(e) {
    return(list(
      estimate = NA,
      p_value = NA,
      r_squared = NA,
      n_samples = nrow(reg_data)
    ))
  })
}

# Function to calculate fold change (using median values)
calculate_fold_change <- function(protein_data, age, age_threshold = 65) {
  # Split data into young and old groups
  young_data <- protein_data[age < age_threshold]
  old_data <- protein_data[age >= age_threshold]
  
  # Remove NA values
  young_data <- young_data[!is.na(young_data)]
  old_data <- old_data[!is.na(old_data)]
  
  if(length(young_data) < 5 || length(old_data) < 5) {
    return(NA)
  }
  
  # Calculate median values
  young_median <- median(young_data, na.rm = TRUE)
  old_median <- median(old_data, na.rm = TRUE)
  
  # Calculate fold change (old/young)
  if(young_median == 0) {
    return(NA)
  }
  
  fold_change <- old_median / young_median
  
  # Convert to log2 fold change
  log2_fold_change <- log2(fold_change)
  
  return(log2_fold_change)
}

# Extract age and sex data
age_data <- full.data.wide$Age2
sex_data <- full.data.wide$Sex

# Get protein columns (exclude SampleName, Sex, Age2)
protein_cols <- colnames(full.data.wide)[!colnames(full.data.wide) %in% c("SampleName", "Sex", "Age2")]

cat("Performing regression analysis for", length(protein_cols), "proteins...\n")

# Initialize results storage
regression_results <- list()
fold_change_results <- numeric(length(protein_cols))

# Process each protein
for(i in 1:length(protein_cols)) {
  protein_name <- protein_cols[i]
  protein_data <- full.data.wide[[protein_name]]
  
  if(i %% 50 == 0) {
    cat("Processed", i, "of", length(protein_cols), "proteins...\n")
  }
  
  # Perform regression
  reg_result <- perform_protein_regression(protein_data, age_data, sex_data)
  
  # Calculate fold change
  fc_result <- calculate_fold_change(protein_data, age_data)
  
  # Store results
  regression_results[[i]] <- reg_result
  fold_change_results[i] <- fc_result
}

cat("Regression analysis completed.\n")

# Convert results to data frame
results_df <- data.frame(
  protein = protein_cols,
  estimate = sapply(regression_results, function(x) x$estimate),
  p_value = sapply(regression_results, function(x) x$p_value),
  r_squared = sapply(regression_results, function(x) x$r_squared),
  n_samples = sapply(regression_results, function(x) x$n_samples),
  log2_fold_change = fold_change_results
)

# Remove rows with missing values
results_df <- results_df[complete.cases(results_df), ]

# Adjust p-values using FDR (Benjamini-Hochberg)
results_df$adj_p_value <- p.adjust(results_df$p_value, method = "BH")

# Add significance levels
results_df$significance <- case_when(
  results_df$adj_p_value < 0.001 ~ "***",
  results_df$adj_p_value < 0.01 ~ "**",
  results_df$adj_p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Add significance categories for plotting
results_df$significance_category <- case_when(
  results_df$adj_p_value < 0.05 & abs(results_df$log2_fold_change) > 0.5 ~ "Significant & High FC",
  results_df$adj_p_value < 0.05 ~ "Significant",
  abs(results_df$log2_fold_change) > 0.5 ~ "High FC",
  TRUE ~ "Not Significant"
)

cat("Results summary:\n")
cat("Total proteins analyzed:", nrow(results_df), "\n")
cat("Significant proteins (adj p < 0.05):", sum(results_df$adj_p_value < 0.05), "\n")
cat("High fold change proteins (|log2FC| > 0.5):", sum(abs(results_df$log2_fold_change) > 0.5), "\n")
cat("Significant & high fold change:", sum(results_df$adj_p_value < 0.05 & abs(results_df$log2_fold_change) > 0.5), "\n")

# Save results
write.csv(results_df, "protein_age_regression/protein_age_regression_results.csv", row.names = FALSE)
cat("Results saved to: protein_age_regression/protein_age_regression_results.csv\n")

# Create volcano plot
cat("Creating volcano plot...\n")

# Define colors for significance categories
color_palette <- c(
  "Significant & High FC" = "#D7261E",
  "Significant" = "#F46036", 
  "High FC" = "#2E294E",
  "Not Significant" = "#C5C3C6"
)

# Create volcano plot - X-axis: estimate (age effect), Y-axis: -log10(adjusted p-value)
volcano_plot <- ggplot(results_df, aes(x = estimate, y = -log10(adj_p_value))) +
  geom_point(aes(color = significance_category, size = abs(log2_fold_change)), alpha = 0.7) +
  scale_color_manual(values = color_palette) +
  scale_size_continuous(range = c(1, 4), name = "|Log2 Fold Change|") +
  
  # Add reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue", alpha = 0.7) +
  
  # Labels and theme
  labs(
    title = "Protein-Age Regression Volcano Plot",
    subtitle = "Age effect (regression coefficient) vs. Statistical significance",
    x = "Age Effect (Regression Coefficient)",
    y = "-Log10(Adjusted P-value)",
    color = "Significance Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save basic volcano plot
ggsave("protein_age_regression/volcano_plot_basic.pdf", volcano_plot, width = 12, height = 8, dpi = 300)
ggsave("protein_age_regression/volcano_plot_basic.png", volcano_plot, width = 12, height = 8, dpi = 300)

cat("Basic volcano plot saved to: protein_age_regression/volcano_plot_basic.pdf and .png\n")

# Create enhanced volcano plot with protein labels and GDF15 highlighting
cat("Creating enhanced volcano plot with protein labels...\n")

# Define criteria for labeling proteins
# Label proteins with high estimates (>10) and significant p-values (<0.05)
label_criteria <- results_df$estimate > 10 & results_df$adj_p_value < 0.05

# Create a subset for labeling
label_data <- results_df[label_criteria, ]

# Add GDF15 to label data if it exists and doesn't meet the criteria
if("GDF15" %in% results_df$protein) {
  gdf15_data <- results_df[results_df$protein == "GDF15", ]
  if(!(gdf15_data$estimate > 10 & gdf15_data$adj_p_value < 0.05)) {
    label_data <- rbind(label_data, gdf15_data)
  }
}

# Add GFAP to label data if it exists and doesn't meet the criteria
if("GFAP" %in% results_df$protein) {
  gfap_data <- results_df[results_df$protein == "GFAP", ]
  if(!(gfap_data$estimate > 10 & gfap_data$adj_p_value < 0.05)) {
    label_data <- rbind(label_data, gfap_data)
  }
}

# Add more proteins that aren't crowded - use a more comprehensive approach
# Include proteins with moderate estimates (>7) and high significance
moderate_estimate_proteins <- results_df %>%
  filter(estimate > 7, adj_p_value < 0.01) %>%
  filter(!protein %in% label_data$protein) %>%
  head(15)  # Limit to avoid overcrowding

# Add to label data
label_data <- rbind(label_data, moderate_estimate_proteins)

# Add some proteins with very high significance regardless of estimate
high_significance_proteins <- results_df %>%
  filter(adj_p_value < 1e-20) %>%
  filter(!protein %in% label_data$protein) %>%
  head(10)  # Limit to avoid overcrowding

# Add to label data
label_data <- rbind(label_data, high_significance_proteins)

# Remove duplicates
label_data <- label_data[!duplicated(label_data$protein), ]

# Create enhanced volcano plot
enhanced_volcano_plot <- ggplot(results_df, aes(x = estimate, y = -log10(adj_p_value))) +
  # Plot all points
  geom_point(aes(color = ifelse(protein == "GDF15", "GDF15", significance_category), 
                   size = abs(log2_fold_change)), alpha = 0.7) +
  
  # Custom color scale including GDF15
  scale_color_manual(values = c(
    "GDF15" = "#0066CC",  # Blue color for GDF15
    "Significant & High FC" = "#D7261E",
    "Significant" = "#F46036", 
    "High FC" = "#2E294E",
    "Not Significant" = "#C5C3C6"
  )) +
  
  scale_size_continuous(range = c(1, 4), name = "|Log2 Fold Change|") +
  
  # Add reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue", alpha = 0.7) +
  
  # Add protein labels
  geom_text_repel(
    data = label_data,
    aes(label = protein),
    size = 3,
    max.overlaps = 100,  # Increased to allow more labels
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.alpha = 0.6
  ) +
  
  # Labels and theme
  labs(
    title = "Protein-Age Regression Volcano Plot (Enhanced)",
    subtitle = "Age effect vs. Statistical significance with protein labels",
    x = "Age Effect (Regression Coefficient)",
    y = "-Log10(Adjusted P-value)",
    color = "Protein Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save enhanced volcano plot
ggsave("protein_age_regression/volcano_plot_enhanced.pdf", enhanced_volcano_plot, width = 14, height = 10, dpi = 300)
ggsave("protein_age_regression/volcano_plot_enhanced.png", enhanced_volcano_plot, width = 14, height = 8, dpi = 300)

cat("Enhanced volcano plot saved to: protein_age_regression/volcano_plot_enhanced.pdf and .png\n")

# Create a focused plot for high-estimate proteins
cat("Creating focused plot for high-estimate proteins...\n")

# Filter for proteins with high estimates (>5) and significant p-values (<0.05)
high_estimate_proteins <- results_df %>%
  filter(estimate > 5, adj_p_value < 0.05) %>%
  arrange(desc(estimate))

# Create focused volcano plot
focused_volcano_plot <- ggplot(high_estimate_proteins, aes(x = estimate, y = -log10(adj_p_value))) +
  geom_point(aes(color = ifelse(protein == "GDF15", "GDF15", "High Estimate"), 
                   size = abs(log2_fold_change)), alpha = 0.8) +
  
  # Custom color scale
  scale_color_manual(values = c(
    "GDF15" = "#0066CC",  # Blue color for GDF15
    "High Estimate" = "#2E86AB"
  )) +
  
  scale_size_continuous(range = c(2, 5), name = "|Log2 Fold Change|") +
  
  # Add protein labels
  geom_text_repel(
    aes(label = protein),
    size = 4,
    max.overlaps = 100,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "grey50",
    segment.alpha = 0.6
  ) +
  
  # Labels and theme
  labs(
    title = "High-Estimate Protein-Age Associations",
    subtitle = "Proteins with estimate > 5 and adj p < 0.05",
    x = "Age Effect (Regression Coefficient)",
    y = "-Log10(Adjusted P-value)",
    color = "Protein Type"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save focused volcano plot
ggsave("protein_age_regression/volcano_plot_focused.pdf", focused_volcano_plot, width = 16, height = 12, dpi = 300)
ggsave("protein_age_regression/volcano_plot_focused.png", focused_volcano_plot, width = 16, height = 12, dpi = 300)

cat("Focused volcano plot saved to: protein_age_regression/volcano_plot_focused.pdf and .png\n")

# Save the list of high-estimate proteins
write.csv(high_estimate_proteins, "protein_age_regression/high_estimate_proteins.csv", row.names = FALSE)
cat("High-estimate proteins list saved to: protein_age_regression/high_estimate_proteins.csv\n")

# Display summary of high-estimate proteins
cat("\n=== HIGH-ESTIMATE PROTEINS (estimate > 5, adj p < 0.05) ===\n")
print(high_estimate_proteins %>% select(protein, estimate, adj_p_value, log2_fold_change, r_squared))

# Create additional visualizations

# 1. P-value distribution
pvalue_plot <- ggplot(results_df, aes(x = p_value)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Distribution of P-values",
    x = "P-value",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

# 2. Adjusted P-value distribution
adj_pvalue_plot <- ggplot(results_df, aes(x = adj_p_value)) +
  geom_histogram(bins = 50, fill = "darkred", alpha = 0.7) +
  labs(
    title = "Distribution of FDR-Adjusted P-values",
    x = "Adjusted P-value",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

# 3. Fold change distribution
fc_plot <- ggplot(results_df, aes(x = log2_fold_change)) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  labs(
    title = "Distribution of Log2 Fold Changes",
    x = "Log2 Fold Change (Old/Young)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

# 4. Age effect vs Fold change correlation
correlation_plot <- ggplot(results_df, aes(x = estimate, y = log2_fold_change)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Age Effect vs. Fold Change Correlation",
    x = "Age Effect (Regression Coefficient)",
    y = "Log2 Fold Change (Old/Young)"
  ) +
  theme_minimal(base_size = 12)

# Combine plots
combined_plots <- (pvalue_plot + adj_pvalue_plot) / (fc_plot + correlation_plot)

# Save combined plots
ggsave("protein_age_regression/distribution_plots.pdf", combined_plots, width = 16, height = 12, dpi = 300)
ggsave("protein_age_regression/distribution_plots.png", combined_plots, width = 16, height = 12, dpi = 300)

cat("Distribution plots saved to: protein_age_regression/distribution_plots.pdf and .png\n")

# Create summary statistics table
summary_stats <- results_df %>%
  summarise(
    total_proteins = n(),
    significant_05 = sum(adj_p_value < 0.05),
    significant_01 = sum(adj_p_value < 0.01),
    significant_001 = sum(adj_p_value < 0.001),
    high_fc_05 = sum(abs(log2_fold_change) > 0.5),
    high_fc_1 = sum(abs(log2_fold_change) > 1),
    sig_and_high_fc = sum(adj_p_value < 0.05 & abs(log2_fold_change) > 0.5),
    mean_age_effect = mean(abs(estimate), na.rm = TRUE),
    median_age_effect = median(abs(estimate), na.rm = TRUE),
    mean_fold_change = mean(abs(log2_fold_change), na.rm = TRUE),
    median_fold_change = median(abs(log2_fold_change), na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Save summary statistics
write.csv(summary_stats, "protein_age_regression/summary_statistics.csv", row.names = FALSE)
cat("Summary statistics saved to: protein_age_regression/summary_statistics.csv\n")

# Display summary
cat("\n=== SUMMARY STATISTICS ===\n")
print(summary_stats)

# Top significant proteins
top_significant <- results_df %>%
  filter(adj_p_value < 0.05) %>%
  arrange(adj_p_value) %>%
  select(protein, estimate, p_value, adj_p_value, log2_fold_change, r_squared) %>%
  head(20)

cat("\n=== TOP 20 SIGNIFICANT PROTEINS ===\n")
print(top_significant)

# Save top significant proteins
write.csv(top_significant, "protein_age_regression/top_significant_proteins.csv", row.names = FALSE)
cat("Top significant proteins saved to: protein_age_regression/top_significant_proteins.csv\n")

cat("\nAnalysis completed successfully!\n")
cat("All results saved in: protein_age_regression/\n")
