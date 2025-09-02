rm(list = ls())

library(tidyverse)
library(ggplot2)
library(gridExtra)

# This script compares the performance of PCA-based models vs. original models
# Run this after executing both the original and PCA-based scripts

cat("PCA vs. Original Model Performance Comparison\n")
cat("============================================\n\n")

# Function to extract results from model files
extract_model_results <- function(model_type, is_pca = FALSE) {
  suffix <- ifelse(is_pca, ".PCA", "")
  if(is_pca) {
    pattern <- paste0("PCA_results/model.iter.12/", model_type, suffix, "/both/")
  } else {
    pattern <- paste0("figures/model.iter.12/", model_type, suffix, "/both/")
  }

  # Look for model result files
  result_files <- list.files(pattern = "*.csv", path = pattern, full.names = TRUE, recursive = TRUE)

  if(length(result_files) == 0) {
    cat("No result files found for", model_type, ifelse(is_pca, "(PCA)", "(Original)"), "\n")
    return(NULL)
  }

  # Read and combine results
  all_results <- list()
  for(file in result_files) {
    results <- read.csv(file)
    fold_num <- str_extract(file, "model\\.(\\d+)\\.csv") %>% str_extract("\\d+")
    results$fold <- fold_num
    results$file <- file
    all_results[[length(all_results) + 1]] <- results
  }

  if(length(all_results) == 0) return(NULL)

  combined_results <- bind_rows(all_results)
  combined_results$model_type <- model_type
  combined_results$is_pca <- is_pca

  return(combined_results)
}

# Extract results for each model type
cat("Extracting results...\n")

# Elastic Net
elastic_net_original <- extract_model_results("elastic.net", FALSE)
elastic_net_pca <- extract_model_results("elastic.net.PCA", TRUE)
elastic_net_pca

# SVM
svm_original <- extract_model_results("svm", FALSE)
svm_pca <- extract_model_results("svm.PCA", TRUE)
svm_pca

# LightGBM
lightgbm_original <- extract_model_results("lightgbm", FALSE)
lightgbm_pca <- extract_model_results("lightgbm.PCA", TRUE)

# Combine all results
all_results <- bind_rows(
  elastic_net_original, elastic_net_pca,
  svm_original, svm_pca,
  lightgbm_original, lightgbm_pca
)

if(!is.null(all_results) && nrow(all_results) > 0) {

  # Create summary statistics
  summary_stats <- all_results %>%
    group_by(model_type, is_pca) %>%
    summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      mean_cor = mean(cor, na.rm = TRUE),
      sd_cor = sd(cor, na.rm = TRUE),
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      n_models = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      model_label = paste0(model_type, ifelse(is_pca, " (PCA)", " (Original)")),
      pca_status = ifelse(is_pca, "PCA", "Original")
    )

  # Print summary
  cat("\nModel Performance Summary:\n")
  cat("==========================\n")
  print(summary_stats)

  # Create comparison plots
  cat("\nCreating comparison plots...\n")

  # RMSE comparison
  p1 <- ggplot(summary_stats, aes(x = model_type, y = mean_rmse, fill = pca_status)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_rmse - sd_rmse, ymax = mean_rmse + sd_rmse),
                  position = position_dodge(0.9), width = 0.25) +
    labs(title = "Model Performance Comparison: RMSE",
         x = "Model Type", y = "Mean RMSE", fill = "Data Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Correlation comparison
  p2 <- ggplot(summary_stats, aes(x = model_type, y = mean_cor, fill = pca_status)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_cor - sd_cor, ymax = mean_cor + sd_cor),
                  position = position_dodge(0.9), width = 0.25) +
    labs(title = "Model Performance Comparison: Correlation",
         x = "Model Type", y = "Mean Correlation", fill = "Data Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # AUC comparison (if available)
  if(all(!is.na(summary_stats$mean_auc))) {
    p3 <- ggplot(summary_stats, aes(x = model_type, y = mean_auc, fill = pca_status)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                    position = position_dodge(0.9), width = 0.25) +
      labs(title = "Model Performance Comparison: AUC",
           x = "Model Type", y = "Mean AUC", fill = "Data Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Combine all plots
    combined_plot <- grid.arrange(p1, p2, p3, ncol = 2)
  } else {
    # Combine RMSE and correlation plots only
    combined_plot <- grid.arrange(p1, p2, ncol = 2)
  }

  # Save plots
  dir.create('PCA_results', recursive = TRUE, showWarnings = FALSE)
  ggsave("PCA_results/PCA_vs_Original_Comparison.png", combined_plot,
         width = 12, height = 8, dpi = 300)

  # Statistical comparison
  cat("\nStatistical Comparison:\n")
  cat("======================\n")

  # Compare PCA vs Original for each model type
  for(model in unique(summary_stats$model_type)) {
    original_data <- summary_stats %>% filter(model_type == model & !is_pca)
    pca_data <- summary_stats %>% filter(model_type == model & is_pca)

    if(nrow(original_data) > 0 && nrow(pca_data) > 0) {
      cat("\n", model, ":\n")
      cat("  RMSE: Original =", round(original_data$mean_rmse, 4),
          "vs PCA =", round(pca_data$mean_rmse, 4))
      if(original_data$mean_rmse > pca_data$mean_rmse) {
        cat(" (PCA better by", round((original_data$mean_rmse - pca_data$mean_rmse) / original_data$mean_rmse * 100, 2), "%)\n")
      } else {
        cat(" (Original better by", round((pca_data$mean_rmse - original_data$mean_rmse) / pca_data$mean_rmse * 100, 2), "%)\n")
      }

      cat("  Correlation: Original =", round(original_data$mean_cor, 4),
          "vs PCA =", round(pca_data$mean_cor, 4))
      if(pca_data$mean_cor > original_data$mean_cor) {
        cat(" (PCA better by", round((pca_data$mean_cor - original_data$mean_cor) / abs(original_data$mean_cor) * 100, 2), "%)\n")
      } else {
        cat(" (Original better by", round((original_data$mean_cor - pca_data$mean_cor) / abs(pca_data$mean_cor) * 100, 2), "%)\n")
      }
    }
  }

  # Save summary results
  write.csv(summary_stats, "PCA_results/PCA_vs_Original_Summary.csv", row.names = FALSE)
  cat("\nResults saved to 'PCA_results/PCA_vs_Original_Summary.csv'\n")

} else {
  cat("No results found. Please run the original and PCA-based model scripts first.\n")
}

cat("\nComparison completed!\n")
