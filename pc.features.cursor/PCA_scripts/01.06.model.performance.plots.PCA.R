#!/usr/bin/env Rscript

# Model Performance Plots for PCA Models
# This script generates performance plots similar to 01.04.model.performance.plots.R
# but adapted for PCA-based models

# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# Set working directory to project root
setwd("/Users/hhuang/Desktop/aging prediction/pc.features.cursor")

cat("Starting PCA Model Performance Plots Generation...\n")

# Function to generate model performance plots
generate_model_performance_plots <- function(model_type, test_fold_i = 1) {
  cat("Generating plots for", model_type, "model, test fold", test_fold_i, "\n")
  
  # Define the path to the model iteration results
  model_path <- paste0('PCA_results/model.iter.12/', model_type, '/both/both.model.', test_fold_i, '.csv')
  
  if (!file.exists(model_path)) {
    cat("Warning: Model file not found at", model_path, "\n")
    return(NULL)
  }
  
  # Read the model results
  iter_model <- read_csv(model_path)
  
  # Clean column names and handle data types
  colnames(iter_model)[1] <- "iteration"
  
  # Convert AUC to numeric if it's logical
  if (is.logical(iter_model$auc)) {
    iter_model$auc <- as.numeric(iter_model$auc)
  }
  
  # Remove rows with missing values
  iter_model <- iter_model %>% filter(!is.na(rmse) & !is.na(auc))
  
  if (nrow(iter_model) == 0) {
    cat("Warning: No valid data found for", model_type, "\n")
    return(NULL)
  }
  
  # Find the best RMSE
  best_rmse_row <- iter_model %>% filter(rmse == min(rmse, na.rm = TRUE))
  cat("Best RMSE:", best_rmse_row$rmse, "\n")
  
  # Create scaled metrics for comparison
  iter_model_scale <- cbind(iter_model, iter_model %>% select(rmse, auc, error) %>% scale()) %>% 
    data.frame() %>% 
    as_tibble()
  
  # Rename scaled columns
  colnames(iter_model_scale)[(ncol(iter_model_scale)-2):ncol(iter_model_scale)] <- c("rmse_scaled", "auc_scaled", "error_scaled")
  
  # Create a combined metric (lower RMSE and higher AUC is better)
  iter_model_scale <- iter_model_scale %>% 
    mutate(sum_metric = -rmse_scaled + auc_scaled)
  
  # Find the best combined metric
  best_combined_row <- iter_model_scale %>% filter(sum_metric == max(sum_metric, na.rm = TRUE))
  cat("Best combined metric:", best_combined_row$sum_metric, "\n")
  
  # Convert to long format for plotting
  iter_model_long <- iter_model %>%
    pivot_longer(cols = c(rmse, auc, error), names_to = 'metric')
  
  # Create performance plots
  plots <- list()
  
  # Plot 1: Performance metrics over iterations
  p1 <- iter_model_long %>%
    ggplot(aes(x = 1:nrow(.), y = value, color = metric)) +
    geom_point(alpha = 0.7) +
    geom_line(alpha = 0.5) +
    labs(title = paste(model_type, "- Performance Metrics Over Iterations"),
         subtitle = paste("Test Fold", test_fold_i),
         x = "Iteration", y = "Value", color = "Metric") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plots[[1]] <- p1
  
  # Plot 2: Scaled metrics comparison
  p2 <- iter_model_scale %>%
    arrange(desc(sum_metric)) %>%
    mutate(rank = 1:n()) %>%
    mutate(auc = auc * 100) %>%
    pivot_longer(cols = c(rmse, auc, error), names_to = 'metric') %>%
    ggplot() +
    geom_point(aes(rank, value, color = metric), alpha = 0.7) +
    labs(title = paste(model_type, "- Scaled Metrics by Rank"),
         subtitle = paste("Test Fold", test_fold_i, "- Ranked by Combined Metric"),
         x = "Rank", y = "Scaled Value", color = "Metric") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plots[[2]] <- p2
  
  # Plot 3: Hyperparameter analysis (if available)
  if ("alpha" %in% colnames(iter_model) && "lambda" %in% colnames(iter_model)) {
    p3 <- iter_model %>%
      ggplot(aes(x = alpha, y = lambda, color = rmse, size = auc)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "red", high = "blue") +
      scale_size_continuous(range = c(1, 5)) +
      labs(title = paste(model_type, "- Hyperparameter Space"),
           subtitle = paste("Test Fold", test_fold_i, "- Color: RMSE, Size: AUC"),
           x = "Alpha", y = "Lambda", color = "RMSE", size = "AUC") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plots[[3]] <- p3
  } else if ("gamma" %in% colnames(iter_model) && "cost" %in% colnames(iter_model)) {
    # For SVM models
    p3 <- iter_model %>%
      ggplot(aes(x = gamma, y = cost, color = rmse, size = auc)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "red", high = "blue") +
      scale_size_continuous(range = c(1, 5)) +
      labs(title = paste(model_type, "- Hyperparameter Space"),
           subtitle = paste("Test Fold", test_fold_i, "- Color: RMSE, Size: AUC"),
           x = "Gamma", y = "Cost", color = "RMSE", size = "AUC") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plots[[3]] <- p3
  } else if ("learning_rate" %in% colnames(iter_model) && "num_leaves" %in% colnames(iter_model)) {
    # For LightGBM models
    p3 <- iter_model %>%
      ggplot(aes(x = learning_rate, y = num_leaves, color = rmse, size = auc)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "red", high = "blue") +
      scale_size_continuous(range = c(1, 5)) +
      labs(title = paste(model_type, "- Hyperparameter Space"),
           subtitle = paste("Test Fold", test_fold_i, "- Color: RMSE, Size: AUC"),
           x = "Learning Rate", y = "Num Leaves", color = "RMSE", size = "AUC") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plots[[3]] <- p3
  }
  
  # Plot 4: RMSE distribution
  p4 <- iter_model %>%
    ggplot(aes(x = rmse)) +
    geom_histogram(bins = 30, fill = "#D7261E", alpha = 0.7) +
    geom_vline(xintercept = best_rmse_row$rmse, color = "blue", linetype = "dashed", size = 1) +
    labs(title = paste(model_type, "- RMSE Distribution"),
         subtitle = paste("Test Fold", test_fold_i, "- Blue line: Best RMSE"),
         x = "RMSE", y = "Count") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plots[[4]] <- p4
  
  # Combine plots
  if (length(plots) == 4) {
    combined_plot <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])
  } else {
    combined_plot <- (plots[[1]] + plots[[2]]) / plots[[4]]
  }
  
  # Save the combined plot
  output_path <- paste0('PCA_results/model.iter.12/', model_type, '/both/performance_plots_', test_fold_i, '.pdf')
  ggsave(output_path, combined_plot, width = 16, height = 12, dpi = 300)
  cat("Performance plots saved to:", output_path, "\n")
  
  return(combined_plot)
}

# Generate plots for each model type
cat("\n=== Generating Performance Plots for All Models ===\n")

# Elastic Net
cat("\n--- Elastic Net Model ---\n")
elastic_net_plots <- generate_model_performance_plots("elastic.net.PCA", test_fold_i = 1)

# SVM
cat("\n--- SVM Model ---\n")
svm_plots <- generate_model_performance_plots("SVM.PCA", test_fold_i = 1)

# LightGBM
cat("\n--- LightGBM Model ---\n")
lightgbm_plots <- generate_model_performance_plots("lightGBM.PCA", test_fold_i = 1)

# Generate summary statistics
cat("\n=== Summary Statistics ===\n")

# Function to get summary stats for each model
get_model_summary <- function(model_type, test_fold_i = 1) {
  model_path <- paste0('PCA_results/model.iter.12/', model_type, '/both/both.model.', test_fold_i, '.csv')
  
  if (!file.exists(model_path)) {
    return(NULL)
  }
  
  iter_model <- read_csv(model_path)
  
  # Clean column names
  colnames(iter_model)[1] <- "iteration"
  
  # Convert AUC to numeric if it's logical
  if (is.logical(iter_model$auc)) {
    iter_model$auc <- as.numeric(iter_model$auc)
  }
  
  # Remove rows with missing values
  iter_model <- iter_model %>% filter(!is.na(rmse) & !is.na(auc))
  
  if (nrow(iter_model) == 0) {
    return(NULL)
  }
  
  summary_stats <- iter_model %>%
    summarise(
      model = model_type,
      min_rmse = min(rmse, na.rm = TRUE),
      max_rmse = max(rmse, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      median_rmse = median(rmse, na.rm = TRUE),
      min_auc = min(auc, na.rm = TRUE),
      max_auc = max(auc, na.rm = TRUE),
      mean_auc = mean(auc, na.rm = TRUE),
      median_auc = median(auc, na.rm = TRUE),
      n_iterations = n()
    ) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  
  return(summary_stats)
}

# Get summaries for all models
model_summaries <- list(
  get_model_summary("elastic.net.PCA"),
  get_model_summary("SVM.PCA"),
  get_model_summary("lightGBM.PCA")
)

# Combine summaries
all_summaries <- bind_rows(model_summaries[!sapply(model_summaries, is.null)])

if (nrow(all_summaries) > 0) {
  cat("\nModel Performance Summary:\n")
  print(all_summaries)
  
  # Save summary to CSV
  summary_path <- 'PCA_results/model.iter.12/model_performance_summary.csv'
  write.csv(all_summaries, summary_path, row.names = FALSE)
  cat("\nSummary saved to:", summary_path, "\n")
}

cat("\n=== Performance Plots Generation Complete ===\n")
cat("All plots have been saved to their respective model folders.\n")
