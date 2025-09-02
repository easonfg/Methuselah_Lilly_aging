#####################################################
## Update 5/21 to include ratios
## pSNCA-129/SNCA, Oligo-SNCA/SNCA, pTau-181/MAPT, 
## pTau-271/MAPT, pTau-231/MAPT, Aβ42/Aβ40 ratios
####################################################
library(NULISAseqR)
library(xlsx)
library(tidyverse)
library(TAPGeneralization)
library(ggpubr)
library(gridExtra)
library(purrr)
library(readr)
library(ggVennDiagram)
library(corrplot)

# setting ggplot theme
# custom ggplot theme

custom_theme <- theme_bw() + theme(
  panel.background = element_rect(fill='white'),
  plot.background = element_rect(fill='transparent', color = NA),
  legend.background = element_rect(fill='transparent'),
  legend.key = element_rect(fill = "transparent", color = NA)
)

theme_set(custom_theme)

############################
## Read Data
############################
## set workingDir
setwd("~/TAP-P000676-UCL-CNS")

source("scripts/heatmap_pca_functions.R")

## process metadata and save a file
metadata_676 <- readxl::read_xlsx("data/TAP_Sample_Info_UCL_CNS_P1.xlsx") %>%
  rename(SampleName = sampleName,
         PlateID = plateID) %>% 
  select(-`Replicate #`, -SampleID) %>%
  filter(!grepl("NC|IPC", SampleName)) %>% 
  mutate(Bridged_sample = ifelse(grepl("Visit", SampleName), 1, 0),
         SampleMatrix = ifelse(grepl("P|Plasma|SC", SampleName), "Plasma", "Serum"),
         SampleID = str_extract(SampleName, "(?<=__).+?(?=(_Serum|_Plasma| P| S|$))"),
         Sex = ifelse(Sex == "F", "Female", "Male"),
         Group = case_when(
           Bridged_sample == 1 ~ NA,
           grepl("SC", SampleName) ~ "SC",
           TRUE ~"Healthy control")) 

bridge_samples <- metadata_676 %>% filter(Bridged_sample == 1) %>% pull(SampleID)

metadata_535 <- readxl::read_xlsx("data/TAP_Sample_Info_P535_UCL-Crick_CNS.xlsx") %>%
  rename(SampleName = sampleName,
         PlateID = plateID) %>% 
  select(-`Replicate #`, -SampleID) %>%
  filter(!grepl("NC|IPC", SampleName)) %>% 
  mutate(string_split = SampleName) %>%
  separate(string_split, into = c("Prefix", "WellNum", "SampleID", "Timepoint", "SampleMatrix"), sep = "_") %>%
  select(-Prefix, -WellNum) %>% 
  # filter(SampleMatrix != "CSF") %>% 
  mutate(SampleID = paste0(SampleID, "_", Timepoint),
         Bridged_sample = ifelse(
           SampleMatrix %in% c("Plasma", "Serum") &
             SampleID %in% bridge_samples, 1, 0))

metadata_535_1 <- read_csv("data/UCL.Aliquot.Information.and.metadata.csv") %>% 
  rename(SampleMatrix = `Sample Type`,
         Age = `Age at Baseline`,
         Sex = Gender,
         Group = Disease) %>% 
  select(SampleMatrix,`Sample ID`, `Visit no.`, Group, Randomisation, Sex, Age,
         `Diagnosis Duration`, Visit, Score) %>% 
  mutate(SampleID = paste0(`Sample ID`, "_", `Visit no.`))

metadata_535 <- metadata_535 %>% 
  left_join(metadata_535_1) %>% 
  mutate(Group = ifelse(grepl("SC", SampleName), "SC", Group),
         SampleMatrix = ifelse(Group == "SC", "Plasma", SampleMatrix),
         Timepoint = ifelse(Group == "SC", NA, Timepoint))

metadata_535 %>% write_csv("data/metadata_535.csv")

# Perform the join to add missing info to metadata_676 from metadata_535
metadata_676 <- metadata_676 %>%
  left_join(
    metadata_535 %>% 
      filter(Bridged_sample == 1) %>% 
      select(SampleID, Group, Sex, Age), 
    by = "SampleID") %>%
  # Combine the information (keep original if no match found)
  mutate(
    Group = coalesce(Group.x, Group.y),
    Sex = coalesce(Sex.x, Sex.y),
    Age = coalesce(Age.x, Age.y)
  ) %>%
  # Remove the temporary columns
  select(-ends_with(".x"), -ends_with(".y")) %>% 
  distinct()

metadata_676 %>% write_csv("data/metadata_676.csv")

## read data
## include SC too
process_data <- function(data_dir, file_name, metadata_file_name){
  
  data_long <- read_csv(file.path(data_dir, file_name)) %>% 
    filter(!grepl("NC|IPC", SampleName)) %>% 
    mutate(SampleName_long = paste0(PlateID, "_", SampleName))
  
  data_wide <- data_long %>% 
    select(Target, SampleName_long, NPQ) %>% 
    distinct(.keep_all = T) %>% 
    pivot_wider(names_from = SampleName_long, values_from = NPQ) %>%  
    column_to_rownames("Target")
  
  data_wide_lod <- data_long %>% 
    select(Target, SampleName_long, LOD) %>% 
    distinct(.keep_all = T) %>% 
    pivot_wider(names_from = SampleName_long, values_from = LOD) %>%  
    column_to_rownames("Target")
  
  metadata <- read_csv(file.path(data_dir, metadata_file_name)) %>% 
    mutate(SampleName_long = paste0(PlateID, "_",SampleName))
  
  bridge_metadata <- metadata %>% filter(Bridged_sample == 1)
  
  data_long <- data_long %>% left_join(metadata)
  
  return(list(data_long = data_long,
              data_wide = data_wide,
              data_wide_lod = data_wide_lod,
              metadata = metadata,
              bridge_metadata = bridge_metadata))
}

project_535 <- process_data(data_dir = "data", 
                            file_name = "20240816_P-000535_UCL-Crick_NULISAseq_TAP_Counts_Report_FULL.csv", 
                            metadata_file_name = "metadata_535.csv")
project_676 <- process_data(data_dir = "data", 
                            file_name = "20241120_P-000676_UCL_NULISAseq_TAP_Counts_Report_FULL.csv", 
                            metadata_file_name = "metadata_676.csv")

############################
## Bridge Normalization
############################
# function to perform bridge normalization from run2 to run1
perform_bridge_norm <- function(run1_data, run2_data, bridge_ids) {
  # Calculate bridge terms
  bridge_terms <- apply(run1_data[, bridge_ids$run1] - run2_data[, bridge_ids$run2], 1, median)
  
  # Apply normalization
  run2_norm <- run2_data + bridge_terms
  
  # Combine results
  list(
    run2_data_wide_bridge = run2_norm,
    data_wide_bridge = bind_cols(run1_data, run2_norm)
  )
}

## make a bridge sample only metadata from both projects
bridge_metadata <- bind_rows(
  project_535$metadata %>% filter(Bridged_sample == 1) %>% mutate(Project = "P-000535"),
  project_676$metadata %>% filter(Bridged_sample == 1) %>% mutate(Project = "P-000676") 
)

# function to process each matrix type
process_matrix_type <- function(matrix_type, bridge_sample_ids) {
  # Get samples for each run
  run1_samples <- project_535$metadata %>%
    filter(SampleMatrix == matrix_type) %>%
    pull(SampleName_long)
  
  run2_samples <- project_676$metadata %>%
    filter(SampleMatrix == matrix_type) %>%
    pull(SampleName_long)
  
  # Get bridge samples
  bridge_pairs <- bridge_metadata %>%
    filter(SampleMatrix == matrix_type,
           SampleID %in% bridge_sample_ids) %>%
    select(SampleID, Project, SampleName_long) %>%
    pivot_wider(names_from = "Project", values_from = "SampleName_long")
  
  # Perform normalization
  npq <- perform_bridge_norm(
    run1_data = project_535$data_wide[, run1_samples],
    run2_data = project_676$data_wide[, run2_samples],
    bridge_ids = list(run1 = bridge_pairs$`P-000535`, 
                      run2 = bridge_pairs$`P-000676`)
  )
  
  lod <- perform_bridge_norm(
    run1_data = project_535$data_wide_lod[, run1_samples],
    run2_data = project_676$data_wide_lod[, run2_samples],
    bridge_ids = list(run1 = bridge_pairs$`P-000535`, 
                      run2 = bridge_pairs$`P-000676`)
  )
  
  return(list(NPQ_bridge = npq,
              LOD_bridge = lod))
}

# Define bridge samples for each matrix type
bridge_samples <- list(
  Plasma = c("EX-048_Visit 6", "EX-067_Visit 2", "EX-116_Visit 10", "EX-190_Visit 2"),
  Serum = c("01-39_Visit 2", "01-20_Visit 2", "01-27_Visit 4")
)

# Process all matrix types at once
bridged_results <- lapply(names(bridge_samples), function(matrix_type) {
  process_matrix_type(matrix_type, bridge_samples[[matrix_type]])
}) %>%
  setNames(names(bridge_samples))

data_wide_bridged <- do.call(cbind, list(bridged_results$Plasma$NPQ_bridge$data_wide_bridge, 
                                         bridged_results$Serum$NPQ_bridge$data_wide_bridge))
data_wide_lod_bridged <- do.call(cbind, list(bridged_results$Plasma$LOD_bridge$data_wide_bridge, 
                                             bridged_results$Serum$LOD_bridge$data_wide_bridge))
############################
## Combine data from two projects
############################
## combine metadata
metadata <- bind_rows(
  project_535$metadata %>% mutate(Project = "P-000535"),
  project_676$metadata %>% mutate(Project = "P-000676") 
) %>% 
  mutate(Project_PlateID = paste(Project, PlateID))

metadata_1setb <- bind_rows(
  project_535$metadata %>% mutate(Project = "P-000535"),
  project_676$metadata %>% filter(Bridged_sample == 0) %>% mutate(Project = "P-000676") 
) %>% 
  mutate(Project_PlateID = paste(Project, PlateID))

# bridged samples from project 676 not in metadata

## combine data wide
data_wide <- bind_cols(project_535$data_wide, project_676$data_wide)

## combine data long
data_long <- bind_rows(
  project_535$data_long %>% mutate(Project = "P-000535"),
  project_676$data_long %>% mutate(Project = "P-000676")
)

project_all <- list(
  data_long = data_long,
  data_wide = data_wide, 
  metadata = metadata,
  metadata_535_bridge_only = metadata_1setb,
  data_wide_bridged = data_wide_bridged
) 

## convert bridged data wide to long 
data_long_bridged <- project_all$data_wide_bridged %>% 
  rownames_to_column("Target") %>% 
  pivot_longer(
    cols = -Target, 
    names_to = "SampleName_long", 
    values_to = "NPQ_bridge_norm"
  ) %>% 
  mutate(SampleType = ifelse(grepl("SC", SampleName_long), "SC", "Sample"))

## convert bridged data wide load to long 
t <- data_wide_lod_bridged %>% 
  rownames_to_column("Target") %>% 
  pivot_longer(
    cols = -Target, 
    names_to = "SampleName_long", 
    values_to = "LOD_bridge_norm"
  ) %>% 
  mutate(SampleType = ifelse(grepl("SC", SampleName_long), "SC", "Sample"))

project_all$data_long_bridged <- data_long_bridged %>% 
  left_join(t) %>% 
  left_join(project_all$data_long)
# select(all_of(names(project_535$data_long))) 

# bridged data does not include CSF
# note: outlier from QC, 01-24_Visit 6_Serum not good in both study

#######################
## Overall Detectability 
#######################
# Calculate overall detectability 
detectability <- project_all$data_long_bridged %>%
  filter(SampleType == "Sample") %>% 
  mutate(above_LOD = NPQ_bridge_norm > LOD_bridge_norm) %>%
  group_by(SampleMatrix) %>%
  mutate(matrix_sample_count = n_distinct(SampleName_long)) %>%
  ungroup() %>%
  group_by(Target, SampleMatrix, matrix_sample_count) %>%
  summarize(
    detectability = round(sum(above_LOD, na.rm = TRUE) / n() * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(
    detectability = ifelse(Target %in% c("APOE", "CRP"), 100.0, detectability), # reverse curve targets
    SampleMatrix = paste0("Overall: ", SampleMatrix, " (n = ", matrix_sample_count, ")")
  ) %>%
  # Convert to wide format
  select(Target, SampleMatrix, detectability) %>%
  pivot_wider(
    names_from = SampleMatrix,
    values_from = detectability,
    values_fill = NA  # Use NA for missing combinations (or 0 if preferred)
  ) 

# project 535
detect_535 <- read_csv("figures/target_detectability_table_535.csv") %>% 
  select(1:4)
colnames(detect_535) <- c("Target", "Project 535: Plasma (n = 483)", "Project 535: CSF (n = 148)", "Project 535: Serum (n = 113)")

# project 676
detect_676 <- read_csv("figures/target_detectability_table.csv")
colnames(detect_676) <- c("Target", "Project 676: Serum (n = 43)", "Project 676: Plasma (n = 43)")

## Combine detectability tables
project_all$detectability <- detectability %>% left_join(detect_535) %>% left_join(detect_676)

#######################
## Detectability 
## cut off >= 50%
#######################
project_all$detectability %>% 
  filter(`Overall: Plasma (n = 526)` >= 50) %>% 
  pull(Target) -> target_passed_detct_qc_p
# 118 targets (124 targets)
# Exclude: APOE4, GDI1, PTN, SNCB, UCHL1, YWHAZ

project_all$detectability %>% 
  filter(`Overall: Serum (n = 156)` >= 50) %>% 
  pull(Target) -> target_passed_detct_qc_s
# 116 targets (124 targets)
# Exclude: APOE4, GDI1, IL1B, PTN, SNCB, UCHL1, YWHAZ, pTDP43-409

project_all$detectability %>% 
  filter(`Project 535: CSF (n = 148)` >= 50, Target != "APOE4") %>% 
  pull(Target) -> target_passed_detct_qc_c
# 102 targets (124 targets)
# Exclude: APOE4, ARSA, BASP1, BDNF, GDNF, HBA1, IFNG, IL17A, IL1B, IL2, IL33, IL4, 
# MME, Oligo-SNCA, PARK7, PDLIM5, S100A12, SAA1, VEGFD, VGF, YWHAZ, pTDP43-409


############################
## Add ratio as artificial targets
############################
ratio_pairs <- list(
  c("pSNCA-129", "SNCA"),
  c("Oligo-SNCA", "SNCA"),
  c("pTau-181", "MAPT"),
  c("pTau-217", "MAPT"),
  c("pTau-231", "MAPT"),
  c("Aβ42", "Aβ40")
)

## wide data format
project_all$data_wide_bridged_ratio <- project_all$data_wide_bridged
# Iterate over ratio_pairs to create new rows
for (pair in ratio_pairs) {
  # Define the new target row name based on the pair
  new_target <- paste(pair[1], "/", pair[2], sep = "")
  
  # Perform the subtraction for each sample
  project_all$data_wide_bridged_ratio[new_target, ] <- project_all$data_wide_bridged[pair[1], ] - project_all$data_wide_bridged[pair[2], ]
}

## long data format
project_all$data_long_bridged_ratio <- project_all$data_long_bridged

# data_long_bridged_ratio <- project_all$data_wide_bridged_ratio %>% 
#   rownames_to_column("Target") %>% 
#   pivot_longer(
#     cols = -Target, 
#     names_to = "SampleName_long", 
#     values_to = "NPQ_bridge_norm"
#   ) %>% 
#   mutate(SampleType = ifelse(grepl("SC", SampleName_long), "SC", "Sample"))

new_rows <- list()
# Iterate over the ratio_pairs
for (pair in ratio_pairs) {
  # Filter for the two targets in the pair
  data_subset <- project_all$data_long_bridged %>%
    filter(Target %in% pair)
  
  # Create the new Target name (e.g., pSNCA-129/SNCA)
  new_target <- paste(pair[1], "/", pair[2], sep = "")
  
  # Perform the subtraction of NPQ_bridge_norm for the pair
  data_new <- data_subset %>%
    group_by(SampleName_long) %>%
    summarise(NPQ_bridge_norm = NPQ_bridge_norm[Target == pair[1]] - NPQ_bridge_norm[Target == pair[2]]) %>%
    mutate(Target = new_target) 
  
  # Append the new rows to the list
  new_rows[[new_target]] <- data_new
}
new_rows_combined <- bind_rows(new_rows)
# Add metadata to the subset data frame
new_rows_combined <- new_rows_combined %>% left_join(project_all$metadata_535_bridge_only)
# Combine the original data with the new rows
project_all$data_long_bridged_ratio <- bind_rows(project_all$data_long_bridged, new_rows_combined)




############################
## Differential Expression
############################
## update from 3/26/25 did not account for different sample matrix in healthy control and updated the analysis 
## update to add ratio as artificial targets
## exclude outlier B_11_01-24_Visit 2_Serum_UCL_P8 serum
project_all$metadata_ex <- project_all$metadata %>% filter(SampleName != "B_11_01-24_Visit 2_Serum_UCL_P8")
project_all$metadata_535_bridge_only_ex <- project_all$metadata_535_bridge_only %>% 
  filter(SampleName != "B_11_01-24_Visit 2_Serum_UCL_P8")

target_ratios <- c("pSNCA-129/SNCA", "Oligo-SNCA/SNCA", "pTau-181/MAPT", "pTau-217/MAPT", 
                   "pTau-231/MAPT", "Aβ42/Aβ40")
## DE
metadata_analysis <- project_all$metadata_535_bridge_only_ex %>% 
  filter(`Visit no.` %in% c(NA, "Visit 2"),
         Group != "SC", 
         SampleMatrix != "CSF")

control <- "Healthy control"
comparison <- list(
  case = c("PD", "MSA-C", "MSA-P", "MSA"),
  sample_matrix = c("Plasma", rep("Serum", 3))
)

for (i in seq_along(comparison$case)) {
  case <- comparison$case[i]
  sample_matrix <- comparison$sample_matrix[i]
  
  if(sample_matrix == "Plasma") {
    targets <- c(target_passed_detct_qc_p, target_ratios)
  } else {
    targets <- c(target_passed_detct_qc_s, target_ratios)
  }
  if(case == "MSA"){
    metadat <- metadata_analysis %>%
      mutate(Group = ifelse(grepl("MSA", Group), "MSA", Group)) %>%
      filter(Group %in% c(control, case),
             SampleMatrix == sample_matrix) %>% 
      mutate(Group = factor(Group, levels = c(control, case)))
    
    dat <- project_all$data_long_bridged_ratio %>%
      mutate(Group = ifelse(grepl("MSA", Group), "MSA", Group)) %>%
      filter(Target %in% targets,
             SampleName_long %in% metadat$SampleName_long) %>% 
      mutate(Group = factor(Group, levels = c(control, case))) %>% 
      select(-NPQ) %>% 
      rename(NPQ = NPQ_bridge_norm) # use bridged norm NPQ
    
  } else {
    metadat <- metadata_analysis %>%
      filter(Group %in% c(control, case),
             SampleMatrix == sample_matrix) %>% 
      mutate(Group = factor(Group, levels = c(control, case)))
    
    dat <- project_all$data_long_bridged_ratio %>%
      filter(Target %in% targets,
             SampleName_long %in% metadat$SampleName_long) %>% 
      mutate(Group = factor(Group, levels = c(control, case))) %>% 
      select(-NPQ) %>% 
      rename(NPQ = NPQ_bridge_norm) # use bridged norm NPQ
  }
  
  metadat %>% filter(Group == case) %>% nrow() -> ncase
  metadat %>% filter(Group == control) %>% nrow() -> ncontrol
  
  cat(control, "=", ncontrol, ",", case, "=", ncase, "\n")
  
  TAPGeneralization::run.models.w.cov(
    full.data = dat,
    full.metadata = metadat,
    sample_col = "SampleName_long",
    interaction.term = FALSE,
    model.type = "linear",
    col.interest = "Group",
    cov.col = "Sex",
    other.cov = "Age",
    subject.id = NULL,
    output.folder = paste0("figures/DE/Updated_Exclude_outlier_with_ratio/", case, "_vs_", control, "_")
  )
}

# Remake volcano plots
indi_volcano_plot <- function(model_type = "linear.with.cov", case, control, suffix = "FDR"){
  file_prefix <- paste0("figures/DE/Updated_Exclude_outlier_with_ratio/", case, "_vs_", control, "_", model_type)
  text <- paste0("Group", gsub("-", "_", case))
  
  output_prefix <- paste0(file_prefix, "/Group Effect/")
  
  de <- read_csv(paste0(output_prefix, "csv/model.res.csv"))
  
  coefs <- de[[paste0(text, "_coef")]]
  p_vals <- de[[paste0(text, "_pval_", suffix)]]
  
  target_labels <- de$target
  sig_label1 <- paste0(suffix, " = 5%")
  if(suffix == "FDR"){ylabel1 = expression("-log"[10] * "(FDR-adjusted p-value)")}
  else{ylabel1 = expression("-log"[10] * "(Unadjusted p-value)")}
  title <- paste0(case, " vs ", control)
  
  pdf(paste0(output_prefix, "volcano.plots/", case, "_pval_", suffix, ".volcano.plot.pdf"), width = 6, height = 4.5)
  print(volcanoPlot(coefs, p_vals, target_labels, title, sig_label = sig_label1, ylabel = ylabel1, plot_title_font_size = 11))
  dev.off()
  
  return(de)
}
de_results <- list()
# Loop through each row of the data frame 
for(case in comparison$case){
  # Call your indi_volcano_plot function
  de <- indi_volcano_plot(model_type = "linear.with.cov", case, control)
  
  # Create a name by concatenating model_type, case, and control
  result_name <- paste(case, control, sep = "_")
  
  # Store the result and assign the concatenated name to de_results
  de_results[[result_name]] <- de
}

############################
## Heatmap & PCA
############################
metadata_analysis_no_ex <- project_all$metadata_535_bridge_only %>% 
  filter(`Visit no.` %in% c(NA, "Visit 2"),
         Group != "SC", 
         SampleMatrix != "CSF") %>% 
  mutate(
           Age_Category = cut(
             Age,
             breaks = seq(33, 90, length.out = 9),
             include.lowest = TRUE,
             labels = paste0(
               seq(33, 90, length.out = 8) %>% round(),   # Start of each range
               "-", 
               seq(33, 90, length.out = 9)[-1] %>% round() # End of each range
             )
           )
         ) 

metadata_analysis <- metadata_analysis %>% mutate(
  Age_Category = cut(
    Age,
    breaks = seq(33, 90, length.out = 9),
    include.lowest = TRUE,
    labels = paste0(
      seq(33, 90, length.out = 8) %>% round(),   # Start of each range
      "-", 
      seq(33, 90, length.out = 9)[-1] %>% round() # End of each range
    )
  )
) 

## By Sample Matrix
for(i in unique(metadata_analysis_no_ex$SampleMatrix)){
  m <- metadata_analysis_no_ex %>% filter(SampleMatrix == i)
  
  # if(i == "Plasma"){targets = target_passed_detct_qc_p} else targets = target_passed_detct_qc_s
  if(i == "Plasma") {
    targets <- c(target_passed_detct_qc_p, target_ratios)
  } else {
    targets <- c(target_passed_detct_qc_s, target_ratios)
  }
  generate_heatmap_pca(figureDir = "figures/heatmap_pca/Updated_with_ratio/",
                       plot_prefix = paste0(i,"_Samples"),
                       data_matrix = project_all$data_wide_bridged_ratio[,m$SampleName_long],
                       metadata = m,
                       targets = targets,
                       sampleName = "SampleName_long",
                       annotate_sample_by = c("Group", "Randomisation", "Sex", "Age_Category"),
                       col_split_by = "Group",
                       pca_annotate_by = "Group",
                       pca_shape_by = "Sex")
}

m <- metadata_analysis %>% filter(SampleMatrix == "Serum")
generate_heatmap_pca(figureDir = "figures/heatmap_pca/Updated_with_ratio/",
                     plot_prefix = "Serum_Samples_exclude_outlier",
                     data_matrix = project_all$data_wide_bridged_ratio[,m$SampleName_long],
                     metadata = m,
                     targets = c(target_passed_detct_qc_s, target_ratios),
                     sampleName = "SampleName_long",
                     annotate_sample_by = c("Group", "Randomisation", "Sex", "Age_Category"),
                     col_split_by = "Group",
                     pca_annotate_by = "Group",
                     pca_shape_by = "Sex")

############################
## Venn Diagram
############################
# function to apply based on the result name
process_de_results <- function(name, df) {
  case <- sub("(.*)_Healthy control", "\\1", name)
  text <- paste0("Group", gsub("-", "_", case))
  return(df %>% filter(get(paste0(text, "_pval_FDR")) < 0.05) %>% pull(target))
}

## Pull DE targets
de_targets <- lapply(names(de_results), function(name) {
  process_de_results(name, de_results[[name]])
})
names(de_targets) <- sapply(names(de_results), function(name) gsub("_", " vs ", name))

# function to generate and save Venn diagrams
generate_venn_diagram <- function(venn_data, folder_prefix) {
  venn_plot <- ggVennDiagram(venn_data) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme(legend.position = "none") + 
    scale_x_continuous(expand = expansion(mult = .2))
  
  # Save the Venn diagram
  output_path <- paste0("figures/DE/Updated_Exclude_outlier_with_ratio/", folder_prefix, "_venn_diagram")
  ggsave(paste0(output_path, ".pdf"), plot = venn_plot, device = "pdf", width = 13, height = 8)
  ggsave(paste0(output_path, ".svg"), plot = venn_plot, device = "svg", width = 13, height = 8)
}

## Make venn diagram for all comparisons
generate_venn_diagram(venn_data = de_targets, folder_prefix = "All_comparisons")
generate_venn_diagram(venn_data = de_targets[-1], folder_prefix = "Serum_comparisons")

common_targets <- Reduce(intersect, de_targets)
Reduce(intersect, de_targets[[1:2]])

# common targets on venn diagram
venn_data <- ggVennDiagram::process_data(Venn(de_targets))
region_data <- venn_data$regionData
intersection_df <- data.frame(
  Intersection = region_data$name,  # Binary intersection codes
  Count = region_data$count,        # Number of items in each intersection
  stringsAsFactors = FALSE
)
intersection_df$Items <- sapply(region_data$item, paste, collapse = ", ")

###### Remake clustering heatmap with DE targets ######
for (i in seq_along(comparison$case)) {
  case <- comparison$case[i]
  sample_matrix <- comparison$sample_matrix[i]
  
  if(sample_matrix == "Plasma") {
    targets <- target_passed_detct_qc_p
  } else {
    targets <- target_passed_detct_qc_s
  }
  
  if(case == "MSA"){
    metadat <- metadata_analysis %>%
      mutate(Group = ifelse(grepl("MSA", Group), "MSA", Group)) %>%
      filter(Group %in% c(control, case), 
             SampleMatrix == sample_matrix) %>% 
      mutate(Group = factor(Group, levels = c(control, case)))
  } else {
    metadat <- metadata_analysis %>%
      filter(Group %in% c(control, case), 
             SampleMatrix == sample_matrix) %>% 
      mutate(Group = factor(Group, levels = c(control, case))) 
  }
  
  data_matrix <- project_all$data_wide_bridged[,metadat$SampleName_long]
  
  figureDir <- paste0("figures/DE/Updated_Exclude_outlier_with_ratio/", case, "_vs_", control, "_linear.with.cov/Group Effect/clustering/")
  plot_prefix_d <- "Group_pval_FDR"
  targets_d <- de_targets[[paste0(case, " vs ", control)]]
  
  pdf(file=file.path(figureDir, paste0(plot_prefix_d,'_Heatmap.pdf')), width = 10, height = 8)
  generate_heatmap(data = data_matrix,
                   sampleInfo = metadat,
                   sampleName_var = "SampleName_long",
                   target_subset = targets_d,
                   annotate_sample_by = c("Group", "Sex", "Age_Category"),
                   row_fontsize = 6,
                   cluster_column_slices = FALSE)
  dev.off()
  
  figureDir <- paste0("figures/DE/Updated_Exclude_outlier_with_ratio/", case, "_vs_", control, "_linear.with.cov/Group Effect/")
  plot_prefix <- "all.targets"
  
  pdf(file=file.path(figureDir, paste0(plot_prefix,'_Heatmap.pdf')), width = 10, height = 8)
  generate_heatmap(data = data_matrix,
                   sampleInfo = metadat,
                   sampleName_var = "SampleName_long",
                   target_subset = targets,
                   annotate_sample_by = c("Group", "Sex", "Age_Category"),
                   row_fontsize = 3,
                   cluster_column_slices = FALSE)
  dev.off()
}

###### Remake boxplot to show y-axis as NPQ bridged norm, also fix text and size ###### 
remake_boxplot <- function(case, control, suffix = "FDR",
                           cat1 = "Group", data_long){
  output_prefix <- paste0("figures/DE/Updated_Exclude_outlier_with_ratio/", case, "_vs_", control, "_linear.with.cov/Group Effect/")
  de <- read_csv(paste0(output_prefix, "csv/model.res.csv"))
  
  coefs <- paste0(cat1, gsub("-", "_",case), "_coef")
  p_vals <- paste0(cat1, gsub("-", "_",case), "_pval_", suffix)
  
  de_targets_up <- de %>% filter(!!sym(p_vals) < 0.05 & !!sym(coefs) > 0)
  de_targets_down <- de %>% filter(!!sym(p_vals) < 0.05 & !!sym(coefs) < 0)
  
  ## Boxplots
  plot_boxplots_by_cat <- function(data_to_plot, filename_prefix, title_text, w, l) {
    if (nrow(data_to_plot) > 0) {
      # Define plot parameters based on number of rows
      # plot_params <- define_plot_parameters(nrow(data_to_plot))
      
      pdf(paste0(output_prefix, filename_prefix, ".pdf"), width = w, height = l)
      ggplot(data_to_plot, aes(!!sym(cat1), NPQ_bridge_norm, fill = !!sym(cat1))) +
        geom_boxplot() +
        facet_wrap(~Target, scales = "free_y") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              strip.text = element_text(size = 13)) +
        labs(
          title = title_text,
          y = "NPQ Bridged Normalized"
        ) -> p
      
      print(p)
      dev.off()
      
      ggsave(filename = paste0(output_prefix, filename_prefix, ".svg"), device = "svg", width = w, height = l)
    }
  }
  
  plot_boxplots_by_cat(data_long %>% 
                         filter(Target %in% de_targets_up$target), 
                       paste0("boxplot/", suffix, "_Upregulated_Hits"), 
                       "Upregulated Targets", w=10, l=8)
  
  plot_boxplots_by_cat(data_long %>% 
                         filter(Target %in% de_targets_down$target), 
                       paste0("boxplot/", suffix, "_Downregulated_Hits"),
                       "Down regulated Targets", w=14, l=9)
  
}

for (i in seq_along(comparison$case)) {
  case <- comparison$case[i]
  sample_matrix <- comparison$sample_matrix[i]
  
  if(sample_matrix == "Plasma") {
    targets <- c(target_passed_detct_qc_p, target_ratios)
  } else {
    targets <- c(target_passed_detct_qc_s, target_ratios)
  }
  
  if(case == "MSA"){
    metadat <- metadata_analysis %>%
      mutate(Group = ifelse(grepl("MSA", Group), "MSA", Group)) %>%
      filter(Group %in% c(control, case),
             SampleMatrix == sample_matrix) %>% 
      mutate(Group = factor(Group, levels = c(control, case)))
    
    dat <- project_all$data_long_bridged_ratio %>%
      mutate(Group = ifelse(grepl("MSA", Group), "MSA", Group)) %>%
      filter(Target %in% targets,
             SampleName_long %in% metadat$SampleName_long) %>% 
      mutate(Group = factor(Group, levels = c(control, case))) 
    
  } else {
    metadat <- metadata_analysis %>%
      filter(Group %in% c(control, case),
             SampleMatrix == sample_matrix) %>% 
      mutate(Group = factor(Group, levels = c(control, case)))
    
    dat <- project_all$data_long_bridged_ratio %>%
      filter(Target %in% targets,
             SampleName_long %in% metadat$SampleName_long) %>% 
      mutate(Group = factor(Group, levels = c(control, case))) 
  }
  
  metadat %>% filter(Group == case) %>% nrow() -> ncase
  metadat %>% filter(Group == control) %>% nrow() -> ncontrol
  
  cat(control, "=", ncontrol, ",", case, "=", ncase, "\n")
  
  remake_boxplot(case = case, control = control, 
                 data_long = dat)
}


