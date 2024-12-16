#####################################################
# AEG microbiome
# 2024-12-14
# Yunzhe Wang
#
# Step 06.Analysing metabolome from EF supernatant
#
#####################################################

setwd("/data/yzwang/project/AEG_seiri/")

library(readxl) # Read excel
library(limma)  # Differential analysis
library(ropls)  # PCA & PLS-DA # BiocManager::install("ropls")
library(dplyr)
library(stats)  # Basic statistical test
library(ggplot2)
library(pheatmap)  # Heatmap visualization

DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TABLE <- "/data/yzwang/project/AEG_seiri/table_infos/metabolome/"

info <- as.data.frame(read_excel(paste0(DIR_TABLE, "intensity_information_en.xlsx")))
rownames(info) <- info$Compound_ID
nrow(info) # 3445
length(grep("C", info$KEGG_ID)) # 992
length(grep("H", info$HMDB_ID)) # 2012
length(info$CAS[info$CAS != "-"]) # 1572

#####################################################
# Function defination
filter_by_cv <- function(data_matrix, group_labels, 
                         cv_threshold = 0.5, min_samples = 3) {
  # # group_labels is a factor that identifies which group each sample belongs to
  group_labels <- as.factor(group_labels)
  
  # Initialize the result list
  filtered_results <- list()
  cv_stats <- list()
  
  # Calculate CV for each group
  for(group in levels(group_labels)) {
    group_data <- data_matrix[, group_labels == group]
    
    # Calculate the mean and standard deviation of each metabolite
    means <- rowMeans(group_data, na.rm = TRUE)
    sds <- apply(group_data, 1, sd, na.rm = TRUE)
    
    # Calculate the CV
    cvs <- sds / means
    
    # Store the CV results for each group
    cv_stats[[group]] <- cvs
  }
  
  cv_df <- do.call(cbind, cv_stats)
  colnames(cv_df) <- paste0("CV_", levels(group_labels))
  
  # Calculate the number of times each metabolite exceeded the threshold in each group
  cv_exceeds <- rowSums(cv_df > cv_threshold, na.rm = TRUE)
  
  # Retain the metabolites with CV below threshold in all groups
  keep_metabolites <- cv_exceeds == 0
  
  # Filter the data matrix
  filtered_matrix <- data_matrix[keep_metabolites, ]
  
  # Prepare statistics
  stats <- data.frame(
    Original_Features = nrow(data_matrix),
    Filtered_Features = nrow(filtered_matrix),
    Removed_Features = nrow(data_matrix) - nrow(filtered_matrix),
    CV_Threshold = cv_threshold
  )
  
  # Prepare detailed information for each metabolite
  metabolite_info <- data.frame(
    cv_df,
    Kept = keep_metabolites,
    row.names = rownames(data_matrix)
  )
  
  filter_results <- list(
    filtered_matrix = filtered_matrix,
    stats = stats,
    metabolite_info = metabolite_info
  )
  
  # Sort CV values into long format
  cv_data <- stack(filter_results$metabolite_info[, grep("^CV_", colnames(filter_results$metabolite_info))])
  
  # Create a boxplot
  boxplot(values ~ ind, data = cv_data,
          main = "CV Distribution Across Groups",
          ylab = "Coefficient of Variation (CV)",
          xlab = "Groups",
          outline = TRUE)
  
  # Add the CV threshold line
  abline(h = filter_results$stats$CV_Threshold, 
         col = "red", 
         lty = 2)
  legend("topright", 
         legend = paste("CV threshold =", filter_results$stats$CV_Threshold),
         col = "red", 
         lty = 2)
  
  # Return results
  return(filter_results)
}

# Differential analysis
analyze_metabolomics <- function(data_matrix, group_labels, 
                                 adj_method = "BH", 
                                 fc_threshold = 1.3,
                                 pval_threshold = 0.05) {
  # The data matrix should be sample in column, metabolite in row
  # group_labels is a factor that identifies which group each sample belongs to
  
  # 1. Performed differential analysis using limma
  data_matrix <- log2(data_matrix + 1)
  design <- model.matrix(~0 + group_labels)
  colnames(design) <- levels(group_labels)
  fit <- lmFit(data_matrix, design)
  
  # Set contrasts
  contrast.str <- paste0(levels(group_labels)[2], "-", levels(group_labels)[1])
  contrast_matrix <- makeContrasts(contrasts = contrast.str, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Get the results of the variance analysis
  results <- topTable(fit2, adjust.method = adj_method, number = Inf)
  
  # 2. Add fold change column
  results$FC <- 2^results$logFC
  
  # 3. Marker for significant differences in metabolites
  results$Significant <- with(results, 
                              adj.P.Val <= pval_threshold & abs(logFC) >= log2(fc_threshold)
  )
  results$Change <- with(results,
                         case_when(
                           adj.P.Val <= pval_threshold & logFC >= log2(fc_threshold) ~ "Up",
                           adj.P.Val <= pval_threshold & logFC <= -log2(fc_threshold) ~ "Dn",
                           TRUE ~ "NS"
                         ))
  
  # 4. Perform a univariate t-test as a supplement
  t_test_results <- apply(data_matrix, 1, function(x) {
    test <- t.test(x[group_labels == levels(group_labels)[1]], 
                   x[group_labels == levels(group_labels)[2]])
    c(t_test_pval = test$p.value)
  })
  results$t_test_pval <- t_test_results
  
  # 5. Multivariable analysis
  # PCA
  pca_model <- opls(t(data_matrix), group_labels, predI = 2)
  
  # PLS-DA
  plsda_model <- opls(t(data_matrix), group_labels, predI = 1, orthoI = NA)
  
  # Return result list
  return(list(
    diff_results = results,
    pca_model = pca_model,
    plsda_model = plsda_model
  ))
}

# Visuliazation
plot_metabolomics_results <- function(analysis_results, data_matrix, 
                                      group_labels,
                                      fc_threshold = 1.3,
                                      pval_threshold = 0.05) {
  # 1. Heat map-shows only significantly different metabolites
  sig_metabolites <- rownames(analysis_results$diff_results)[
    analysis_results$diff_results$Significant
  ]
  if(length(sig_metabolites) > 0) {
    sig_data <- data_matrix[sig_metabolites,]
    pheatmap(sig_data,
             annotation_col = data.frame(
               Group = group_labels,
               row.names = colnames(data_matrix)
             ),
             scale = "row",
             main = "Significant Metabolites Heatmap")
  }
  
  logFC_max <- max(abs(analysis_results$diff_results$logFC))
  draw_label <- analysis_results$diff_results[analysis_results$diff_results$Significant, ]
  draw_label <- draw_label[order(draw_label$logFC), ]
  draw_label <- draw_label[c(1:3, (nrow(draw_label)-2):nrow(draw_label)), ]
  draw_label$Name <- info[rownames(draw_label), "Name"]
  
  # 2. Volcano plot
  ggplot(analysis_results$diff_results, 
         aes(x = logFC, y = -log10(adj.P.Val), colour = Change)) + 
    ggtitle("Treatment v.s. Control") + 
    geom_point(size = 1, alpha = .8) + 
    xlim(-logFC_max, logFC_max) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    geom_vline(xintercept = log2(fc_threshold), linetype = "dashed") +
    geom_vline(xintercept = -log2(fc_threshold), linetype = "dashed") +
    geom_point(draw_label, mapping = aes(x = logFC, 
                                         y = -log10(adj.P.Val),
                                         color = Change), 
               size = 3) + 
    scale_color_manual(values = c("#126CAA", "grey", "#9A342C")) +
    ggrepel::geom_label_repel(data = draw_label,
                              aes(x = logFC, y = -log10(adj.P.Val), 
                                  label = Name),
                              color="grey27",
                              #nudge_x = 1,
                              #nudge_y = -1,
                              box.padding = 1.5,
                              alpha = .8) +
    theme_test() +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1)) 
}

#####################################################
# Perform analysis
mtx_all <- as.matrix(info[ , (ncol(info) - 6):ncol(info)])
colnames(mtx_all)
rownames(mtx_all) <- info$Compound_ID

group_labels <- factor(c(rep("Treatment", 4), rep("Control", 3)),
                       levels = c("Control", "Treatment"))

# Filter by: 1. whether from human; 2. CV
info_rmhm <- info[info$HMDB_ID == "-", ]
mtx_rmhm <- mtx_all[rownames(mtx_all) %in% info_rmhm$Compound_ID, ]
dim(mtx_rmhm)
res_rmhm <- filter_by_cv(mtx_rmhm, group_labels)
dim(res_rmhm$filtered_matrix)

result_rmhm <- analyze_metabolomics(res_rmhm$filtered_matrix, group_labels)
plot_metabolomics_results(result_rmhm, res_rmhm$filtered_matrix, group_labels)

write.csv(result_rmhm$diff_results, "./results/EF_metab/Diff_removeHMDB.csv")

sig_rmhm <- result_rmhm$diff_results[result_rmhm$diff_results$Significant, ]
dim(sig_rmhm)
sig_rmhm$Compound_ID <- rownames(sig_rmhm)
sig_rmhm <- merge(sig_rmhm, info_rmhm, by = "Compound_ID", all.x = TRUE)
sig_rmhm <- sig_rmhm[order(sig_rmhm$logFC, decreasing = T), ]
write.csv(sig_rmhm, "./results/EF_metab/DiffSig_removeHMDB.csv")

# Filter by CV only
res_filt <- filter_by_cv(mtx_all, group_labels)
dim(res_filt$filtered_matrix)
print(res_filt$stats)
head(res_filt$metabolite_info)
mtx_filt <- res_filt$filtered_matrix

results <- analyze_metabolomics(mtx_filt, group_labels)
plot_metabolomics_results(results, mtx_filt, group_labels)

write.csv(results$diff_results, "./results/EF_metab/Diff_all.csv")

sig <- results$diff_results[results$diff_results$Significant, ]
dim(sig)
sig$Compound_ID <- rownames(sig)
sig <- merge(sig, info, by = "Compound_ID", all.x = TRUE)
sig <- sig[order(sig$logFC, decreasing = T), ]
write.csv(sig, "./results/EF_metab/DiffSig_all.csv")
