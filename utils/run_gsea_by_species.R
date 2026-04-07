run_gsea_by_species <- function(species_name, de_results, dir_rds) {
  DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"
  source(file.path(DIR_TOOL, "sort_gene.R"))
  source(file.path(DIR_TOOL, "gsea_enrich.R"))
  
  de_df <- de_results[[species_name]]$de_df
  if (is.null(de_df) || nrow(de_df) == 0) {
    warning(paste("No DE result for:", species_name))
    return(NULL)
  }
  
  gene_list <- sort_gene(de_df)
  safe_name <- gsub("[^A-Za-z0-9_]", "_", species_name)
  
  run_one <- function(level, label) {
    res <- tryCatch(
      gsea_enrich(level, gene_list, "Hs"),
      error = function(e) {
        warning(paste(label, "GSEA failed for:", species_name, "-", e$message))
        NULL
      }
    )
    if (!is.null(res))
      saveRDS(res, file.path(dir_rds, paste0("gsea_", label, "_", safe_name, ".rds")))
    res
  }
  
  list(
    hallmark = run_one("h.all",        "hallmark"),
    kegg     = run_one("c2.cp.kegg",   "kegg")
  )
}