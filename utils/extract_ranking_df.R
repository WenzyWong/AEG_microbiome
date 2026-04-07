extract_ranking_df <- function(gsea_list, gene_set = c("hallmark", "kegg"),
                               padj_thresh = 0.05) {
  gene_set <- match.arg(gene_set)
  
  bind_rows(lapply(names(gsea_list), function(sp) {
    obj <- gsea_list[[sp]][[gene_set]]
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    
    res       <- obj@result
    gene_rank <- obj@geneList
    sig_ids   <- res %>% filter(p.adjust < padj_thresh) %>% pull(ID)
    if (length(sig_ids) == 0) return(NULL)
    
    bind_rows(lapply(sig_ids, function(pid) {
      row_idx    <- which(res$ID == pid)
      if (length(row_idx) != 1) return(NULL)
      core_str   <- res$core_enrichment[row_idx]
      if (is.na(core_str) || nchar(core_str) == 0) return(NULL)
      matched    <- intersect(strsplit(core_str, "/")[[1]], names(gene_rank))
      if (length(matched) == 0) return(NULL)
      
      data.frame(pathway = pid, species = sp,
                 NES = res$NES[row_idx], stat = gene_rank[matched],
                 stringsAsFactors = FALSE)
    }))
  }))
}