sort_gene <- function(dedf) {
  sortList <- dedf$logFC
  names(sortList) <- dedf$gene
  sortList <- sort(sortList, decreasing = T)
  return(sortList)
}
