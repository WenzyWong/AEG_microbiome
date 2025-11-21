sort_gene <- function(dedf) {
  sortList <- dedf$log2FC
  names(sortList) <- dedf$symbol
  sortList <- sort(sortList, decreasing = T)
  return(sortList)
}
