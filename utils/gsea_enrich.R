gsea_enrich <- function(level, geneList, species) {
    set.seed(129)
    library(clusterProfiler)
    library(pathview)
    library(dplyr)
    geneList <- sort(geneList, decreasing = T)
    if (species == "Mm") {gmt <- read.gmt(paste0("/data/yzwang/reference/gsea/", level, ".v2025.1.Mm.symbols.gmt"))}
    if (species == "Hs") {gmt <- read.gmt(paste0("/data/yzwang/reference/GSEA/", level, ".v2023.1.Hs.symbols.gmt"))}
    gsea <- GSEA(geneList = geneList, TERM2GENE = gmt,
                 minGSSize = 1, verbose = F, pvalueCutoff = 0.05)
    return(gsea)
}
