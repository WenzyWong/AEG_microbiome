wilcox_diff <- function(dt, treat, control, ifpair, P.cutoff, LFC.cutoff){
  library(dplyr)
  Pvalue <- c(rep(0, nrow(dt)))
  log2FC <- c(rep(0, nrow(dt)))
  for(i in 1:nrow(dt)){
    if(sd(treat[i, ]) == 0 | sd(control[i, ]) == 0){
      Pvalue[i] <- 1
      log2FC[i] <- 0
    } else{
      y <- wilcox.test(as.numeric(treat[i, ]),
                       as.numeric(control[i, ]),
                       paired = ifpair, exact = F)
      Pvalue[i] <- y$p.value
      log2FC[i] <- log2((mean(as.numeric(treat[i, ])) + 1)/
                          (mean(as.numeric(control[i, ])) + 1))
    }
  }
  diff <- data.frame(ID = rownames(dt),
                     log2FC = as.numeric(log2FC),
                     Pvalue = as.numeric(Pvalue))
  diff$P.adj <- p.adjust(diff$Pvalue, "BH")
  diff$Change <- case_when(diff$log2FC >= LFC.cutoff & diff$P.adj <= P.cutoff ~ "Up",
                           diff$log2FC <= -LFC.cutoff & diff$P.adj <= P.cutoff ~ "Dn",
                           abs(diff$log2FC) < LFC.cutoff | diff$P.adj > P.cutoff ~ "NS")
  return(diff)
}
