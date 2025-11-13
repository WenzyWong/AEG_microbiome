#' Calculate hazard ratios from Cox proportional hazards models
#'
#' @param variable A character vector of variable names to test in Cox models
#' @param hrInput A data frame containing survival data with columns 'time' and 'state',
#'   plus the variables specified in the 'variable' parameter
#' @param cut.up Upper threshold for hazard ratio to classify as "Increase" (default threshold)
#' @param cut.dn Lower threshold for hazard ratio to classify as "Decrease" (default threshold)
#' @param cut.P P-value threshold for statistical significance (default threshold)
#'
#' @return A data frame with the following columns:
#'   - HR: Hazard ratio (exp(coef))
#'   - HR.95L: Lower bound of 95% confidence interval
#'   - HR.95H: Upper bound of 95% confidence interval
#'   - Pvalue: Wald test p-value (rounded to 2 significant digits)
#'   - Risk: Classification as "Increase", "Decrease", or "NS" (not significant)
#'     based on HR thresholds and p-value cutoff
#'
#' @details
#' For each variable, the function:
#' 1. Fits a univariate Cox proportional hazards model: Surv(time, state) ~ variable
#' 2. Extracts hazard ratio, 95% CI, and p-value
#' 3. Classifies risk as:
#'    - "Increase": HR > cut.up AND Pvalue < cut.P
#'    - "Decrease": HR < cut.dn AND Pvalue < cut.P
#'    - "NS": Otherwise (not significant or HR within thresholds)
#'
#' @examples
#' # Calculate HRs for multiple genes
#' hr_results <- hr_calc(
#'   variable = c("gene1", "gene2", "gene3"),
#'   hrInput = survival_data,
#'   cut.up = 1.5,
#'   cut.dn = 0.67,
#'   cut.P = 0.05
#' )
#'
#' @importFrom survival coxph Surv
#' @importFrom dplyr case_when
#' @export
hr_calc <- function(variable, input, cut.up, cut.dn, cut.P) {
  library(survival)
  library(survminer)
  library(dplyr)
  formulas <- sapply(variable,
                     function(x) as.formula(paste("Surv(time, state)~", x)))

  models <- lapply(formulas, function(x) {coxph(x, data = input)})

  res <- lapply(models,
                function(x){
                  x <- summary(x)
                  p.value <- signif(x$wald["pvalue"], digits = 2)
                  HR <- x$conf.int[, "exp(coef)"]
                  HR.95L <- signif(x$conf.int[,"lower .95"], 3)
                  HR.95H <- signif(x$conf.int[,"upper .95"], 3)
                  res <- c(HR, HR.95L, HR.95H, p.value)
                  names(res) <- c("HR", "HR.95L", "HR.95H", "Pvalue")
                  return(res)
                })
  res <- as.data.frame(t(as.data.frame(res)))
  res$Risk <- case_when(
    res$HR < cut.dn & res$Pvalue < cut.P ~ "Decrease",
    res$HR > cut.up & res$Pvalue < cut.P ~ "Increase",
    res$Pvalue >= cut.P | (res$HR >= cut.dn & res$HR <= cut.up) ~ "NS"
  )

  return(res)
}
