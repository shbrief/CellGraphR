#' Test factor associations with sample groups
#'
#' Tests whether MOFA factor scores are significantly associated with different
#' sample groups using the Kruskal-Wallis test. Grouping variables can be
#' specified explicitly or selected automatically from the sample metadata.
#'
#' @importFrom stats kruskal.test p.adjust
#'
#' @param object A trained MOFA2 object containing factor scores and sample
#'   metadata.
#' @param factor Character or integer(1). The MOFA factor to test (required).
#' @param groups Character vector. Names of metadata columns to use as grouping
#'   variables. If \code{NULL} (default), the function automatically selects
#'   columns with 2–10 unique non-\code{NA} values.
#'
#' @return A \code{data.frame} with one row per tested grouping variable and
#'   columns:
#'   \describe{
#'     \item{group}{Name of the grouping variable.}
#'     \item{statistic}{Kruskal-Wallis H statistic.}
#'     \item{p_value}{Raw p-value.}
#'     \item{df}{Degrees of freedom.}
#'     \item{p_adjusted}{Benjamini-Hochberg adjusted p-value.}
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts factor scores for the specified factor from the MOFA object.
#'   \item Merges factor scores with sample metadata.
#'   \item Selects grouping variables (auto or user-specified).
#'   \item Runs Kruskal-Wallis tests and applies Benjamini-Hochberg correction.
#' }
#'
#' When \code{groups = NULL}, columns with 2–10 unique non-\code{NA} values are
#' selected (these typically represent categorical clinical variables such as
#' Gleason grade, tumour stage, or treatment arm).
#'
#' @examples
#' \dontrun{
#' # Auto-select grouping variables
#' results <- testFactorAssociations(mofa_object, factor = 1)
#'
#' # Test specific grouping variables
#' results <- testFactorAssociations(mofa_object,
#'                                   factor = "Factor1",
#'                                   groups = c("gleason_grade", "pathologic_stage"))
#'
#' significant <- results[results$p_adjusted < 0.05, ]
#' }
#'
#' @seealso \code{\link[stats]{kruskal.test}}, \code{\link[stats]{p.adjust}}
#'
#' @export
testFactorAssociations <- function(object,
                                   factor = NULL,
                                   groups = NULL) {

    if (is.null(factor)) {
        stop("Specify the MOFA Factor to run the test.")
    }

    ## Extract factor scores (Z matrix)
    factor_scores <- get_factors(object,
                                 factors    = factor,
                                 groups     = "all",
                                 as.data.frame = TRUE)

    ## Get sample metadata
    metadata <- samples_metadata(object)

    ## Merge factor scores with metadata
    df_factors <- merge(factor_scores, metadata, by = "sample")

    ## Auto-select groups when not specified (fix: null check before validity check)
    if (is.null(groups)) {
        n_unique <- apply(metadata, 2L, function(x) length(unique(stats::na.omit(x))))
        groups   <- names(n_unique[n_unique >= 2L & n_unique <= 10L])
    } else {
        ## Validate user-specified groups
        bad <- groups[!groups %in% colnames(metadata)]
        if (length(bad) > 0L) {
            stop("The following group(s) are not in the metadata: ",
                 paste(bad, collapse = ", "))
        }
    }

    if (length(groups) == 0L) {
        warning("No grouping variables selected. Check metadata columns.")
        return(data.frame(group = character(), statistic = numeric(),
                          p_value = numeric(), df = integer(),
                          p_adjusted = numeric()))
    }

    ## Run Kruskal-Wallis test for each group
    kw_results <- lapply(groups, function(g) {
        formula_str <- paste("value ~", g)
        kw <- stats::kruskal.test(stats::as.formula(formula_str), data = df_factors)
        data.frame(group     = g,
                   statistic = kw$statistic,
                   p_value   = kw$p.value,
                   df        = kw$parameter,
                   stringsAsFactors = FALSE)
    })

    result_df <- do.call(rbind, kw_results)
    result_df$p_adjusted <- stats::p.adjust(result_df$p_value, method = "BH")
    rownames(result_df)  <- NULL

    result_df
}
