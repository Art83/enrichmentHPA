#' Summarize expression levels and categorize by strength
#'
#' Computes mean pTPM and expression category for each gene across groups.
#'
#' @param data A data.frame with HPA expression values (pTPM).
#' @param group_col Column to group by (e.g. "Immune.cell", "Region").
#' @param gene_col Column for gene names (default: "Gene.name").
#' @param expr_col Column for expression values (default: "pTPM").
#' @param thresholds Numeric vector of length 2: c(medium_cutoff, high_cutoff)
#'
#' @return A data.frame with gene stats and expression category.
#' @export
summarize_expression <- function(data,
                                 group_col,
                                 gene_col = "Gene.name",
                                 expr_col = "pTPM",
                                 thresholds = c(1, 10)) {
  # Defensive checks
  if (!all(c(group_col, gene_col, expr_col) %in% colnames(data))) {
    stop("Missing one or more required columns.")
  }

  data_summary <- data %>%
    dplyr::group_by(.data[[gene_col]]) %>%
    dplyr::summarise(
      mean_expr = mean(.data[[expr_col]], na.rm = TRUE),
      sd_expr   = sd(.data[[expr_col]], na.rm = TRUE),
      cv_expr   = sd_expr / mean_expr,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      expr_category = dplyr::case_when(
        mean_expr <= thresholds[1]        ~ "low",
        mean_expr <= thresholds[2]        ~ "medium",
        mean_expr >  thresholds[2]        ~ "high",
        TRUE                              ~ NA_character_
      )
    )

  return(data_summary)
}
