#' Subset HPA RNA-seq reference data by gene and group
#'
#' This function filters a Human Protein Atlas (HPA) RNA-seq reference dataset by gene symbols
#' and group values (e.g., immune cell types, brain regions, or cluster IDs).
#' It supports flexible column selection to adapt across different dataset types.
#'
#' @param data A data.frame returned by `load_reference()`.
#' @param genes Character vector of gene names (matching `Gene.name` column).
#' @param groups Character vector of group values (e.g., cell types, regions).
#' @param group_col Name of the grouping column (e.g. `Immune.cell`).
#'
#' @return A filtered data.frame containing only matching genes and groups.
#' @export
subset_reference <- function(data,
                             genes,
                             groups,
                             group_col,
                             ...) {
  # Input validation
  if (!"Gene.name" %in% colnames(data)) {
    stop("The input data must contain a 'Gene.name' column.")
  }
  if (!group_col %in% colnames(data)) {
    stop("Column '", group_col, "' not found in the input data.")
  }
  dots <- rlang::enquos(...)

  named_dots <- dots[nzchar(names(dots))]
  unnamed_dots <- dots[!nzchar(names(dots))]

  # Convert named args into proper == expressions
  named_filters <- purrr::imap(named_dots, function(val, name) {
    rlang::new_quosure(
      rlang::expr(.data[[!!name]] == !!val),
      env = rlang::caller_env()
    )
  })

  # Combine and forcibly remove all names to avoid dplyr::filter() complaint
  all_filters <- unname(c(unnamed_dots, named_filters))


  filtered <- dplyr::filter(
    data,
    Gene.name %in% genes,
    .data[[group_col]] %in% groups,
    !!!all_filters
  )

  return(filtered)
}

