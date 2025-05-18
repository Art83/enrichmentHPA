#' Subset HPA RNA-seq reference data by gene and group
#'
#' This function filters a Human Protein Atlas (HPA) RNA-seq reference dataset by gene symbols
#' and group values (e.g., immune cell types, brain regions, or cluster IDs).
#' It supports flexible column selection to adapt across different dataset types.
#'
#' @param data A data.frame returned by `load_reference()`.
#' @return A filtered data.frame containing only matching genes and groups.
#' @export
subset_reference <- function(data, ...) {
  if (!"Gene.name" %in% colnames(data)) {
    stop("The input data must contain a 'Gene.name' column.")
  }

  dots <- rlang::enquos(...)
  named_dots <- dots[nzchar(names(dots))]
  unnamed_dots <- dots[!nzchar(names(dots))]

  named_filters <- purrr::imap(named_dots, function(val, name) {
    rlang::new_quosure(
      rlang::expr(.data[[!!name]] %in% !!val),
      env = rlang::caller_env()
    )
  })

  all_filters <- unname(c(unnamed_dots, named_filters))

  dplyr::filter(data, !!!all_filters)
}

