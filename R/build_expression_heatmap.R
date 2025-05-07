#' Build a heatmap of normalized RNA expression for selected proteins and cell types
#'
#' This function subsets the RNA-seq reference data by user-defined proteins and immune cell types,
#' performs min-max normalization of expression values (pTPM), and visualizes the result as a heatmap.
#'
#' @param data A data.frame as returned by `load_reference()`, containing RNA expression data.
#' @param proteins A character vector of gene/protein names to include (must match `Gene.name` column).
#' @param cells_of_interest A character vector of immune cell types to include (must match `Immune.cell` column).
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{data}{A normalized, wide-format data.frame with scaled pTPM values.}
#'   \item{plot}{A ggplot2 heatmap object.}
#' }
#' @export
build_expression_heatmap <- function(data, proteins, cells_of_interest) {
  # Check required columns
  required_cols <- c("Immune.cell", "Gene.name", "pTPM")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Subset and reshape to wide format
  expression_wide <- data %>%
    dplyr::select(-Gene, -TPM, -nTPM, dplyr::everything()) %>%
    dplyr::filter(
      Immune.cell %in% cells_of_interest,
      Gene.name %in% proteins
    ) %>%
    tidyr::pivot_wider(names_from = Gene.name, values_from = pTPM)

  # Min-max normalize each numeric column
  normalize <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(NA, length(x)))
    (x - rng[1]) / diff(rng)
  }

  scaled <- expression_wide
  scaled[,-1] <- as.data.frame(apply(expression_wide[,-1], 2, normalize))

  # Drop protein columns with all NA values
  protein_data <- scaled[, -1]
  keep_cols <- !apply(protein_data, 2, function(x) all(is.na(x)))
  scaled_clean <- cbind(Immune.cell = scaled[[1]], protein_data[, keep_cols, drop = FALSE])

  # Reshape for plotting
  long_df <- scaled_clean %>%
    tidyr::pivot_longer(-Immune.cell, names_to = "protein", values_to = "value")

  # Generate heatmap
  heatmap <- ggplot2::ggplot(long_df, ggplot2::aes(x = Immune.cell, y = protein, fill = value)) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_gradient(low = "white", high = "blue", na.value = "grey80") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = NULL, fill = "Scaled pTPM") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(list(data = scaled_clean, plot = heatmap))
}
