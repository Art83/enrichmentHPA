#' Prepare HPA-style reference data for enrich_by_region()
#'
#' This helper standardises HPA expression data and optionally drops
#' genes that are zero across all groups.
#'
#' @param df_ref Data frame in long format with at least:
#'   - a gene column (default "Gene.name")
#'   - a grouping column (e.g. "Region", "Immune.cell")
#'   - an expression column (e.g. "pTPM")
#' @param gene_col Column containing gene identifiers.
#' @param group_col Column containing group labels (regions/cell types).
#' @param expr_col Column containing numeric expression values.
#' @param drop_all_zero Logical; if TRUE (default), drop genes that are zero
#'   (or NA) across all groups.
#'
#' @return A cleaned data.frame with the same columns as input (plus any extras),
#'   ready to be passed to `enrich_by_region()`.
#' @export
prep_hpa_atlas <- function(df_ref,
                             gene_col  = Gene.name,
                             group_col = Region,
                             expr_col  = pTPM,
                             drop_all_zero = TRUE) {
    stopifnot(is.data.frame(df_ref))

    gene_col_sym  <- rlang::ensym(gene_col)
    group_col_sym <- rlang::ensym(group_col)
    expr_col_sym  <- rlang::ensym(expr_col)

    # Standardise core columns for the operations
    core <- df_ref |>
      dplyr::mutate(
        .gene  = as.character(!!gene_col_sym),
        .group = as.character(!!group_col_sym),
        .expr  = suppressWarnings(as.numeric(!!expr_col_sym))
      ) |>
      dplyr::select(.gene, .group, .expr)
    if (drop_all_zero) {
      # Identify genes that are zero across all groups
      gene_zero_flag <- core |>
        dplyr::group_by(.gene) |>
        dplyr::summarise(
          all_zero = all((.expr) == 0 | is.na(.expr)),
          .groups  = "drop"
        )

      keep_genes <- gene_zero_flag$.gene[!gene_zero_flag$all_zero]

      core <- core |>
        dplyr::filter(.gene %in% keep_genes)
    }

    core
}
