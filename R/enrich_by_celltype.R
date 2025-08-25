#' Cell-type-specific enrichment from single-cell reference
#'
#' @param gene_list A vector of genes of interest (e.g. DEGs)
#' @param expr_df A long-form dataframe with columns: Gene, CellType, Expression
#' @param threshold Minimum expression to consider a gene expressed (default = 1 on log-scale)
#' @param method Statistical method: "hypergeometric" or "gsea"
#' @param universe Optional background gene universe
#' @return A data frame with enrichment results per cell type
#' @export
enrich_by_celltype <- function(gene_list,
                               expr_df,
                               threshold = 1,
                               method = c("hypergeometric", "gsea"),
                               universe = NULL) {
  method <- match.arg(method)

  expr_df <- expr_df %>% dplyr::filter(!is.na(Expression))

  results <- expr_df %>%
    dplyr::group_by(CellType) %>%
    dplyr::group_map(~{
      cell_data <- .x
      cell_type <- unique(cell_data$CellType)

      expressed_genes <- cell_data %>%
        dplyr::filter(Expression >= threshold) %>%
        dplyr::pull(Gene) %>%
        unique()

      if (length(expressed_genes) == 0) return(NULL)

      gene_stats <- NULL
      if (method == "gsea") {
        gene_stats <- cell_data %>%
          dplyr::distinct(Gene, Expression) %>%
          dplyr::arrange(desc(Expression)) %>%
          tibble::deframe()
      }

      enrichmentHPA::run_enrichment(
        gene_list = gene_list,
        enriched_genes = expressed_genes,
        gene_stats = gene_stats,
        method = method,
        universe = universe
      ) %>%
        dplyr::mutate(CellType = cell_type)
    }) %>%
    dplyr::bind_rows()

  return(results)
}
