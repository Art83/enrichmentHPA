#' Perform enrichment analysis across regions or cell types
#'
#' Applies enrichment testing (GSEA or hypergeometric) across multiple groups (e.g., brain regions or immune cell types).
#'
#' @param gene_list Character vector of gene symbols.
#' @param df_ref Expression reference table returned by `load_reference()`.
#' @param group_col Name of the grouping column (e.g. "Region", "Immune.cell").
#' @param expression_col Column holding expression values (e.g. "pTPM").
#' @param method "gsea" or "hypergeometric".
#' @param threshold Numeric; pTPM cutoff for defining enrichment for hypergeometric test.
#' @param universe_type "global" or "group"; defines if the universe is shared or per-group.
#' @param n_perm Number of permutations for GSEA p-value (default = 100).
#' @param gsea_weight GSEA running sum exponent (default = 1).
#' @param use_parallel Use parallel computation for permutations (GSEA only).
#' @param seed Optional random seed.
#'
#' @return A data.frame with enrichment results per region/cell type.
#' @export
enrich_by_region <- function(gene_list,
                             df_ref,
                             group_col = "Region",
                             expression_col = "pTPM",
                             method = c("gsea", "hypergeometric"),
                             threshold = 10,
                             universe_type = c("global", "group"),
                             n_perm = 100,
                             gsea_weight = 1,
                             use_parallel = FALSE,
                             seed = NULL) {
  method <- match.arg(method)
  universe_type <- match.arg(universe_type)

  groups <- unique(df_ref[[group_col]])
  results <- list()

  # Define global universe if needed
  global_universe <- unique(df_ref$Gene.name)

  for (grp in groups) {
    df_grp <- dplyr::filter(df_ref, .data[[group_col]] == grp)

    # Define universe
    universe <- if (universe_type == "global") global_universe else unique(df_grp$Gene.name)

    if (method == "hypergeometric") {
      mean_expr <- df_grp %>%
        dplyr::group_by(Gene.name) %>%
        dplyr::summarise(mean_expr = mean(.data[[expression_col]], na.rm = TRUE), .groups = "drop")

      enriched_genes <- mean_expr$Gene.name[mean_expr$mean_expr > threshold]

      res <- run_enrichment(
        gene_list = gene_list,
        universe = universe,
        enriched_genes = enriched_genes,
        method = "hypergeometric"
      )
    }

    if (method == "gsea") {
      gene_stats <- tapply(df_grp[[expression_col]], df_grp$Gene.name, mean, na.rm = TRUE)

      res <- run_enrichment(
        gene_list = gene_list,
        universe = universe,
        gene_stats = gene_stats,
        method = "gsea",
        gsea_weight = gsea_weight,
        n_perm = n_perm,
        seed = seed
      )
    }

    res$group <- grp
    results[[grp]] <- res
  }

  dplyr::bind_rows(results, .id = "group_id")
}
