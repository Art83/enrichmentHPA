#' Enrichment Analysis by Tissue
#'
#' Runs GSEA or hypergeometric enrichment across tissues using a reference expression dataset.
#'
#' @param gene_list Character vector of genes of interest.
#' @param reference_df Data frame containing columns: Gene.name, Tissue, and expression values.
#' @param method Enrichment method: "gsea" or "hypergeo".
#' @param expr_col Column name for expression values (e.g., "nTPM").
#' @param universe Optional character vector of background genes. If NULL, all genes in reference are used.
#' @param expression_threshold Expression cutoff for defining "expressed" genes.
#' @param n_perm Number of permutations (only for GSEA).
#' @param min_genes Minimum number of expressed genes per tissue to include in analysis.
#' @param verbose If TRUE, prints progress messages.
#'
#' @return A data.frame with enrichment results per tissue.
#' @export
enrich_by_tissue <- function(gene_list,
                             reference_df,
                             method = c("gsea", "hypergeo"),
                             expr_col = "nTPM",
                             universe = NULL,
                             expression_threshold = 1,
                             n_perm = 1000,
                             min_genes = 10,
                             verbose = FALSE) {
  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }
  method <- match.arg(method)

  if (!all(c("Gene.name", "Tissue", expr_col) %in% colnames(reference_df))) {
    stop("Reference must include 'Gene.name', 'Tissue', and selected expression column.")
  }

  reference_df <- dplyr::filter(reference_df, !is.na(.data[[expr_col]]))
  if(is.null(universe)){
    universe <- unique(reference_df$Gene.name)
  }

  gene_list <- intersect(gene_list, universe)

  results <- reference_df %>%
    dplyr::group_by(Tissue) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(df_tissue) {
      tissue <- unique(df_tissue$Tissue)

      df_tissue <- dplyr::filter(df_tissue, Gene.name %in% universe)
      if (nrow(df_tissue) < min_genes) return(NULL)

      if (method == "gsea") {
        ranked_genes <- df_tissue %>%
          dplyr::distinct(Gene.name, .keep_all = TRUE) %>%
          dplyr::arrange(dplyr::desc(.data[[expr_col]])) %>%
          dplyr::select(Gene.name, !!rlang::sym(expr_col))

        stats <- setNames(ranked_genes[[expr_col]], ranked_genes$Gene.name)

        es <- enrichmentHPA::run_enrichment(
          gene_list = gene_list,
          gene_stats = stats,
          method = "gsea",
          n_perm = n_perm,
          universe = universe
        )

        es$Tissue <- tissue
        return(es)
      } else {
        expressed_genes <- df_tissue %>%
          dplyr::filter(.data[[expr_col]] >= expression_threshold) %>%
          dplyr::pull(Gene.name) %>%
          unique()

        if (length(expressed_genes) < min_genes) return(NULL)

        res <- enrichmentHPA::run_enrichment(
          gene_list = gene_list,
          enriched_genes = expressed_genes,
          method = "hypergeometric",
          universe = universe
        )

        res$Tissue <- tissue
        return(res)
      }
    })

  if (!tibble::is_tibble(results)) {
    return(tibble::tibble(results))
  }


  return(results)
}
