#' Enrichment analysis across cell types using single-cell reference
#'
#' This function tests whether a gene list is significantly enriched in any cell type
#' using single-cell expression data (e.g., Tabula Sapiens).
#'
#' @param gene_list Character vector of genes of interest (e.g., DEGs).
#' @param reference A long-format data frame with columns: Gene, Cell, Expression.
#' @param method Enrichment method: "hypergeometric" (default) or "gsea".
#' @param expr_col Name of the column in `reference` with expression values (default = "Expression").
#' @param expression_threshold Minimum expression to define a gene as "expressed" in a cell type (default = 0.3).
#' @param universe Optional character vector of background genes.
#' @param gsea_weight Weighting parameter for GSEA statistic (default = 1).
#' @param n_perm Number of permutations for GSEA p-value estimation (default = 0 = skip).
#' @param seed Random seed for reproducibility (optional).
#'
#' @return A data.frame of enrichment results per cell type.
#' @export
enrich_by_celltype <- function(gene_list,
                               reference,
                               method = c("hypergeometric", "gsea"),
                               expr_col = "Expression",
                               gene_stats = NULL,
                               expression_threshold = 0.3,
                               universe = NULL,
                               universe_global = T,
                               gsea_weight = 1,
                               n_perm = 0,
                               seed = NULL) {
  method <- match.arg(method)

  # Validate required columns
  required_cols <- c("Gene", "CellType", expr_col)
  missing_cols <- setdiff(required_cols, names(reference))
  if (length(missing_cols) > 0) {
    stop("Missing columns in reference: ", paste(missing_cols, collapse = ", "))
  }

  reference <- reference %>%
    dplyr::filter(!is.na(.data[[expr_col]]))


  if(is.null(universe)) {
    global_universe <- reference %>%
      dplyr::filter(.data[[expr_col]] >= expression_threshold) %>%
      dplyr::pull(Gene) %>%
      unique()
  }

  # Group by CellType and perform enrichment for each
  results <- reference %>%
    dplyr::group_by(CellType) %>%
    dplyr::group_map(~{
      cell_df <- .x
      cell_type <- unique(cell_df$CellType)
      print(cell_type)
      # Identify expressed genes in this cell type
      expressed_genes <- cell_df %>%
        dplyr::filter(.data[[expr_col]] >= expression_threshold) %>%
        dplyr::pull(Gene) %>%
        unique()

      if (length(expressed_genes) == 0) return(NULL)

      universe <- if(universe_global) global_universe else unique(cell_df$Gene)

      # GSEA gene stats handling
      if (method == "gsea" && is.null(gene_stats)) {
        local_gene_stats <- tapply(cell_df[[expr_col]], cell_df$Gene, mean, na.rm = TRUE)
      } else {
        local_gene_stats <- gene_stats[names(gene_stats) %in% universe]
      }



      res <- enrichmentHPA::run_enrichment(
        gene_list = gene_list,
        enriched_genes = expressed_genes,
        gene_stats = local_gene_stats,
        method = method,
        universe = universe,
        gsea_weight = gsea_weight,
        n_perm = n_perm,
        seed = seed
      )
      if (is.null(res) || length(res) == 0) return(NULL)
      res_df <- as.data.frame(res, stringsAsFactors = FALSE)
      res_df$CellType <- cell_type
      res_df
    }, .keep=T) %>%
    dplyr::bind_rows()



  return(results)
}
