#' Run enrichment analysis using hypergeometric test or GSEA-like score
#'
#' @param gene_list Character vector of user-supplied gene symbols.
#' @param universe Character vector of all genes in the background set.
#' @param method Enrichment method to use: "hypergeometric" or "gsea".
#' @param enriched_genes Required if method = "hypergeometric".
#' @param gene_stats Required if method = "gsea". Named numeric vector of gene scores (e.g., mean pTPM).
#' @param gsea_weight Exponent for GSEA running sum (default = 1).
#' @param n_perm Number of permutations for empirical p-value (GSEA only).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named list with enrichment results.
#' @export
run_enrichment <- function(gene_list,
                           universe,
                           method = c("hypergeometric", "gsea"),
                           enriched_genes = NULL,
                           gene_stats = NULL,
                           gsea_weight = 1,
                           n_perm = 0,
                           seed = NULL) {
  method <- match.arg(method)
  gene_list <- intersect(gene_list, universe)

  if (method == "hypergeometric") {
    if (is.null(enriched_genes)) stop("You must provide 'enriched_genes' for hypergeometric test.")
    enriched_genes <- intersect(enriched_genes, universe)

    N <- length(universe)
    K <- length(enriched_genes)
    n <- length(gene_list)
    k <- sum(gene_list %in% enriched_genes)

    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

    return(list(
      method = "hypergeometric",
      p_value = pval,
      overlap = k,
      input_set_size = n,
      enriched_set_size = K,
      universe_size = N
    ))
  }

  # GSEA
  if (is.null(gene_stats)) stop("You must provide 'gene_stats' for GSEA.")
  gene_stats <- gene_stats[names(gene_stats) %in% universe]
  ranked_genes <- sort(gene_stats, decreasing = TRUE)

  hits <- names(ranked_genes) %in% gene_list
  Nh <- sum(hits)
  N <- length(ranked_genes)

  if (Nh == 0) {
    warning("No overlap between gene_list and ranked gene_stats.")
    return(list(method = "gsea", enrichment_score = 0, p_value = NA, overlap = 0))
  }

  # Weighting
  weights <- abs(ranked_genes)^gsea_weight
  weights <- weights / sum(weights[hits])

  step_hit <- numeric(N)
  step_miss <- numeric(N)
  step_hit[hits] <- weights[hits]
  step_miss[!hits] <- 1 / (N - Nh)

  running_score <- cumsum(step_hit - step_miss)
  ES <- max(running_score, na.rm = TRUE)

  # Permutation-based p-value
  pval <- NA
  if (n_perm > 0) {
    if (!is.null(seed)) set.seed(seed)

    perm_ES <- replicate(n_perm, {
      random_genes <- sample(universe, Nh)
      random_hits <- names(ranked_genes) %in% random_genes

      step_hit <- numeric(N)
      step_miss <- numeric(N)
      step_hit[random_hits] <- weights[random_hits]
      step_miss[!random_hits] <- 1 / (N - Nh)

      max(cumsum(step_hit - step_miss), na.rm = TRUE)
    })

    pval <- (sum(perm_ES >= ES) + 1) / (n_perm + 1)
  }

  return(list(
    method = "gsea",
    enrichment_score = ES,
    p_value = pval,
    overlap = Nh,
    input_set_size = length(gene_list),
    universe_size = N,
    leading_edge = names(ranked_genes)[which.max(running_score)]
  ))
}
