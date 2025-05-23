% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_expression.R
\name{summarize_expression}
\alias{summarize_expression}
\title{Summarize expression levels and categorize by strength}
\usage{
summarize_expression(
  data,
  group_col,
  gene_col = "Gene.name",
  expr_col = "pTPM",
  thresholds = c(1, 10)
)
}
\arguments{
\item{data}{A data.frame with HPA expression values (pTPM).}

\item{group_col}{Column to group by (e.g. "Immune.cell", "Region").}

\item{gene_col}{Column for gene names (default: "Gene.name").}

\item{expr_col}{Column for expression values (default: "pTPM").}

\item{thresholds}{Numeric vector of length 2: c(medium_cutoff, high_cutoff)}
}
\value{
A data.frame with gene stats and expression category.
}
\description{
Computes mean pTPM and expression category for each gene across groups.
}
