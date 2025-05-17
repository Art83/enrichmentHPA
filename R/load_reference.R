#' Download and load HPA RNA-seq data
#'
#' Downloads and reads a compressed RNA-seq dataset from the Human Protein Atlas.
#'
#' @param file Name of the file to be saved locally (e.g., "rna_ic.zip").
#' @param type Type of dataset to download. Options are:
#'   \itemize{
#'     \item \code{"hpa_brain"} – RNA across 13 brain regions
#'     \item \code{"rna_sn_clusters_brain"} – RNA across 34 single nuclei cluster types in 11 brain regions
#'     \item \code{"rna_ic"} – RNA expression in immune cell types
#'     \item \code{"rna_ic_monaco"} - RNA expression in subsets of immune cell types
#'     \item \code{"rna_sn_clusters_organs"} - RNA expression across clusters in different organs
#'   }
#' @param dest_dir Local directory to store the downloaded file (default: \code{tempdir()}).
#' @param overwrite Logical; if \code{TRUE}, the file will be re-downloaded even if it exists.
#'
#' @return A \code{data.frame} containing the unzipped RNA-seq dataset.
#' @export
load_reference <- function(file,
                           type = "hpa_brain",
                           dest_dir = tempdir(),
                           overwrite = FALSE) {
  urls_list <- list(
    hpa_brain = "https://www.proteinatlas.org/download/tsv/rna_brain_region_hpa.tsv.zip",
    rna_sn_clusters_brain = "https://www.proteinatlas.org/download/tsv/rna_single_nuclei_cluster_type.tsv.zip",
    rna_ic = "https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip",
    rna_ic_monaco = "https://www.proteinatlas.org/download/tsv/rna_immune_cell_monaco.tsv.zip",
    rna_sn_clusters_organs = "https://www.proteinatlas.org/download/tsv/rna_single_cell_cluster.tsv.zip"
  )

  # Check for valid type
  if (!type %in% names(urls_list)) {
    stop("Invalid 'type' argument. Choose from: ", paste(names(urls_list), collapse = ", "))
  }

  base_url <- urls_list[[type]]
  dest_file <- file.path(dest_dir, file)

  # Download if needed
  if (!file.exists(dest_file) || overwrite) {
    message("Downloading data from: ", base_url)
    download.file(url = base_url, destfile = dest_file, mode = "wb")
  } else {
    message("Using cached file: ", dest_file)
  }

  # Unzip and read first TSV file
  unzipped_file <- unzip(dest_file, exdir = dest_dir)
  data <- read.delim(unzipped_file[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  return(data)
}
