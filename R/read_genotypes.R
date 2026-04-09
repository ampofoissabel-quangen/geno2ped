#' Read genotypes and optional sample metadata
#'
#' @param geno_file Path to genotype matrix CSV/TSV. First column must be `ID`.
#'   SNPs must be coded as 0, 1, 2, or NA.
#' @param sample_metadata Optional CSV/TSV with columns `ID`, and optionally
#'   `Sex` and `BirthYear`.
#' @param sep Field separator. If NULL, guessed from file extension.
#'
#' @return A list with:
#' \describe{
#'   \item{G}{integer matrix of genotypes coded 0/1/2 with rownames = IDs}
#'   \item{ids}{character vector of individual IDs}
#'   \item{map}{data.frame with SNP column names}
#'   \item{meta}{data.frame of sample metadata or NULL}
#' }
#' @export
read_genotypes <- function(geno_file, sample_metadata = NULL, sep = NULL) {
  if (missing(geno_file) || is.null(geno_file)) {
    stop("Please provide 'geno_file'.", call. = FALSE)
  }
  if (!file.exists(geno_file)) {
    stop("Genotype file not found: ", geno_file, call. = FALSE)
  }

  if (is.null(sep)) {
    sep <- if (grepl("\\.tsv$|\\.txt$", geno_file, ignore.case = TRUE)) "\t" else ","
  }

  dat <- utils::read.table(
    geno_file,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"ID" %in% names(dat)) {
    stop("Genotype file must contain a column named 'ID'.", call. = FALSE)
  }

  if (anyDuplicated(dat$ID)) {
    dup <- unique(dat$ID[duplicated(dat$ID)])
    stop("Duplicate IDs found in genotype file: ", paste(dup, collapse = ", "), call. = FALSE)
  }

  snp_cols <- setdiff(names(dat), "ID")
  if (!length(snp_cols)) {
    stop("No SNP columns found in genotype file.", call. = FALSE)
  }

  G <- as.matrix(dat[, snp_cols, drop = FALSE])
  mode(G) <- "numeric"

  bad <- !(is.na(G) | G %in% c(0, 1, 2))
  if (any(bad)) {
    stop("Genotypes must be coded as 0, 1, 2, or NA.", call. = FALSE)
  }

  G <- matrix(as.integer(G), nrow = nrow(G), ncol = ncol(G),
              dimnames = list(dat$ID, snp_cols))

  meta <- NULL
  if (!is.null(sample_metadata)) {
    if (!file.exists(sample_metadata)) {
      stop("Metadata file not found: ", sample_metadata, call. = FALSE)
    }

    sep_meta <- if (grepl("\\.tsv$|\\.txt$", sample_metadata, ignore.case = TRUE)) "\t" else ","
    meta <- utils::read.table(
      sample_metadata,
      header = TRUE,
      sep = sep_meta,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (!"ID" %in% names(meta)) {
      stop("Metadata file must contain a column named 'ID'.", call. = FALSE)
    }

    if (anyDuplicated(meta$ID)) {
      dup <- unique(meta$ID[duplicated(meta$ID)])
      stop("Duplicate IDs found in metadata file: ", paste(dup, collapse = ", "), call. = FALSE)
    }

    meta <- meta[match(dat$ID, meta$ID), , drop = FALSE]
  }

  map <- data.frame(
    SNP = snp_cols,
    stringsAsFactors = FALSE
  )

  list(
    G = G,
    ids = dat$ID,
    map = map,
    meta = meta
  )
}
