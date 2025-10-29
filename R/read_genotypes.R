#' Read genotypes and sample metadata
#'
#' @param bed_prefix Path prefix to PLINK bed/bim/fam (optional)
#' @param vcf Path to VCF/VCF.GZ (optional)
#' @param sample_metadata CSV with columns: ID, Sex (M/F), BirthYear (optional)
#' @return A list with $G (matrix 0/1/2), $ids (character), $map (data.frame), $meta (data.frame)
#' @export
read_genotypes <- function(bed_prefix = NULL, vcf = NULL, sample_metadata = NULL) {
  if (is.null(bed_prefix) && is.null(vcf)) {
    stop("Provide either 'bed_prefix' or 'vcf'.")
  }
  meta <- if (!is.null(sample_metadata)) utils::read.csv(sample_metadata, stringsAsFactors = FALSE) else NULL

  # Minimal mock reader: expect an RDS at bed_prefix.rds or vcf.rds for prototype use
  # Replace with real PLINK/VCF readers (e.g., SNPRelate, bigsnpr, VariantAnnotation) as needed.
  if (!is.null(bed_prefix)) {
    rds <- paste0(bed_prefix, ".rds")
    if (!file.exists(rds)) stop("For the prototype, place a 0/1/2 genotype matrix at ", rds)
    obj <- readRDS(rds)
  } else {
    rds <- paste0(vcf, ".rds")
    if (!file.exists(rds)) stop("For the prototype, place a 0/1/2 genotype matrix at ", rds)
    obj <- readRDS(rds)
  }

  # Expect obj: list(G = matrix, ids = character, map = data.frame(chr, pos, snp))
  if (is.null(obj$G) || is.null(obj$ids) || is.null(obj$map)) stop("Malformed input RDS: needs G, ids, map.")
  list(G = obj$G, ids = obj$ids, map = obj$map, meta = meta)
}