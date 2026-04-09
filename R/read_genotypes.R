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

  # Load sample metadata if provided
  meta <- if (!is.null(sample_metadata)) {
    utils::read.csv(sample_metadata, stringsAsFactors = FALSE)
  } else NULL

  # ── PLINK bed/bim/fam ────────────────────────────────────────────────────────
  if (!is.null(bed_prefix)) {

    # Check all three PLINK files exist
    for (ext in c(".bed", ".bim", ".fam")) {
      f <- paste0(bed_prefix, ext)
      if (!file.exists(f)) stop("Missing PLINK file: ", f)
    }

    if (!requireNamespace("bigsnpr", quietly = TRUE)) {
      stop("Package 'bigsnpr' is required for PLINK files. Install with:\n",
           "  install.packages('bigsnpr')")
    }

    # Read PLINK files using bigsnpr
    rds_out <- paste0(bed_prefix, "_bigsnpr.rds")
    snp_obj <- bigsnpr::snp_readBed(paste0(bed_prefix, ".bed"), backingfile = sub("\\.rds$", "", rds_out))
    snp_obj <- bigsnpr::snp_attach(rds_out)

    G_raw   <- bigsnpr::snp_cor(snp_obj$genotypes)  # FBM object
    G       <- as.matrix(snp_obj$genotypes[])        # Convert to 0/1/2 matrix

    # Build map from .bim file
    bim <- utils::read.table(paste0(bed_prefix, ".bim"),
                             col.names = c("chr", "snp", "cm", "pos", "A1", "A2"),
                             stringsAsFactors = FALSE)
    map <- bim[, c("chr", "pos", "snp")]

    # Get sample IDs from .fam file
    fam <- utils::read.table(paste0(bed_prefix, ".fam"),
                             col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"),
                             stringsAsFactors = FALSE)
    ids <- fam$IID

    rownames(G) <- ids
    colnames(G) <- bim$snp

    # ── VCF / VCF.GZ ─────────────────────────────────────────────────────────────
  } else {

    if (!file.exists(vcf)) stop("VCF file not found: ", vcf)

    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
      stop("Package 'VariantAnnotation' is required for VCF files. Install with:\n",
           "  BiocManager::install('VariantAnnotation')")
    }

    vcf_obj  <- VariantAnnotation::readVcf(vcf)
    geno_raw <- VariantAnnotation::geno(vcf_obj)$GT  # Genotype calls (e.g. "0/1")

    # Convert GT strings to 0/1/2 dosage
    gt_to_dosage <- function(gt) {
      gt <- gsub("\\|", "/", gt)  # Handle phased genotypes
      sapply(gt, function(g) {
        if (is.na(g) || g %in% c("./.", ".")) return(NA_integer_)
        alleles <- strsplit(g, "/")[[1]]
        sum(alleles != "0", na.rm = TRUE)
      })
    }

    G <- t(apply(geno_raw, 1, gt_to_dosage))  # SNPs x Samples → transpose to Samples x SNPs
    G <- t(G)                                  # Final: Samples x SNPs
    storage.mode(G) <- "integer"

    # Build map from VCF rowRanges
    rd  <- SummarizedExperiment::rowRanges(vcf_obj)
    map <- data.frame(
      chr = as.character(GenomicRanges::seqnames(rd)),
      pos = GenomicRanges::start(rd),
      snp = rownames(geno_raw),
      stringsAsFactors = FALSE
    )

    ids <- colnames(geno_raw)
    rownames(G) <- ids
    colnames(G) <- map$snp
  }

  # ── Return ────────────────────────────────────────────────────────────────────
  list(G = G, ids = ids, map = map, meta = meta)
}
