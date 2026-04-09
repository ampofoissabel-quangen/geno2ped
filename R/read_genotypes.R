#' Read genotypes and sample metadata
#'
#' @param bed_prefix Path prefix to PLINK bed/bim/fam (optional)
#' @param vcf Path to VCF/VCF.GZ (optional)
#' @param csv Path to CSV file with sample IDs in first column and SNPs in remaining columns (optional)
#' @param sample_metadata CSV with columns: ID, Sex (M/F), BirthYear (optional)
#' @return A list with $G (matrix 0/1/2), $ids (character), $map (data.frame), $meta (data.frame)
#' @importFrom utils read.csv read.table
#' @importFrom GenomicRanges seqnames start
#' @importFrom SummarizedExperiment rowRanges
#' @export

read_genotypes <- function(bed_prefix = NULL, vcf = NULL,
                           csv = NULL, sample_metadata = NULL) {

  if (is.null(bed_prefix) && is.null(vcf) && is.null(csv)) {
    stop("Provide either 'bed_prefix', 'vcf', or 'csv'.")
  }

  meta <- if (!is.null(sample_metadata)) {
    utils::read.csv(sample_metadata, stringsAsFactors = FALSE)
  } else NULL

  # ── CSV input ────────────────────────────────────────────────────────────────
  if (!is.null(csv)) {
    if (!file.exists(csv)) stop("CSV file not found: ", csv)
    raw  <- utils::read.csv(csv, stringsAsFactors = FALSE)
    ids  <- raw[, 1]
    G    <- as.matrix(raw[, -1])
    storage.mode(G) <- "integer"
    rownames(G) <- ids
    map <- data.frame(
      chr = NA_character_,
      pos = NA_integer_,
      snp = colnames(G),
      stringsAsFactors = FALSE
    )
    return(list(G = G, ids = ids, map = map, meta = meta))
  }

  # ── PLINK bed/bim/fam ────────────────────────────────────────────────────────
  if (!is.null(bed_prefix)) {

    # Allow .rds shortcut for testing/prototype
    rds <- paste0(bed_prefix, ".rds")
    if (file.exists(rds)) {
      obj <- readRDS(rds)
      if (is.null(obj$G) || is.null(obj$ids) || is.null(obj$map))
        stop("Malformed input RDS: needs G, ids, map.")
      return(list(G = obj$G, ids = obj$ids, map = obj$map, meta = meta))
    }

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
    snp_obj <- bigsnpr::snp_readBed(
      paste0(bed_prefix, ".bed"),
      backingfile = sub("\\.rds$", "", rds_out)
    )
    snp_obj <- bigsnpr::snp_attach(rds_out)
    G       <- as.matrix(snp_obj$genotypes[])

    # Build map from .bim file
    bim <- utils::read.table(
      paste0(bed_prefix, ".bim"),
      col.names = c("chr", "snp", "cm", "pos", "A1", "A2"),
      stringsAsFactors = FALSE
    )
    map <- bim[, c("chr", "pos", "snp")]

    # Get sample IDs from .fam file
    fam <- utils::read.table(
      paste0(bed_prefix, ".fam"),
      col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"),
      stringsAsFactors = FALSE
    )
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
    geno_raw <- VariantAnnotation::geno(vcf_obj)$GT

    # Convert GT strings to 0/1/2 dosage
    gt_to_dosage <- function(gt) {
      gt <- gsub("\\|", "/", gt)
      sapply(gt, function(g) {
        if (is.na(g) || g %in% c("./.", ".")) return(NA_integer_)
        alleles <- strsplit(g, "/")[[1]]
        sum(alleles != "0", na.rm = TRUE)
      })
    }

    G <- t(apply(geno_raw, 1, gt_to_dosage))
    G <- t(G)
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
