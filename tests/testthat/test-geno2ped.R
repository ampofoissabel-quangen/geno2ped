# tests/testthat/test-geno2ped.R

library(testthat)
library(geno2ped)

# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

make_geno <- function(n_samples = 8, n_snps = 50, seed = 42) {
  set.seed(seed)
  ids <- paste0("IND", seq_len(n_samples))
  G   <- matrix(
    sample(0:2, n_samples * n_snps, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
    nrow = n_samples, ncol = n_snps,
    dimnames = list(ids, paste0("SNP", seq_len(n_snps)))
  )
  map <- data.frame(
    chr = 1L,
    pos = seq(1000L, by = 500L, length.out = n_snps),
    snp = paste0("SNP", seq_len(n_snps)),
    stringsAsFactors = FALSE
  )
  list(G = G, ids = ids, map = map, meta = NULL)
}

make_meta <- function(ids) {
  n <- length(ids)
  data.frame(
    ID        = ids,
    Sex       = rep(c("M", "F"), length.out = n),
    BirthYear = c(rep(1990L, n %/% 2), rep(2010L, n - n %/% 2)),
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# read_genotypes()
# ══════════════════════════════════════════════════════════════════════════════

test_that("read_genotypes() errors with no input", {
  expect_error(read_genotypes(), "Provide either")
})

test_that("read_genotypes() errors when PLINK files are missing", {
  expect_error(
    read_genotypes(bed_prefix = "nonexistent/data"),
    "Missing PLINK file"
  )
})

test_that("read_genotypes() errors when VCF file does not exist", {
  expect_error(
    read_genotypes(vcf = "nonexistent.vcf.gz"),
    "VCF file not found"
  )
})

test_that("read_genotypes() loads metadata CSV correctly", {
  tmp      <- tempdir()
  mock     <- make_geno()
  meta     <- make_meta(mock$ids)
  meta_csv <- file.path(tmp, "meta.csv")
  write.csv(meta, meta_csv, row.names = FALSE)

  saveRDS(mock, file.path(tmp, "meta_test.rds"))

  result <- read_genotypes(
    bed_prefix      = file.path(tmp, "meta_test"),
    sample_metadata = meta_csv
  )

  expect_false(is.null(result$meta))
  expect_equal(nrow(result$meta), length(mock$ids))
  expect_true(all(c("ID", "Sex", "BirthYear") %in% colnames(result$meta)))
})

test_that("read_genotypes() returns NULL meta when not provided", {
  tmp  <- tempdir()
  mock <- make_geno()
  saveRDS(mock, file.path(tmp, "no_meta.rds"))

  result <- read_genotypes(bed_prefix = file.path(tmp, "no_meta"))
  expect_null(result$meta)
})

# ══════════════════════════════════════════════════════════════════════════════
# build_pedigree()
# ══════════════════════════════════════════════════════════════════════════════

test_that("build_pedigree() returns correct list structure", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("pedigree", "assignments", "trios", "trios_all", "summary", "settings"))
})

test_that("build_pedigree() pedigree has correct columns", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  expect_true(all(c("ID", "Sire", "Dam", "assignment_status") %in% colnames(res$pedigree)))
})

test_that("build_pedigree() pedigree contains all input individuals", {
  geno <- make_geno(n_samples = 8)
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  expect_equal(nrow(res$pedigree), 8)
  expect_setequal(res$pedigree$ID, geno$ids)
})

test_that("build_pedigree() assignment_status is only complete or unassigned", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  expect_true(all(res$pedigree$assignment_status %in% c("complete", "unassigned")))
})

test_that("build_pedigree() summary has expected fields", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  expect_true(all(c("n_individuals", "n_complete_assignments",
                    "assignment_rate", "mean_ME_accepted") %in% names(res$summary)))
})

test_that("build_pedigree() assignment_rate is between 0 and 1", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  expect_gte(res$summary$assignment_rate, 0)
  expect_lte(res$summary$assignment_rate, 1)
})

test_that("build_pedigree() high_precision preset sets correct thresholds", {
  geno      <- make_geno()
  geno$meta <- make_meta(geno$ids)  # needed because high_precision uses use_age = TRUE
  res       <- build_pedigree(geno, preset = "high_precision", verbose = FALSE)

  expect_equal(res$settings$s_threshold, 0.85)
  expect_equal(res$settings$me_max, 0.002)
  expect_equal(res$settings$top_k, 2)
  expect_true(res$settings$use_age)
})

test_that("build_pedigree() high_coverage preset sets correct thresholds", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "high_coverage", verbose = FALSE)

  expect_equal(res$settings$s_threshold, 0.70)
  expect_equal(res$settings$oh_max, 50)
})

test_that("build_pedigree() respects sex filtering with metadata", {
  geno       <- make_geno(n_samples = 8)
  geno$meta  <- make_meta(geno$ids)
  res        <- build_pedigree(geno, preset = "none", use_sex = TRUE, verbose = FALSE)

  if (nrow(res$assignments)) {
    sire_ids <- res$assignments$Parent[res$assignments$Role == "Sire"]
    dam_ids  <- res$assignments$Parent[res$assignments$Role == "Dam"]
    males    <- geno$meta$ID[tolower(geno$meta$Sex) == "m"]
    females  <- geno$meta$ID[tolower(geno$meta$Sex) == "f"]

    expect_true(all(sire_ids %in% males))
    expect_true(all(dam_ids  %in% females))
  }
})

test_that("build_pedigree() errors on invalid geno object", {
  expect_error(build_pedigree(list(G = NULL)), regexp = NULL)
})

test_that("build_pedigree() trios ME values are between 0 and 1", {
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  if (nrow(res$trios_all)) {
    me_vals <- res$trios_all$ME[is.finite(res$trios_all$ME)]
    expect_true(all(me_vals >= 0 & me_vals <= 1))
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# write_outputs()
# ══════════════════════════════════════════════════════════════════════════════

test_that("write_outputs() creates expected files", {
  geno   <- make_geno()
  res    <- build_pedigree(geno, preset = "none", verbose = FALSE)
  outdir <- file.path(tempdir(), "test_outputs")

  write_outputs(res, outdir = outdir)

  expect_true(file.exists(file.path(outdir, "pedigree.csv")))
  expect_true(file.exists(file.path(outdir, "assignments.csv")))
  expect_true(file.exists(file.path(outdir, "trios.csv")))
  expect_true(file.exists(file.path(outdir, "summary.txt")))
  expect_true(file.exists(file.path(outdir, "settings.txt")))
})

test_that("write_outputs() respects file prefix", {
  geno   <- make_geno()
  res    <- build_pedigree(geno, preset = "none", verbose = FALSE)
  outdir <- file.path(tempdir(), "test_prefix_outputs")

  write_outputs(res, outdir = outdir, prefix = "run1")

  expect_true(file.exists(file.path(outdir, "run1_pedigree.csv")))
  expect_true(file.exists(file.path(outdir, "run1_assignments.csv")))
})

test_that("write_outputs() pedigree CSV has correct columns", {
  geno   <- make_geno()
  res    <- build_pedigree(geno, preset = "none", verbose = FALSE)
  outdir <- file.path(tempdir(), "test_csv_check")

  write_outputs(res, outdir = outdir)

  ped <- read.csv(file.path(outdir, "pedigree.csv"))
  expect_true(all(c("ID", "Sire", "Dam", "assignment_status") %in% colnames(ped)))
})

test_that("write_outputs() errors on invalid result object", {
  expect_error(
    write_outputs(list(pedigree = NULL)),
    "does not look like output from build_pedigree"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# plot functions (smoke tests)
# ══════════════════════════════════════════════════════════════════════════════

test_that("plot_kinship() returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", s_threshold = 0.1, verbose = FALSE)

  if (nrow(res$assignments)) {
    p <- plot_kinship(res)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_me() returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", s_threshold = 0.1, verbose = FALSE)

  if (nrow(res$trios_all)) {
    p <- plot_me(res)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_assignment_status() returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  geno <- make_geno()
  res  <- build_pedigree(geno, preset = "none", verbose = FALSE)

  p <- plot_assignment_status(res)
  expect_s3_class(p, "ggplot")
})
