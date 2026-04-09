# tests/testthat/test-basic-workflow.R

library(testthat)
library(geno2ped)

test_that("basic geno2ped workflow runs with toy data", {
  geno_file <- system.file("extdata", "toy_genotypes.csv", package = "geno2ped")
  meta_file <- system.file("extdata", "toy_metadata.csv", package = "geno2ped")

  skip_if(
    !nzchar(geno_file) || !file.exists(geno_file),
    "toy_genotypes.csv not found in inst/extdata — skipping"
  )

  # Use csv argument (not positional) since read_genotypes expects bed_prefix or vcf
  geno <- read_genotypes(csv = geno_file, sample_metadata = meta_file)

  expect_type(geno, "list")
  expect_named(geno, c("G", "ids", "map", "meta"))
  expect_true(nrow(geno$G) > 0)

  res <- build_pedigree(
    geno,
    preset  = "high_precision",
    verbose = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("pedigree", "assignments", "trios", "summary", "settings") %in% names(res)))
  expect_true(nrow(res$pedigree) == nrow(geno$G))  # one row per individual
  expect_true("assignment_status" %in% names(res$pedigree))
  expect_true(is.list(res$summary))
  expect_true("assignment_rate" %in% names(res$summary))
  expect_gte(res$summary$assignment_rate, 0)
  expect_lte(res$summary$assignment_rate, 1)
})

