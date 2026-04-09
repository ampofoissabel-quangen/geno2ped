test_that("basic geno2ped workflow runs", {
  geno_file <- system.file("extdata", "toy_genotypes.csv", package = "geno2ped")
  meta_file <- system.file("extdata", "toy_metadata.csv", package = "geno2ped")

  geno <- read_genotypes(geno_file, sample_metadata = meta_file)

  res <- build_pedigree(
    geno,
    preset = "high_precision",
    verbose = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("pedigree", "assignments", "trios", "summary", "settings") %in% names(res)))

  expect_true(nrow(res$pedigree) == 6)
  expect_true("assignment_status" %in% names(res$pedigree))

  expect_true(is.list(res$summary))
  expect_true("assignment_rate" %in% names(res$summary))
})
