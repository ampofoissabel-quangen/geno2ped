# tests/testthat/test-read_genotypes.R

library(testthat)
library(geno2ped)

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

test_that("read_genotypes() returns correct output structure", {
  tmp  <- tempdir()
  mock <- list(
    G   = matrix(sample(0:2, 40, replace = TRUE), nrow = 8,
                 dimnames = list(paste0("IND", 1:8), paste0("SNP", 1:5))),
    ids = paste0("IND", 1:8),
    map = data.frame(chr = 1, pos = 1:5, snp = paste0("SNP", 1:5))
  )
  saveRDS(mock, file.path(tmp, "structure_test.rds"))

  result <- read_genotypes(bed_prefix = file.path(tmp, "structure_test"))

  expect_type(result, "list")
  expect_named(result, c("G", "ids", "map", "meta"))
  expect_true(is.matrix(result$G))
  expect_type(result$ids, "character")
  expect_true(is.data.frame(result$map))
})

test_that("read_genotypes() genotype values are only 0, 1, 2, or NA", {
  tmp  <- tempdir()
  mock <- list(
    G   = matrix(sample(0:2, 40, replace = TRUE), nrow = 8,
                 dimnames = list(paste0("IND", 1:8), paste0("SNP", 1:5))),
    ids = paste0("IND", 1:8),
    map = data.frame(chr = 1, pos = 1:5, snp = paste0("SNP", 1:5))
  )
  saveRDS(mock, file.path(tmp, "values_test.rds"))

  result <- read_genotypes(bed_prefix = file.path(tmp, "values_test"))
  expect_true(all(result$G %in% c(0, 1, 2, NA)))
})

test_that("read_genotypes() returns NULL meta when not provided", {
  tmp  <- tempdir()
  mock <- list(
    G   = matrix(0L, 4, 4, dimnames = list(paste0("S", 1:4), paste0("SNP", 1:4))),
    ids = paste0("S", 1:4),
    map = data.frame(chr = 1, pos = 1:4, snp = paste0("SNP", 1:4))
  )
  saveRDS(mock, file.path(tmp, "null_meta.rds"))

  result <- read_genotypes(bed_prefix = file.path(tmp, "null_meta"))
  expect_null(result$meta)
})

test_that("read_genotypes() loads metadata and returns correct columns", {
  tmp  <- tempdir()
  ids  <- paste0("IND", 1:6)
  mock <- list(
    G   = matrix(sample(0:2, 30, replace = TRUE), nrow = 6,
                 dimnames = list(ids, paste0("SNP", 1:5))),
    ids = ids,
    map = data.frame(chr = 1, pos = 1:5, snp = paste0("SNP", 1:5))
  )
  saveRDS(mock, file.path(tmp, "with_meta.rds"))

  meta_csv <- file.path(tmp, "meta.csv")
  write.csv(data.frame(ID = ids, Sex = rep(c("M","F"), 3),
                       BirthYear = 1980:1985), meta_csv, row.names = FALSE)

  result <- read_genotypes(
    bed_prefix      = file.path(tmp, "with_meta"),
    sample_metadata = meta_csv
  )

  expect_false(is.null(result$meta))
  expect_true(all(c("ID", "Sex", "BirthYear") %in% colnames(result$meta)))
  expect_equal(nrow(result$meta), 6)
})
