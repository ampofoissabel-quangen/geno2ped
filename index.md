# geno2ped

## Overview

`geno2ped` is an R package for rapid pedigree construction from SNP
genotype data. It supports PLINK (bed/bim/fam) and VCF/VCF.GZ input
formats and infers familial relationships from genomic similarity.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ampofoissabel-quangen/geno2ped")
```

## Dependencies

Install required Bioconductor and CRAN packages:

``` r
# For PLINK files
install.packages("bigsnpr")

# For VCF files
install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
```

## Usage

### Read Genotypes from PLINK files

``` r
library(geno2ped)

data <- read_genotypes(
  bed_prefix = "path/to/mydata",        # points to mydata.bed/.bim/.fam
  sample_metadata = "path/to/meta.csv"  # optional: columns ID, Sex, BirthYear
)

# Access results
data$G    # genotype matrix (samples x SNPs), values 0/1/2
data$ids  # sample IDs
data$map  # SNP map (chr, pos, snp)
data$meta # sample metadata
```

### Read Genotypes from VCF

``` r
data <- read_genotypes(
  vcf = "path/to/mydata.vcf.gz",
  sample_metadata = "path/to/meta.csv"
)
```

### Sample Metadata Format

Your metadata CSV should look like this:

| ID        | Sex | BirthYear |
|-----------|-----|-----------|
| Sample001 | M   | 1985      |
| Sample002 | F   | 1990      |

## Input File Requirements

| Format | Required Files         |
|--------|------------------------|
| PLINK  | `.bed`, `.bim`, `.fam` |
| VCF    | `.vcf` or `.vcf.gz`    |

## Output

[`read_genotypes()`](https://ampofoissabel-quangen.github.io/geno2ped/reference/read_genotypes.md)
returns a list with: - `$G` — numeric matrix of genotypes (0/1/2), rows
= samples, cols = SNPs - `$ids` — character vector of sample IDs -
`$map` — data.frame with columns `chr`, `pos`, `snp` - `$meta` —
data.frame of sample metadata (or `NULL` if not provided)

## License

MIT © 2025 Issabelle Ampofo
