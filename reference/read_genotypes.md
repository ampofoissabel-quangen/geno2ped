# Read genotypes and sample metadata

Read genotypes and sample metadata

## Usage

``` r
read_genotypes(
  bed_prefix = NULL,
  vcf = NULL,
  csv = NULL,
  sample_metadata = NULL
)
```

## Arguments

- bed_prefix:

  Path prefix to PLINK bed/bim/fam (optional)

- vcf:

  Path to VCF/VCF.GZ (optional)

- csv:

  Path to CSV file with sample IDs in first column and SNPs in remaining
  columns (optional)

- sample_metadata:

  CSV with columns: ID, Sex (M/F), BirthYear (optional)

## Value

A list with \$G (matrix 0/1/2), \$ids (character), \$map (data.frame),
\$meta (data.frame)
