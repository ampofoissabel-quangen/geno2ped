# Introduction to geno2ped

## Overview

`geno2ped` is an R package for rapid pedigree construction from SNP
genotype data. It infers parent-offspring relationships using genomic
similarity and Mendelian error validation, and supports PLINK, VCF, and
CSV input formats.

## Installation

``` r
devtools::install_github("ampofoissabel-quangen/geno2ped")
```

## Loading the Package

``` r
library(geno2ped)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
```

## Reading Genotypes

### From a CSV file

The simplest input format is a CSV where the first column contains
sample IDs and the remaining columns contain SNP genotypes (0/1/2):

``` r
geno <- read_genotypes(
  csv             = "path/to/genotypes.csv",
  sample_metadata = "path/to/metadata.csv"
)
```

### From PLINK files

``` r
geno <- read_genotypes(
  bed_prefix      = "path/to/mydata",
  sample_metadata = "path/to/metadata.csv"
)
```

### From VCF files

``` r
geno <- read_genotypes(
  vcf             = "path/to/mydata.vcf.gz",
  sample_metadata = "path/to/metadata.csv"
)
```

### The genotype object

[`read_genotypes()`](https://ampofoissabel-quangen.github.io/geno2ped/reference/read_genotypes.md)
returns a list with four elements:

| Element | Description                             |
|---------|-----------------------------------------|
| `$G`    | Genotype matrix (samples x SNPs, 0/1/2) |
| `$ids`  | Sample ID vector                        |
| `$map`  | SNP map (chr, pos, snp)                 |
| `$meta` | Sample metadata (ID, Sex, BirthYear)    |

### Sample metadata format

Your metadata CSV should look like this:

| ID        | Sex | BirthYear |
|-----------|-----|-----------|
| Sample001 | M   | 1985      |
| Sample002 | F   | 1990      |
| Sample003 | M   | 2005      |

## Building a Pedigree

### Using a preset

The easiest way is to use a built-in preset:

``` r
res <- build_pedigree(
  geno,
  preset  = "high_precision",
  verbose = TRUE
)
```

Three presets are available:

| Preset           | Use case                 |
|------------------|--------------------------|
| `none`           | Manual threshold control |
| `high_precision` | Minimise false positives |
| `high_coverage`  | Maximise assignments     |

### Custom thresholds

``` r
res <- build_pedigree(
  geno,
  s_threshold    = 0.75,
  me_max         = 0.003,
  oh_max         = 20,
  use_sex        = TRUE,
  use_age        = TRUE,
  min_parent_gap = 2,
  top_k          = 3
)
```

### The result object

[`build_pedigree()`](https://ampofoissabel-quangen.github.io/geno2ped/reference/build_pedigree.md)
returns a list with:

| Element        | Description                            |
|----------------|----------------------------------------|
| `$pedigree`    | Final pedigree (ID, Sire, Dam, status) |
| `$assignments` | All candidate parent scores            |
| `$trios`       | Best trio per offspring                |
| `$trios_all`   | All evaluated trios                    |
| `$summary`     | Assignment statistics                  |
| `$settings`    | Parameters used                        |

### Viewing results

``` r
# View the pedigree
res$pedigree

# View summary statistics
res$summary
```

## Visualising Results

``` r
# Similarity score distribution
plot_kinship(res)

# Mendelian error rate distribution
plot_me(res)

# Assignment outcome counts
plot_assignment_status(res)
```

## Saving Results

``` r
write_outputs(
  res,
  outdir = "geno2ped_results",
  prefix = "my_run"
)
```

This creates the following files in `geno2ped_results/`:

- `my_run_pedigree.csv` — final pedigree
- `my_run_assignments.csv` — candidate scores
- `my_run_trios.csv` — best trios
- `my_run_summary.txt` — run statistics
- `my_run_settings.txt` — parameters used

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] geno2ped_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-4                gtable_0.3.6               
#>  [3] jsonlite_2.0.0              compiler_4.5.3             
#>  [5] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [7] GenomicRanges_1.62.1        jquerylib_0.1.4            
#>  [9] scales_1.4.0                systemfonts_1.3.2          
#> [11] IRanges_2.44.0              Seqinfo_1.0.0              
#> [13] textshaping_1.0.5           yaml_2.3.12                
#> [15] fastmap_1.2.0               lattice_0.22-9             
#> [17] ggplot2_4.0.2               R6_2.6.1                   
#> [19] XVector_0.50.0              S4Arrays_1.10.1            
#> [21] generics_0.1.4              knitr_1.51                 
#> [23] BiocGenerics_0.56.0         DelayedArray_0.36.1        
#> [25] desc_1.4.3                  MatrixGenerics_1.22.0      
#> [27] RColorBrewer_1.1-3          bslib_0.10.0               
#> [29] rlang_1.2.0                 cachem_1.1.0               
#> [31] xfun_0.57                   S7_0.2.1                   
#> [33] fs_2.0.1                    sass_0.4.10                
#> [35] SparseArray_1.10.10         cli_3.6.5                  
#> [37] pkgdown_2.2.0               digest_0.6.39              
#> [39] grid_4.5.3                  lifecycle_1.0.5            
#> [41] vctrs_0.7.2                 S4Vectors_0.48.1           
#> [43] glue_1.8.0                  evaluate_1.0.5             
#> [45] farver_2.1.2                ragg_1.5.2                 
#> [47] abind_1.4-8                 stats4_4.5.3               
#> [49] rmarkdown_2.31              matrixStats_1.5.0          
#> [51] tools_4.5.3                 htmltools_0.5.9
```
