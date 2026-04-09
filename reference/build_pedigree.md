# Build a pedigree from SNP genotypes with tunable thresholds

Build a pedigree from SNP genotypes with tunable thresholds

## Usage

``` r
build_pedigree(
  geno,
  s_threshold = 0.65,
  me_max = 0.005,
  oh_max = Inf,
  use_sex = TRUE,
  use_age = FALSE,
  min_parent_gap = 1,
  top_k = Inf,
  preset = c("none", "high_precision", "high_coverage"),
  verbose = TRUE
)
```

## Arguments

- geno:

  List from
  [`read_genotypes()`](https://ampofoissabel-quangen.github.io/geno2ped/reference/read_genotypes.md)
  containing G, ids, map, meta.

- s_threshold:

  Minimum pairwise similarity for candidate parent-offspring.

- me_max:

  Maximum Mendelian error rate to accept a trio.

- oh_max:

  Maximum opposing homozygotes per candidate.

- use_sex:

  If TRUE, restrict sires to Sex == "M" and dams to Sex == "F".

- use_age:

  If TRUE, require parents be older than offspring.

- min_parent_gap:

  Minimum years older than offspring.

- top_k:

  Keep top-k candidates per offspring per role before trio check.

- preset:

  One of "none", "high_precision", "high_coverage".

- verbose:

  Print progress.

## Value

A list with pedigree, assignments, trios, summary, and settings.
