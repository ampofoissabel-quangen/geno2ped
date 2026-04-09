library(geno2ped)

geno_file <- system.file("extdata", "toy_genotypes.csv", package = "geno2ped")
meta_file <- system.file("extdata", "toy_metadata.csv", package = "geno2ped")

geno <- read_genotypes(geno_file, sample_metadata = meta_file)

res <- build_pedigree(
  geno,
  preset = "high_precision",
  verbose = TRUE
)

print(res$pedigree)
print(res$summary)

plot_kinship(res)
plot_me(res)
plot_assignment_status(res)

write_outputs(res, outdir = "toy_results", prefix = "toy")
