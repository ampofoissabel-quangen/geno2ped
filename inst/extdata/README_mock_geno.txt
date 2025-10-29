This package prototype expects an RDS at inst/extdata/mock_geno.rds
containing a list with elements:
- G: integer matrix (individuals x SNPs) coded 0/1/2 with NAs allowed
- ids: character vector of row names (length = nrow(G))
- map: data.frame with columns 'chr','pos','snp' (length = ncol(G))

You can create a toy file in R like:

set.seed(1)
n <- 50; p <- 2000
G <- matrix(sample(c(0,1,2,NA), n*p, TRUE, prob = c(0.49,0.02,0.49,0.00)), n, p)
ids <- paste0("ID", seq_len(n))
map <- data.frame(chr = sample(1:30, p, TRUE), pos = sample(1e6:5e6, p), snp = paste0("rs", seq_len(p)))
saveRDS(list(G=G, ids=ids, map=map), file="inst/extdata/mock_geno.rds")