# One-page PDF with kinship-like similarity and ME% histograms (using mock values)
pdf("inst/paper/figures_geno2ped_mock.pdf", width = 8.5, height = 11)

par(mfrow=c(2,1), mar=c(5,5,2,1))

set.seed(1)
sim <- pmin(pmax(rnorm(500, mean=0.7, sd=0.08), 0), 1)
hist(sim, breaks=40, main="", xlab="Pairwise similarity (proxy for kinship)", ylab="Count")

me <- pmax(rnorm(200, mean=0.004, sd=0.002), 0)
hist(me, breaks=30, main="", xlab="Mendelian error rate", ylab="Trios")

dev.off()
cat("Wrote inst/paper/figures_geno2ped_mock.pdf\n")