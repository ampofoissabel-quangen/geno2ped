#' Build pedigree from genotypes
#'
#' @param gobj Output from read_genotypes()
#' @param error_rate Allowed genotyping error rate (default 0.005)
#' @param max_candidates Maximum parent candidates to retain per calf
#' @return A list with pedigree data.frame and QC tables
#' @export
build_pedigree <- function(gobj, error_rate = 0.005, max_candidates = 10) {
  if (is.null(gobj$G)) stop("gobj must contain $G")
  G <- gobj$G; ids <- gobj$ids; meta <- gobj$meta
  n <- length(ids)
  # Simple IBS-based similarity as a stand-in for kinship (prototype)
  # Similarity = 1 - mean absolute genotype difference / 2
  sim <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in i:n) {
      gi <- G[i, ]; gj <- G[j, ]
      ok <- !is.na(gi) & !is.na(gj)
      if (!any(ok)) {
        s <- NA_real_
      } else {
        s <- 1 - mean(abs(gi[ok] - gj[ok]))/2
      }
      sim[i, j] <- sim[j, i] <- s
    }
  }

  # Heuristic: parent-offspring likely in a similarity band (e.g., 0.6-0.8 for mock data)
  # In real package, replace with KING-kinship and IBD profile.
  candidates <- lapply(seq_len(n), function(ci){
    s <- sim[ci, ]
    ord <- order(s, decreasing = TRUE, na.last = TRUE)
    ord <- ord[ord != ci][1:min(max_candidates, length(ord)-1)]
    data.frame(child = ids[ci], cand = ids[ord], sim = s[ord], stringsAsFactors = FALSE)
  })
  cand_df <- do.call(rbind, candidates)

  # Duo scoring by opposing homozygotes (OH)
  oh <- function(a, b) oh_count(G[which(ids==a), ], G[which(ids==b), ])
  cand_df$oh <- mapply(oh, cand_df$child, cand_df$cand)

  # Pick best two parents per child (just top-2 by similarity + low OH)
  chosen <- do.call(rbind, lapply(split(cand_df, cand_df$child), function(dd){
    dd <- dd[order(dd$oh, dd$sim, decreasing = c(FALSE, TRUE)), ]
    head(dd, 2)
  }))

  # Trios if two parents chosen
  trios <- do.call(rbind, lapply(split(chosen, chosen$child), function(dd){
    if (nrow(dd) < 2) return(NULL)
    sire <- dd$cand[1]; dam <- dd$cand[2]
    me <- me_rate_trio(G[ids==dd$child[1], ], G[ids==sire, ], G[ids==dam, ])
    data.frame(child = dd$child[1], sire = sire, dam = dam, ME = me, stringsAsFactors = FALSE)
  }))
  if (is.null(trios)) trios <- data.frame(child=character(), sire=character(), dam=character(), ME=numeric())

  # Decision: accept if ME <= 0.015 (3x error_rate just for prototype)
  trios$decision <- ifelse(!is.na(trios$ME) & trios$ME <= 3*error_rate, "pass",
                           ifelse(!is.na(trios$ME) & trios$ME <= 6*error_rate, "marginal", "fail"))

  pedigree <- data.frame(ID = ids, Sire = NA_character_, Dam = NA_character_, stringsAsFactors = FALSE)
  rownames(pedigree) <- pedigree$ID
  if (nrow(trios)) {
    for (i in seq_len(nrow(trios))) {
      ch <- trios$child[i]
      if (trios$decision[i] != "fail") {
        pedigree[ch, "Sire"] <- trios$sire[i]
        pedigree[ch, "Dam"]  <- trios$dam[i]
      }
    }
  }

  list(
    pedigree = unname(pedigree),
    candidates = cand_df,
    trios = trios,
    similarity = sim
  )
}

#' Write standard outputs to a directory
#' @param res Object from build_pedigree()
#' @param out_dir Directory
#' @export
write_outputs <- function(res, out_dir = "geno2ped_out") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  utils::write.table(res$pedigree, file = file.path(out_dir, "pedigree.csv"),
                     sep = ",", row.names = FALSE, quote = TRUE)
  utils::write.table(res$candidates, file = file.path(out_dir, "assignments.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(res$trios, file = file.path(out_dir, "trios.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
}