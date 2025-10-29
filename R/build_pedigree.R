#' Build a pedigree from SNP genotypes with tunable thresholds
#'
#' @param geno List from `read_genotypes()` containing G (0/1/2), ids, map, meta.
#' @param s_threshold Minimum pairwise similarity for candidate PO (default 0.65).
#' @param me_max Trio Mendelian-error maximum to accept a trio (default 0.005).
#' @param oh_max Optional cap on opposing homozygotes for PO candidates (default Inf).
#' @param use_sex If TRUE, restrict sire candidates to Sex=="M" and dams to "F" (if available).
#' @param use_age If TRUE, require parent BirthYear <= offspring BirthYear - min_parent_gap.
#' @param min_parent_gap Minimum years older than offspring for parents (default 1).
#' @param top_k Candidates to keep per offspring per sex before ME check (default Inf = all).
#' @param preset Character; one of c("none","high_precision","high_coverage"). Applies recommended thresholds unless overridden by explicit args.
#' @param verbose Logical; print progress.
#'
#' @return list with elements:
#'   \describe{
#'     \item{pedigree}{data.frame(ID, Sire, Dam)}
#'     \item{trios}{data.frame with trio-level QC (ID, Sire, Dam, ME, decision, s_sire, s_dam, OH_sire, OH_dam)}
#'     \item{assignments}{candidate-level table before trio validation}
#'   }
#' @export
build_pedigree <- function(
  geno,
  s_threshold = 0.65,
  me_max      = 0.005,
  oh_max      = Inf,
  use_sex     = TRUE,
  use_age     = FALSE,
  min_parent_gap = 1,
  top_k       = Inf,
  preset      = c("none", "high_precision", "high_coverage"),
  verbose     = TRUE
) {
  preset <- match.arg(preset)

  # Apply presets unless user explicitly changed args from defaults
  if (preset != "none") {
    if (preset == "high_precision") {
      if (missing(s_threshold) || identical(s_threshold, 0.65)) s_threshold <- 0.85
      if (missing(me_max)      || identical(me_max, 0.005))      me_max      <- 0.002
      if (missing(oh_max)      || is.infinite(oh_max))           oh_max      <- 10
      if (missing(top_k)       || is.infinite(top_k))            top_k       <- 2
      use_sex <- TRUE
      use_age <- TRUE
    } else if (preset == "high_coverage") {
      if (missing(s_threshold) || identical(s_threshold, 0.65)) s_threshold <- 0.70
      if (missing(me_max)      || identical(me_max, 0.005))      me_max      <- 0.005
      if (missing(oh_max)      || is.infinite(oh_max))           oh_max      <- 50
      if (missing(top_k)       || is.infinite(top_k))            top_k       <- Inf
      use_sex <- TRUE
      # age optional for coverage
    }
  }

  G   <- geno$G
  ids <- geno$ids
  meta <- geno$meta %||% data.frame(ID = ids)

  stopifnot(length(ids) == nrow(G))
  if (is.null(rownames(G))) rownames(G) <- ids

  # helper: opposing homozygotes & similarity
  opp_hom <- function(x, y) sum((x == 0L & y == 2L) | (x == 2L & y == 0L), na.rm = TRUE)
  shared  <- function(x, y) sum(!is.na(x) & !is.na(y))
  sim_fun <- function(x, y) {
    m <- shared(x, y); if (m == 0) return(NA_real_)
    1 - opp_hom(x, y) / m
  }

  # derive roles (simple heuristic from Sex if present)
  df <- data.frame(ID = ids, stringsAsFactors = FALSE)
  if (!is.null(meta$Sex)) df$Sex <- as.character(meta$Sex) else df$Sex <- NA_character_
  if (!is.null(meta$BirthYear)) df$BirthYear <- suppressWarnings(as.integer(meta$BirthYear)) else df$BirthYear <- NA_integer_

  # offspring = all individuals by default
  off <- df$ID

  # candidate pools
  sires <- df$ID
  dams  <- df$ID
  if (use_sex && "Sex" %in% names(df)) {
    sires <- df$ID[!is.na(df$Sex) & tolower(df$Sex) %in% c("m","male","1","2")] # tolerant input
    dams  <- df$ID[!is.na(df$Sex) & tolower(df$Sex) %in% c("f","female","1","2")]
    # if coded 1/2 but unknown mapping, we still allow both pools to avoid false negatives
  }

  # precompute similarities & OH against all candidates (vectorized by apply for clarity)
  if (verbose) message("Computing similarities...")
  sim_sire <- lapply(off, function(o) {
    xo <- G[match(o, ids), ]
    sapply(sires, function(s) {
      xs <- G[match(s, ids), ]
      m <- shared(xo, xs); if (m == 0) return(c(s = NA_real_, oh = NA_integer_))
      oh <- opp_hom(xo, xs); s <- 1 - oh/m
      c(s = s, oh = oh)
    })
  })
  names(sim_sire) <- off

  sim_dam <- lapply(off, function(o) {
    xo <- G[match(o, ids), ]
    sapply(dams, function(d) {
      xd <- G[match(d, ids), ]
      m <- shared(xo, xd); if (m == 0) return(c(s = NA_real_, oh = NA_integer_))
      oh <- opp_hom(xo, xd); s <- 1 - oh/m
      c(s = s, oh = oh)
    })
  })
  names(sim_dam) <- off

  # to long tables
  to_tbl <- function(lst, pool) {
    do.call(rbind, lapply(names(lst), function(o) {
      mat <- t(lst[[o]])
      data.frame(Offspring = o, Parent = rownames(mat),
                 s = as.numeric(mat[, "s"]), OH = as.numeric(mat[, "oh"]),
                 Role = pool, stringsAsFactors = FALSE)
    }))
  }
  cand_sire <- to_tbl(sim_sire, "Sire")
  cand_dam  <- to_tbl(sim_dam,  "Dam")

  cand <- rbind(cand_sire, cand_dam)
  cand <- cand[is.finite(cand$s) & is.finite(cand$OH), , drop = FALSE]
  cand <- subset(cand, s >= s_threshold & OH <= oh_max)

  # optional age filter
  if (use_age && "BirthYear" %in% names(df)) {
    by <- df$BirthYear; names(by) <- df$ID
    cand <- subset(cand, !is.na(by[Offspring]) & !is.na(by[Parent]) &
                         by[Parent] <= by[Offspring] - min_parent_gap)
  }

  # top-k pruning per offspring/role
  if (is.finite(top_k)) {
    cand <- cand |>
      split(list(cand$Offspring, cand$Role)) |>
      lapply(function(d) d[order(-d$s, d$OH), ][seq_len(min(nrow(d), top_k)), ]) |>
      do.call(rbind, .)
  }

  assignments <- cand

  # trio validation: for each Offspring choose best sire x dam pair by lowest ME; accept if <= me_max
  if (verbose) message("Validating trios...")
  tri_rows <- list()
  for (o in unique(assignments$Offspring)) {
    s_cand <- subset(assignments, Offspring == o & Role == "Sire")
    d_cand <- subset(assignments, Offspring == o & Role == "Dam")
    if (nrow(s_cand) == 0 || nrow(d_cand) == 0) next
    xo <- G[match(o, ids), ]
    best <- NULL
    for (i in seq_len(nrow(s_cand))) {
      xs <- G[match(s_cand$Parent[i], ids), ]
      for (j in seq_len(nrow(d_cand))) {
        xd <- G[match(d_cand$Parent[j], ids), ]
        # Mendelian set consistency
        me_cnt <- 0L; m_cnt <- 0L
        for (k in seq_len(ncol(G))) {
          go <- xo[k]; gs <- xs[k]; gd <- xd[k]
          if (is.na(go) || is.na(gs) || is.na(gd)) next
          m_cnt <- m_cnt + 1L
          # incompatible cases
          if ((gs == 0L && gd == 0L && go == 2L) ||
              (gs == 2L && gd == 2L && go == 0L)) {
            me_cnt <- me_cnt + 1L
          }
        }
        if (m_cnt == 0L) next
        me <- me_cnt / m_cnt
        row <- data.frame(
          ID = o,
          Sire = s_cand$Parent[i],
          Dam  = d_cand$Parent[j],
          s_sire = s_cand$s[i],
          OH_sire = s_cand$OH[i],
          s_dam  = d_cand$s[j],
          OH_dam = d_cand$OH[j],
          ME = me,
          stringsAsFactors = FALSE
        )
        if (is.null(best) || me < best$ME) best <- row
      }
    }
    if (!is.null(best)) tri_rows[[length(tri_rows) + 1L]] <- best
  }

  trios <- if (length(tri_rows)) do.call(rbind, tri_rows) else
    data.frame(ID=character(), Sire=character(), Dam=character(),
               s_sire=numeric(), OH_sire=numeric(), s_dam=numeric(),
               OH_dam=numeric(), ME=numeric(), stringsAsFactors = FALSE)

  trios$decision <- ifelse(is.finite(trios$ME) & trios$ME <= me_max, "accept",
                           ifelse(is.finite(trios$ME) & trios$ME <= me_max*1.5, "marginal", "reject"))

  ped <- subset(trios, decision == "accept", select = c("ID","Sire","Dam"))
  # ensure every ID appears (fill missing parents with 0)
  all_ids <- data.frame(ID = ids, stringsAsFactors = FALSE)
  ped <- merge(all_ids, ped, by = "ID", all.x = TRUE)
  ped$Sire[is.na(ped$Sire)] <- 0
  ped$Dam[is.na(ped$Dam)]   <- 0

  list(pedigree = ped, trios = trios, assignments = assignments)
}

