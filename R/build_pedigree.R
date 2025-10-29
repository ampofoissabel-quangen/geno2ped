#' Build a pedigree from SNP genotypes with tunable thresholds
#'
#' @param geno List from `read_genotypes()` containing G (0/1/2), ids, map, meta.
#' @param s_threshold Minimum pairwise similarity for candidate PO (default 0.65).
#' @param me_max Trio Mendelian-error maximum to accept a trio (default 0.005).
#' @param oh_max Maximum opposing homozygotes per candidate (default Inf).
#' @param use_sex Logical; if TRUE, restrict sires to Sex=="M" and dams to "F".
#' @param use_age Logical; if TRUE, restrict parents to be older than offspring.
#' @param min_parent_gap Minimum years older than offspring (default 1).
#' @param top_k Keep only top-k candidates per offspring per sex (default Inf).
#' @param preset One of "none", "high_precision", or "high_coverage".
#' @param verbose Logical; print progress.
#' @return list(pedigree, trios, assignments)
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

  # apply preset values
  if (preset != "none") {
    if (preset == "high_precision") {
      s_threshold <- 0.85; me_max <- 0.002; oh_max <- 10; top_k <- 2
      use_sex <- TRUE; use_age <- TRUE
    }
    if (preset == "high_coverage") {
      s_threshold <- 0.70; me_max <- 0.005; oh_max <- 50; top_k <- Inf
      use_sex <- TRUE
    }
  }

  G <- geno$G
  ids <- geno$ids
  meta <- geno$meta %||% data.frame(ID = ids)
  if (is.null(rownames(G))) rownames(G) <- ids

  # Opposing homozygotes + similarity
  opp_hom <- function(x, y) sum((x == 0L & y == 2L) | (x == 2L & y == 0L), na.rm = TRUE)
  shared  <- function(x, y) sum(!is.na(x) & !is.na(y))
  sim_fun <- function(x, y) {
    m <- shared(x, y)
    if (m == 0) return(NA_real_)
    1 - opp_hom(x, y)/m
  }

  df <- data.frame(ID = ids, stringsAsFactors = FALSE)
  if (!is.null(meta$Sex)) df$Sex <- meta$Sex else df$Sex <- NA
  if (!is.null(meta$BirthYear)) df$BirthYear <- as.integer(meta$BirthYear)

  off <- df$ID
  sires <- df$ID; dams <- df$ID
  if (use_sex && "Sex" %in% names(df)) {
    sires <- df$ID[!is.na(df$Sex) & tolower(df$Sex) %in% c("m","male","2")]
    dams  <- df$ID[!is.na(df$Sex) & tolower(df$Sex) %in% c("f","female","1")]
  }

  # compute similarity & OH
  if (verbose) message("Computing similarities ...")
  to_tbl <- function(off, pool, ids_pool) {
    do.call(rbind, lapply(off, function(o) {
      xo <- G[o, , drop = FALSE]
      res <- sapply(ids_pool, function(p) {
        xp <- G[p, , drop = FALSE]
        m <- shared(xo, xp)
        if (m == 0) return(c(s = NA, oh = NA))
        oh <- opp_hom(xo, xp)
        s <- 1 - oh/m
        c(s = s, oh = oh)
      })
      data.frame(Offspring = o, Parent = ids_pool,
                 s = res["s",], OH = res["oh",], Role = pool,
                 stringsAsFactors = FALSE)
    }))
  }

  cand_sire <- to_tbl(off, "Sire", sires)
  cand_dam  <- to_tbl(off, "Dam",  dams)
  cand <- rbind(cand_sire, cand_dam)
  cand <- subset(cand, s >= s_threshold & OH <= oh_max)

  # optional age filter
  if (use_age && "BirthYear" %in% names(df)) {
    by <- df$BirthYear; names(by) <- df$ID
    cand <- subset(cand, !is.na(by[Offspring]) & !is.na(by[Parent]) &
                         by[Parent] <= by[Offspring] - min_parent_gap)
  }

  # top-k per offspring/role
  if (is.finite(top_k)) {
    cand <- cand |>
      split(list(cand$Offspring, cand$Role)) |>
      lapply(function(d) d[order(-d$s, d$OH), ][seq_len(min(nrow(d), top_k)), ]) |>
      do.call(rbind, .)
  }

  assignments <- cand

  # trio validation
  if (verbose) message("Validating trios ...")
  trios <- list()
  for (o in unique(assignments$Offspring)) {
    s_cand <- subset(assignments, Offspring == o & Role == "Sire")
    d_cand <- subset(assignments, Offspring == o & Role == "Dam")
    if (nrow(s_cand) == 0 || nrow(d_cand) == 0) next
    xo <- G[o, ]
    best <- NULL
    for (i in seq_len(nrow(s_cand))) {
      xs <- G[s_cand$Parent[i], ]
      for (j in seq_len(nrow(d_cand))) {
        xd <- G[d_cand$Parent[j], ]
        me_cnt <- 0L; m_cnt <- 0L
        for (k in seq_len(ncol(G))) {
          go <- xo[k]; gs <- xs[k]; gd <- xd[k]
          if (is.na(go) || is.na(gs) || is.na(gd)) next
          m_cnt <- m_cnt + 1L
          if ((gs == 0L && gd == 0L && go == 2L) ||
              (gs == 2L && gd == 2L && go == 0L)) me_cnt <- me_cnt + 1L
        }
        if (m_cnt == 0L) next
        me <- me_cnt / m_cnt
        row <- data.frame(
          ID = o, Sire = s_cand$Parent[i], Dam = d_cand$Parent[j],
          s_sire = s_cand$s[i], s_dam = d_cand$s[j], ME = me,
          stringsAsFactors = FALSE
        )
        if (is.null(best) || me < best$ME) best <- row
      }
    }
    if (!is.null(best)) trios[[length(trios)+1L]] <- best
  }
  trios <- if (length(trios)) do.call(rbind, trios) else
    data.frame(ID=character(),Sire=character(),Dam=character(),
               s_sire=numeric(),s_dam=numeric(),ME=numeric())

  trios$decision <- ifelse(trios$ME <= me_max, "accept", "reject")

  ped <- subset(trios, decision == "accept", select = c("ID","Sire","Dam"))
  all_ids <- data.frame(ID = ids)
  ped <- merge(all_ids, ped, by = "ID", all.x = TRUE)
  ped$Sire[is.na(ped$Sire)] <- 0
  ped$Dam[is.na(ped$Dam)]   <- 0

  list(pedigree = ped, trios = trios, assignments = assignments)
}

