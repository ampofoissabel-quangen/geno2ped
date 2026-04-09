#' Build a pedigree from SNP genotypes with tunable thresholds
#'
#' @param geno List from `read_genotypes()` containing G, ids, map, meta.
#' @param s_threshold Minimum pairwise similarity for candidate parent-offspring.
#' @param me_max Maximum Mendelian error rate to accept a trio.
#' @param oh_max Maximum opposing homozygotes per candidate.
#' @param use_sex If TRUE, restrict sires to Sex == "M" and dams to Sex == "F".
#' @param use_age If TRUE, require parents be older than offspring.
#' @param min_parent_gap Minimum years older than offspring.
#' @param top_k Keep top-k candidates per offspring per role before trio check.
#' @param preset One of "none", "high_precision", "high_coverage".
#' @param verbose Print progress.
#'
#' @return A list with pedigree, assignments, trios, summary, and settings.
#' @export
build_pedigree <- function(
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
) {
  validate_geno_object(geno)
  preset <- match.arg(preset)

  if (preset == "high_precision") {
    s_threshold <- 0.85
    me_max <- 0.002
    oh_max <- 10
    top_k <- 2
    use_sex <- TRUE
    use_age <- TRUE
  } else if (preset == "high_coverage") {
    s_threshold <- 0.70
    me_max <- 0.005
    oh_max <- 50
    top_k <- Inf
    use_sex <- TRUE
  }

  G <- geno$G
  ids <- geno$ids
  meta <- geno$meta

  if (is.null(meta)) {
    meta <- data.frame(ID = ids, stringsAsFactors = FALSE)
  }

  if (!"ID" %in% names(meta)) {
    meta$ID <- ids
  }

  meta <- meta[match(ids, meta$ID), , drop = FALSE]

  df <- data.frame(
    ID = ids,
    Sex = if ("Sex" %in% names(meta)) as.character(meta$Sex) else NA_character_,
    BirthYear = if ("BirthYear" %in% names(meta)) suppressWarnings(as.integer(meta$BirthYear)) else NA_integer_,
    stringsAsFactors = FALSE
  )

  sires <- ids
  dams  <- ids

  if (use_sex) {
    sires <- df$ID[!is.na(df$Sex) & tolower(df$Sex) %in% c("m", "male", "2")]
    dams  <- df$ID[!is.na(df$Sex) & tolower(df$Sex) %in% c("f", "female", "1")]
    if (!length(sires)) sires <- ids
    if (!length(dams)) dams <- ids
  }

  if (verbose) message("Computing pairwise candidate metrics ...")

  eval_pool <- function(off_ids, pool_ids, role) {
    out <- vector("list", length(off_ids))
    for (ii in seq_along(off_ids)) {
      o <- off_ids[ii]
      xo <- G[o, ]
      rows <- lapply(pool_ids, function(p) {
        xp <- G[p, ]
        n_shared <- shared_count(xo, xp)
        sim <- similarity_score(xo, xp)
        oh  <- oh_count(xo, xp)

        data.frame(
          Offspring = o,
          Parent = p,
          Role = role,
          similarity = sim,
          OH = oh,
          n_shared = n_shared,
          stringsAsFactors = FALSE
        )
      })
      out[[ii]] <- do.call(rbind, rows)
    }
    do.call(rbind, out)
  }

  cand_sire <- eval_pool(ids, sires, "Sire")
  cand_dam  <- eval_pool(ids, dams, "Dam")
  candidates <- rbind(cand_sire, cand_dam)

  candidates <- candidates[
    is.finite(candidates$similarity) &
      is.finite(candidates$OH) &
      candidates$Offspring != candidates$Parent,
    ,
    drop = FALSE
  ]

  candidates$passes_similarity <- candidates$similarity >= s_threshold
  candidates$passes_oh <- candidates$OH <= oh_max
  candidates$passes_pair_filter <- candidates$passes_similarity & candidates$passes_oh
  candidates$reject_reason <- ifelse(
    !candidates$passes_similarity, "low_similarity",
    ifelse(!candidates$passes_oh, "too_many_opposing_homozygotes", NA_character_)
  )

  candidates <- candidates[candidates$passes_pair_filter, , drop = FALSE]

  if (use_age) {
    by <- df$BirthYear
    names(by) <- df$ID

    keep <- !is.na(by[candidates$Offspring]) &
      !is.na(by[candidates$Parent]) &
      by[candidates$Parent] <= by[candidates$Offspring] - min_parent_gap

    candidates <- candidates[keep, , drop = FALSE]
  }

  if (is.finite(top_k)) {
    grp <- split(candidates, list(candidates$Offspring, candidates$Role), drop = TRUE)
    grp <- lapply(grp, function(d) {
      d <- d[order(-d$similarity, d$OH), , drop = FALSE]
      d[seq_len(min(nrow(d), top_k)), , drop = FALSE]
    })
    candidates <- do.call(rbind, grp)
    rownames(candidates) <- NULL
  }

  if (verbose) message("Validating trios ...")

  trio_rows <- list()
  idx <- 0L

  for (o in ids) {
    s_c <- candidates[candidates$Offspring == o & candidates$Role == "Sire", , drop = FALSE]
    d_c <- candidates[candidates$Offspring == o & candidates$Role == "Dam", , drop = FALSE]

    if (nrow(s_c) == 0 || nrow(d_c) == 0) next

    xo <- G[o, ]
    best <- NULL

    for (i in seq_len(nrow(s_c))) {
      xs <- G[s_c$Parent[i], ]
      for (j in seq_len(nrow(d_c))) {
        xd <- G[d_c$Parent[j], ]
        me <- me_rate_trio(xo, xs, xd)
        n_trio_shared <- sum(!is.na(xo) & !is.na(xs) & !is.na(xd))

        row <- data.frame(
          ID = o,
          Sire = s_c$Parent[i],
          Dam = d_c$Parent[j],
          s_sire = s_c$similarity[i],
          OH_sire = s_c$OH[i],
          s_dam = d_c$similarity[j],
          OH_dam = d_c$OH[j],
          trio_shared = n_trio_shared,
          ME = me,
          stringsAsFactors = FALSE
        )

        idx <- idx + 1L
        trio_rows[[idx]] <- row

        if (is.null(best) || (!is.na(me) && me < best$ME)) {
          best <- row
        }
      }
    }
  }

  trios_all <- if (length(trio_rows)) {
    do.call(rbind, trio_rows)
  } else {
    data.frame(
      ID = character(),
      Sire = character(),
      Dam = character(),
      s_sire = numeric(),
      OH_sire = numeric(),
      s_dam = numeric(),
      OH_dam = numeric(),
      trio_shared = integer(),
      ME = numeric(),
      stringsAsFactors = FALSE
    )
  }

  best_trios <- if (nrow(trios_all)) {
    split_trios <- split(trios_all, trios_all$ID)
    do.call(rbind, lapply(split_trios, function(d) d[which.min(d$ME), , drop = FALSE]))
  } else {
    trios_all
  }

  if (nrow(best_trios)) {
    best_trios$decision <- ifelse(
      is.finite(best_trios$ME) & best_trios$ME <= me_max,
      "accept",
      "reject"
    )
    best_trios$decision_reason <- ifelse(best_trios$decision == "accept", "passed", "high_trio_mendelian_error")
  } else {
    best_trios$decision <- character()
    best_trios$decision_reason <- character()
  }

  ped <- best_trios[best_trios$decision == "accept", c("ID", "Sire", "Dam"), drop = FALSE]
  all_ids <- data.frame(ID = ids, stringsAsFactors = FALSE)
  ped <- merge(all_ids, ped, by = "ID", all.x = TRUE)
  ped$Sire[is.na(ped$Sire)] <- 0
  ped$Dam[is.na(ped$Dam)] <- 0
  ped$assignment_status <- ifelse(ped$Sire != 0 & ped$Dam != 0, "complete", "unassigned")

  summary <- list(
    n_individuals = length(ids),
    n_complete_assignments = sum(ped$assignment_status == "complete"),
    assignment_rate = mean(ped$assignment_status == "complete"),
    mean_ME_accepted = if (any(best_trios$decision == "accept")) {
      mean(best_trios$ME[best_trios$decision == "accept"], na.rm = TRUE)
    } else {
      NA_real_
    }
  )

  settings <- list(
    preset = preset,
    s_threshold = s_threshold,
    me_max = me_max,
    oh_max = oh_max,
    use_sex = use_sex,
    use_age = use_age,
    min_parent_gap = min_parent_gap,
    top_k = top_k
  )

  list(
    pedigree = ped,
    assignments = candidates,
    trios = best_trios,
    trios_all = trios_all,
    summary = summary,
    settings = settings
  )
}
