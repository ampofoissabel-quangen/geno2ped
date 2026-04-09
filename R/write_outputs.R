#' Write geno2ped outputs to disk
#'
#' @param res Result from `build_pedigree()`
#' @param outdir Output directory
#' @param prefix Optional file prefix
#'
#' @return Invisibly returns a character vector of written file paths
#' @export
write_outputs <- function(res, outdir = "geno2ped_results", prefix = NULL) {
  if (is.null(res$pedigree) || is.null(res$assignments) || is.null(res$trios)) {
    stop("Result object does not look like output from build_pedigree().", call. = FALSE)
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  pref <- if (!is.null(prefix) && nzchar(prefix)) paste0(prefix, "_") else ""

  pedigree_file    <- file.path(outdir, paste0(pref, "pedigree.csv"))
  assignments_file <- file.path(outdir, paste0(pref, "assignments.csv"))
  trios_file       <- file.path(outdir, paste0(pref, "trios.csv"))
  trios_all_file   <- file.path(outdir, paste0(pref, "trios_all.csv"))
  summary_file     <- file.path(outdir, paste0(pref, "summary.txt"))
  settings_file    <- file.path(outdir, paste0(pref, "settings.txt"))

  utils::write.csv(res$pedigree, pedigree_file, row.names = FALSE)
  utils::write.csv(res$assignments, assignments_file, row.names = FALSE)
  utils::write.csv(res$trios, trios_file, row.names = FALSE)

  if (!is.null(res$trios_all)) {
    utils::write.csv(res$trios_all, trios_all_file, row.names = FALSE)
  }

  con <- file(summary_file, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("geno2ped summary", con)
  writeLines("================", con)
  writeLines("", con)

  if (!is.null(res$summary)) {
    for (nm in names(res$summary)) {
      writeLines(paste0(nm, ": ", res$summary[[nm]]), con)
    }
  }

  writeLines("", con)
  writeLines("settings", con)
  writeLines("========", con)

  if (!is.null(res$settings)) {
    for (nm in names(res$settings)) {
      writeLines(paste0(nm, ": ", res$settings[[nm]]), con)
    }
  }

  if (!is.null(res$settings)) {
    settings_df <- data.frame(
      parameter = names(res$settings),
      value = unlist(res$settings),
      stringsAsFactors = FALSE
    )
    utils::write.table(
      settings_df,
      settings_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  files_written <- c(
    pedigree_file,
    assignments_file,
    trios_file,
    if (!is.null(res$trios_all)) trios_all_file,
    summary_file,
    settings_file
  )

  invisible(files_written)
}
