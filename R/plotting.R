#' Plot candidate-parent similarity distribution
#'
#' @param res Result of build_pedigree()
#' @export
plot_kinship <- function(res) {
  package_assert("ggplot2")

  df <- res$assignments
  if (is.null(df) || !nrow(df)) {
    stop("No assignment data found in result object.", call. = FALSE)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = similarity, fill = Role)) +
    ggplot2::geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
    ggplot2::labs(
      x = "Candidate parent-offspring similarity",
      y = "Count",
      title = "Distribution of candidate-parent similarity scores"
    ) +
    ggplot2::theme_minimal()
}

#' Plot Mendelian error rate histogram for tested trios
#'
#' @param res Result of build_pedigree()
#' @export
plot_me <- function(res) {
  package_assert("ggplot2")

  df <- res$trios_all
  if (is.null(df) || !nrow(df)) {
    stop("No trio validation data found in result object.", call. = FALSE)
  }

  df <- df[is.finite(df$ME), , drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(x = ME)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::labs(
      x = "Mendelian error rate",
      y = "Trios",
      title = "Distribution of trio Mendelian error rates"
    ) +
    ggplot2::theme_minimal()
}

#' Plot assignment status counts
#'
#' @param res Result of build_pedigree()
#' @export
plot_assignment_status <- function(res) {
  package_assert("ggplot2")

  df <- res$pedigree
  if (is.null(df) || !nrow(df)) {
    stop("No pedigree data found in result object.", call. = FALSE)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = assignment_status)) +
    ggplot2::geom_bar() +
    ggplot2::labs(
      x = "Assignment status",
      y = "Count",
      title = "Pedigree assignment outcomes"
    ) +
    ggplot2::theme_minimal()
}
