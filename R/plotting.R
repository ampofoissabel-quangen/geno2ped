#' Plot kinship/similarity histogram
#' @param res Result of build_pedigree()
#' @export
plot_kinship <- function(res) {
  package_assert("ggplot2")
  m <- res$similarity
  v <- m[upper.tri(m)]
  df <- data.frame(sim = v)
  p <- ggplot2::ggplot(df, ggplot2::aes(sim)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::labs(x = "Pairwise similarity (proxy for kinship)", y = "Count")
  p
}

#' Plot Mendelian error rate histogram for accepted trios
#' @param res Result of build_pedigree()
#' @export
plot_me <- function(res) {
  package_assert("ggplot2")
  df <- res$trios
  df <- df[is.finite(df$ME), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(ME)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::labs(x = "Mendelian error rate", y = "Trios")
  p
}