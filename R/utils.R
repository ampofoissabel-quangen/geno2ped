#' @keywords internal
package_assert <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}

#' Internal: Opposing homozygotes count between two genotype vectors
#' @param g1 Integer vector coded 0/1/2 (AA/AB/BB) with NA allowed
#' @param g2 Integer vector coded 0/1/2 (AA/AB/BB) with NA allowed
#' @return Integer count
#' @keywords internal
oh_count <- function(g1, g2) {
  ok <- !is.na(g1) & !is.na(g2)
  sum( (g1[ok] == 0 & g2[ok] == 2) | (g1[ok] == 2 & g2[ok] == 0) )
}

#' Internal: Mendelian error rate for a trio
#' @param child,sire,dam integer genotype vectors (0/1/2) with NA allowed
#' @return numeric ME rate in [0,1]
#' @keywords internal
me_rate_trio <- function(child, sire, dam) {
  ok <- !is.na(child) & !is.na(sire) & !is.na(dam)
  if (!any(ok)) return(NA_real_)
  c <- child[ok]; s <- sire[ok]; d <- dam[ok]
  # Allowed child genotypes given parents (simplified)
  bad <- rep(FALSE, length(c))
  for (i in seq_along(c)) {
    ps <- sort(c(s[i], d[i]))
    allowed <- switch(paste(ps, collapse = ""),
      "00" = c(0,1),
      "02" = c(0,1,2),
      "22" = c(1,2),
      "01" = c(0,1,2),
      "12" = c(0,1,2),
      "11" = c(0,1,2),
      c(0,1,2)
    )
    bad[i] <- !(c[i] %in% allowed)
  }
  mean(bad)
}