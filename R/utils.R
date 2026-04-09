# R/utils.R

#' @noRd
validate_geno_object <- function(geno) {
  if (!is.list(geno))
    stop("'geno' must be a list returned by read_genotypes().", call. = FALSE)
  if (is.null(geno$G) || !is.matrix(geno$G))
    stop("'geno$G' must be a numeric matrix.", call. = FALSE)
  if (is.null(geno$ids) || !is.character(geno$ids))
    stop("'geno$ids' must be a character vector.", call. = FALSE)
  if (is.null(geno$map) || !is.data.frame(geno$map))
    stop("'geno$map' must be a data.frame.", call. = FALSE)
  invisible(TRUE)
}

#' @noRd
shared_count <- function(x, y) {
  sum(!is.na(x) & !is.na(y))
}

#' @noRd
similarity_score <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) == 0) return(NA_real_)
  1 - mean(abs(x[keep] - y[keep]) / 2)
}

#' @noRd
oh_count <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  sum((x[keep] == 0 & y[keep] == 2) | (x[keep] == 2 & y[keep] == 0))
}

#' @noRd
me_rate_trio <- function(xo, xs, xd) {
  keep <- !is.na(xo) & !is.na(xs) & !is.na(xd)
  if (sum(keep) == 0) return(NA_real_)
  xo <- xo[keep]; xs <- xs[keep]; xd <- xd[keep]
  errors <- mapply(function(o, s, d) {
    possible <- outer(c(s %/% 2, s %% 2 + s %/% 2),
                      c(d %/% 2, d %% 2 + d %/% 2), "+")
    !any(possible == o)
  }, xo, xs, xd)
  mean(errors)
}

#' @noRd
package_assert <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install with: install.packages('%s')", pkg, pkg),
         call. = FALSE)
  }
}
