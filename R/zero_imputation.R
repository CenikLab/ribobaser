#' Replace zero counts with pseudocounts
#'
#' Many compositional workflows require strictly positive counts before
#' applying log-ratio transforms. This utility either performs a simple
#' zero-imputation step (replace zeros with the minimum non-zero value) or
#' a multiplicative replacement using `zCompositions::cmultRepl()`.
#'
#' @param RIBO numeric matrix or data.frame of counts with genes in rows and
#'   samples in columns.
#' @param RNA numeric matrix or data.frame of counts with genes in rows and
#'   samples in columns.
#' @param method imputation strategy to use.
#'   - `"simple"` replaces zeros with either the matrix minimum value if it's <1
#'       or 1
#'   - `"multiplicative"` delegates to `zCompositions::cmultRepl()`.
#' @param multiplicative_method method argument passed to
#'   `zCompositions::cmultRepl()` when `method = "multiplicative"`.
#' @param output output argument passed to `zCompositions::cmultRepl()` when
#'   `method = "multiplicative"`.
#' @param transpose whether to transpose matrices before multiplicative
#'   imputation. `cmultRepl()` expects samples in rows and features in columns,
#'   so the default `TRUE` works with matrices that store genes in rows.
#'
#' @return list with elements
#'   - ribo: `RIBO` converted to a matrix with imputed values.
#'   - rna:  `RNA` converted to a matrix with imputed values.
#'
#' Samples or genes removed internally during `cmultRepl()` are realigned by
#' keeping the common rows and columns between the imputed RIBO and RNA matrices.
#'
#' @examples
#' RIBO <- matrix(c(0, 5, 10,
#'                  3, 0,  7),
#'                nrow = 2,
#'                dimnames = list(c("g1", "g2"), paste0("s", 1:3)))
#' RNA  <- matrix(c(2, 0, 4,
#'                  0, 6, 8),
#'                nrow = 2,
#'                dimnames = list(c("g1", "g2"), paste0("s", 1:3)))
#'
#' zero_imputation(RIBO, RNA)
#' zero_imputation(RIBO, RNA, method = "multiplicative")
#'
#' @export
zero_imputation <- function(
    RIBO,
    RNA,
    method = c("simple", "multiplicative"),
    multiplicative_method = "GBM",
    output = "p-count",
    transpose = TRUE) {
  method <- match.arg(method)

  ribo_mat <- as.matrix(RIBO)
  rna_mat  <- as.matrix(RNA)

  if (method == "simple") {
    ribo_min <- min(min(ribo_mat[ribo_mat > 0]), 1)
    rna_min <- min(min(rna_mat[rna_mat > 0]), 1)

    ribo_mat[ribo_mat == 0] <- ribo_min
    rna_mat[rna_mat == 0] <- rna_min

  } else {
    impute_multiplicative <- function(mat) {
      target <- if (transpose) t(mat) else mat
      replaced <- zCompositions::cmultRepl(
        target,
        method = multiplicative_method,
        output = output
      )
      if (transpose) {
        replaced <- t(replaced)
      }
      replaced
    }

    ribo_mat <- impute_multiplicative(ribo_mat)
    rna_mat  <- impute_multiplicative(rna_mat)
  }

  common_genes   <- intersect(rownames(ribo_mat), rownames(rna_mat))
  common_samples <- intersect(colnames(ribo_mat), colnames(rna_mat))

  ribo_mat <- ribo_mat[common_genes, common_samples, drop = FALSE]
  rna_mat  <- rna_mat[common_genes, common_samples, drop = FALSE]

  list(
    ribo = ribo_mat,
    rna  = rna_mat
  )
}
