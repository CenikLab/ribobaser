#' Select subset of samples for downstream processing
#'
#' Assumptions
#' 1) Inputs are genes by samples
#' 2) Both RNA and Ribo are calculated for same set of samples
#' 3) The data represents read counts and is not cpm/clr normalized.
#' 4) The Ribo and RNA matrices have undergone zero imputation.
#'
#' Existing features
#' a) Support removing samples where dummy genes are a large proportion of the
#' remaining counts
#' b) Support removing samples with poor correspondence between RNA and Ribo
#' Samples are selected for removal based on r^2 from lm(RIBO~RNA) after
#' CLR-normalization
#'
#' Optional features to Implement later
#' a) Support only ribo or RNA
#'
#' @param RIBO numeric matrix or data.frame with genes in rows and samples in columns
#' @param RNA numeric matrix or data.frame with genes in rows and samples in columns
#' @param high_dummy_percentage numeric between 0 and 1. Remove samples where `DUMMY_REMOVED`
#'   exceeds this fraction of the remaining counts in RIBO or RNA. Set to 0 to keep all samples.
#' @param min_r2 numeric between 0 and 1 specifying the minimum coefficient of determination
#'   (R-squared) required between CLR-normalised RIBO and RNA counts within each sample. Samples
#'   with `R_squared < min_r2` are removed. Set to 0 to skip this filter.
#' @param min_periodicity numeric between 0 and 1 specifying the minimum periodicity score
#'   required for a sample based on `Ribobase_QC_non_dedup_data$\`Periodicity score\``. Set to 0
#'   to keep all samples regardless of periodicity.
#'
#' @return list with elements:
#'   - ribo: filtered RIBO matrix (genes by samples)
#'   - rna:  filtered RNA matrix (genes by samples)
#' @examples
#' set.seed(3)
#' # two genes plus the dummy row, three samples
#' RIBO <- matrix(
#'   c(100, 120, 80,    # g1
#'     150, 130, 90,    # g2
#'      10, 200,  5),   # DUMMY_REMOVED
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("g1","g2","DUMMY_REMOVED"), paste0("s", 1:3))
#' )
#' RNA <- matrix(
#'   c(200, 210, 190,   # g1
#'     300, 290, 280,   # g2
#'      15, 250,  8),   # DUMMY_REMOVED
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("g1","g2","DUMMY_REMOVED"), paste0("s", 1:3))
#' )
#'
#' # Remove samples where DUMMY_REMOVED > 40 percent of remaining counts in either modality
#' res <- select_samples(RIBO, RNA, high_dummy_percentage = 0.4)
#' colnames(res$ribo)  # samples kept
#'
#' @importFrom compositions clr
#' @importFrom stats setNames
#' @export

select_samples <- function(RIBO, RNA,
                           high_dummy_percentage = 0,
                           min_r2 = 0,
                           min_periodicity = 0) {
  RIBO <- as.matrix(RIBO)
  RNA  <- as.matrix(RNA)

  if (!is.numeric(high_dummy_percentage) || length(high_dummy_percentage) != 1L ||
      is.na(high_dummy_percentage) || high_dummy_percentage < 0 || high_dummy_percentage > 1) {
    stop("high_dummy_percentage must be a number between 0 and 1")
  }

  if (!is.numeric(min_r2) || length(min_r2) != 1L ||
      is.na(min_r2) || min_r2 < 0 || min_r2 > 1) {
    stop("min_r2 must be a number between 0 and 1")
  }

  if (!is.numeric(min_periodicity) || length(min_periodicity) != 1L ||
      is.na(min_periodicity) || min_periodicity < 0 || min_periodicity > 1) {
    stop("min_periodicity must be a number between 0 and 1")
  }

  if (!identical(colnames(RIBO), colnames(RNA))) {
    stop("RIBO and RNA must share identical sample columns")
  }

  if (!identical(rownames(RIBO), rownames(RNA))) {
    stop("RIBO and RNA must share identical gene rows")
  }

  if(0 %in% RIBO || 0 %in% RNA) {
    stop("Ribo/RNA matrices contain zeroes. Please impute.")
  }

  # Early exit if no filtering requested
  if (high_dummy_percentage == 0 && min_r2 == 0 && min_periodicity == 0) {
    return(list(ribo = RIBO, rna = RNA))
  }

  # Track current view of matrices and apply filters sequentially
  current_ribo <- RIBO
  current_rna  <- RNA

  if (high_dummy_percentage > 0) {
    dummy_row <- "DUMMY_REMOVED"
    if (!(dummy_row %in% rownames(current_ribo) && dummy_row %in% rownames(current_rna))) {
      stop("DUMMY_REMOVED row must exist in both RIBO and RNA when high_dummy_percentage > 0")
    }

    ribo_dummy <- as.numeric(current_ribo[dummy_row, ])
    ribo_den   <- colSums(current_ribo[setdiff(rownames(current_ribo), dummy_row), , drop = FALSE])
    ribo_prop  <- ribo_dummy / ribo_den

    rna_dummy <- as.numeric(current_rna[dummy_row, ])
    rna_den   <- colSums(current_rna[setdiff(rownames(current_rna), dummy_row), , drop = FALSE])
    rna_prop  <- rna_dummy / rna_den

    remove_samples <- (ribo_prop > high_dummy_percentage) | (rna_prop > high_dummy_percentage)

    if (any(remove_samples)) {
      current_ribo <- current_ribo[, !remove_samples, drop = FALSE]
      current_rna  <- current_rna[,  !remove_samples, drop = FALSE]
    }
  }

  if (min_periodicity > 0 && ncol(current_ribo) > 0) {
    data("Ribobase_QC_non_dedup_data")
    meta_subset <- Ribobase_QC_non_dedup_data[!duplicated(Ribobase_QC_non_dedup_data$Experiment),
                            c("Experiment", "Periodicity score"),
                            drop = FALSE]
    periodicity_lookup <- stats::setNames(meta_subset$`Periodicity score`, meta_subset$Experiment)
    sample_periodicity <- periodicity_lookup[colnames(current_ribo)]

    if (any(is.na(sample_periodicity))) {
      missing <- colnames(current_ribo)[is.na(sample_periodicity)]
      stop("Missing periodicity scores for samples: ", paste(unique(missing), collapse = ", "))
    }

    keep_periodicity <- sample_periodicity >= min_periodicity
    if (any(!keep_periodicity)) {
      current_ribo <- current_ribo[, keep_periodicity, drop = FALSE]
      current_rna  <- current_rna[, keep_periodicity, drop = FALSE]
    }
  }

  if (min_r2 > 0 && ncol(current_ribo) > 0) {
    ribo_clr <- t(compositions::clr(t(current_ribo)))
    rna_clr <- t(compositions::clr(t(current_rna)))

    sample_r2 <- vapply(seq_len(ncol(ribo_clr)), function(i) {
      fit <- stats::lm(ribo_clr[, i] ~ rna_clr[, i])
      summary(fit)$r.squared
    }, numeric(1))

    low_correspondence <- is.na(sample_r2) | (sample_r2 < min_r2)
    if (any(low_correspondence)) {
      current_ribo <- current_ribo[, !low_correspondence, drop = FALSE]
      current_rna  <- current_rna[,  !low_correspondence, drop = FALSE]
    }
  }

  list(
    ribo = current_ribo,
    rna  = current_rna
  )
}
