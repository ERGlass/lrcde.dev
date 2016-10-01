#' decon.cell.props
#'
#' Function to estimate cell proportions given hetergenous gene expression matrix and 
#' cell type-specific gene expression signatures.
#' @param het.sub heterogeneous gene expression matrix, _n_ samples by _m_ genes 
#' (rows by columns). Required. Note that only a subset of heterogeneous measres
#' is used, with the same number and order of genes as in the cell signature matrix
#' @param cell.sigs cell-specific gene expression signatures, _p_ cell types by 
#' _m_ genes (rows by columns). Required. Note that both heterogeneous matrix and 
#' cell signature matrix should be on the same scale, e.g. log2 or decimal
#' @return estimated cell proportions, _n_ samples by _m_ cell types (rows by columns)
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @keywords Deconvolution cell type-specific differential expression detection power analysis

#' @export
#' @examples 
#' 
#' \dontrun{
#' # Loads the heterogeneous matrix (subset) and cell signatures
#' data("lrcde.demo") # Included in the LRCDE package
#' # Estimate cell proportions. 
#' cell.props <- decon.cell.props(het.mtx, cell.sigs)
#' }

decon.cell.props <- function(het.sub, cell.sigs) {
  # Deconvolve cell proportions (should have an sample by cell type matrix):
  # Orient the het.sub and cell.sigs matrices correctly: het.sub is genes by samples.
  # cell.sigs is genes by cells
  # Force linear regression through the origin - zero cell proportions should contribute zero gene expression.
  fit = lm(t(het.sub) ~ 0 + t(cell.sigs))  
  cell.props = t(fit$coefficients)
  return( cell.props )
}
