#' decon.cell.props
#'
#' Call this function to estimate cell proportions given hetergenous expressions and cell signatures.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param het.sub  Should be samples by genomic site (rows by columns).  The samples by genomic measures heterogeneous observations matrix.
#' @param cell.sigs  Should be cell types by genes (rows by columns).
#' @param LOGGED  A boolean indicating whether the heterogeneous observations have been log transformed.
#' @return cell.props The deconvolved (estimated cell proportions)
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details
#' decon.cell.props is an auxillary function meant to be called directly by user in order to estimate cell proportions given generic cell gene signatures.
#' This assumes that G and cell.sigs have been previously made conformal with identical gene names/labels.
#' @export

decon.cell.props = function( het.sub, cell.sigs , LOGGED ) {

###############################################################################
if( LOGGED ) {
  # Convert signatures to log2 space:
  cell.sigs = log2( cell.sigs + 1 )
}
###############################################################################

###############################################################################
# Now deconvolve cell proportions (should have an sample by cell type matrix):
# Orient the het.sub and cell.sigs matrices correctly: het.sub is samples by genes.  cell.sigs is cells by genes:
fit = lm( t( het.sub ) ~ 0 + t( cell.sigs ) )  # Through the origin.  Otherwise there is an intercept term with what biological meaning?
cell.props = t( fit$coefficients )
###############################################################################

return( cell.props )

} # End function
