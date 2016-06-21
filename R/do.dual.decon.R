#' do.dual.decon
#'
#' Performs the core functionality of LRCDE.  This function does the actual linear regression deconvolution per group.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param het.sub The samples by genomic measures heterogeneous observations matrix
#' @param cell.props  Should be samples by cell types (rows by columns).  The relative cell proportions per sample.
#' @param groups Vector of 1's and 2's indicating group membership (controls and cases).
#' @param nonNeg Boolean indicating whether to force cell type-specific expression estiamtes to be non-negative.
#' @param GET.FOLDS Boolean idicating whether to take the absolute difference between group-wise cell type-specific estimate or to take ratio between cases and controls (for fold change estimate).
#' @return Returns a list containing group-wise difference estimates, regression residuals, lm objects per group, and standard errors of cell type-specific expression estimates.  Also returns a two item list containing results of Breusch-Pagan test.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details
#' Here is where the actual group-wise linear regressions are performed.
#' This function is designed to handle a single genomic site at a time; e.g.: a single vector containing heterogeneous observations for controls and cases.
#' This is the core functionality of linear regression cell type-specific differential expression detection.
#' @export
#' @examples
#' do.dual.decon( het.obs, cell.props, groups, medCntr, stdz, nonNeg, GET.FOLDS=FALSE )

do.dual.decon <- function( het, cell.props, groups, medCntr, stdz, nonNeg , GET.FOLDS=FALSE )
{

  ###########################################################################
  unique.groups = unique( groups )
  n = group.wise.sample.size( groups )
  n1 = n[1]
  n2 = n[2]
  ###########################################################################

  numgene  = ncol( het )
  numcell  = ncol( cell.props )
  deconv <- list()
  ########################################
  # Group-wise linear regression:
  for ( curset in unique.groups ){ # print(curset)
    deconv[[ curset ]] = lm( het[ groups==curset, ] ~ 0 + cell.props[ groups==curset, ] ) # through the origin
  }
  ########################################

  # Cell type specific expression estimates per group:
  cell.expr.ests.1 = deconv[[ 1 ]]$coefficients   # Cell-type specific expression estimates from group 1
  cell.expr.ests.2 = deconv[[ 2 ]]$coefficients   # Cell-type specific expression estimates from group 2

  # We use the residuals to calculate MSE per group per regression:
  resids1 = deconv[[1]]$residuals
  resids2 = deconv[[2]]$residuals
  resids = rbind( resids1, resids2 )

  se1 = ls.diag(deconv[[ 1 ]])$std.err           # std err of Cell-type specific expression estimates group 1
  se2 = ls.diag(deconv[[ 2 ]])$std.err           # std err of Cell-type specific expression estimates group 2

  if (nonNeg) { # Force cell-type specific expression estimates to be non-negative (and non-zero if we're estimating fold change)?
    cell.expr.ests.1[cell.expr.ests.1 < 0] = 0.000000001 # Non-zero avoids division by zero if GET.FOLDS
    cell.expr.ests.2[cell.expr.ests.2 < 0] = 0.000000001
  }

  # Get group difference / fold.x estimates
  fold.diff.ests  = array( dim = c( numcell, numgene ))

  if( GET.FOLDS ) {
    fold.diff.ests = cell.expr.ests.2 / cell.expr.ests.1   # fold.x ests
  } else {
    fold.diff.ests = cell.expr.ests.2 - cell.expr.ests.1   # Raw difference ests
  }

  if ( stdz ) {
    se = ((n1 * se1^2 + n2 * se2^2)/(n1 + n2))^(1/2)
    s0.r = apply(se, 1, quantile, 0.5, na.rm = TRUE)
    for (i in 1:numcell) {
      # fold.diff.ests[i ] = fold.diff.ests[i ]/(se[i] + s0.r[i])
      fold.diff.ests[i, ] = fold.diff.ests[i, ]/(se[i, ] + s0.r[i])
  } }

  if ( medCntr ) {
    for (i in 1:numcell) {
      # fold.diff.ests[i] = fold.diff.ests[i] - median(fold.diff.ests[i])
      fold.diff.ests[i, ] = fold.diff.ests[i, ] - median(fold.diff.ests[i, ])
  } }

  colnames( fold.diff.ests ) = colnames(het); rownames( fold.diff.ests ) = colnames( cell.props )
  # colnames( resids         ) = colnames(het); rownames( resids         ) = rownames( het        )

  return( list( fold.diff.ests, resids, deconv, se1, se2 ) )
} # end do.dual.decon
###############################################################################
###############################################################################
