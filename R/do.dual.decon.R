#' do.dual.decon
#'
#' Performs the core functionality of LRCDE.  This function does the actual linear regression deconvolution per group.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param het.sub The samples by genomic measures het.suberogeneous observations matrix
#' @param cell.props  Should be samples by cell types (rows by columns).  The relative cell proportions per sample.
#' @param groups Vector of 1's and 2's indicating group membership (controls and cases).
#' @param medCntr Boolean indicatinng whether to median center difference estimates.
#' @param stdz Boolean indicating whetehr to standardize difference estimates.
#' @param nonNeg Boolean indicating whet.subher to force cell type-specific expression estiamtes to be non-negative.
#' @return Returns a list containing group-wise difference estimates, regression residuals, lm objects per group, and standard errors of cell type-specific expression estimates.  Also returns a two item list containing results of Breusch-Pagan test.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details
#' Here is where the actual group-wise linear regressions are performed.
#' This function is designed to handle a single genomic site at a time; e.g.: a single vector containing het.suberogeneous observations for controls and cases.
#' This is the core functionality of linear regression cell type-specific differential expression detection.
#' @export
#' @examples
#' do.dual.decon( het.sub.obs, cell.props, groups, medCntr, stdz, nonNeg  )

do.dual.decon <- function( het.sub, cell.props, groups, medCntr, stdz, nonNeg  )
{

  ###########################################################################
  unique.groups = unique( groups )
  n = group.wise.sample.size( groups )
  n1 = n[1]
  n2 = n[2]
  ###########################################################################

  numgene  = ncol( het.sub )
  numcell  = ncol( cell.props )
  deconv <- list()
  ########################################
  # Group-wise linear regression:
  for ( curset in unique.groups ){ # print(curset)
    deconv[[ curset ]] = lm( het.sub[ groups==curset, ] ~ 0 + cell.props[ groups==curset, ] ) # through the origin
  }
  ########################################

  # Cell type specific expression estimates per group:
  cell.expr.ests.1 = deconv[[ 1 ]]$coefficients   # Cell-type specific expression estimates from group 1
  cell.expr.ests.2 = deconv[[ 2 ]]$coefficients   # Cell-type specific expression estimates from group 2

  # We use the residuals to calculate MSE per group per regression:
  resids1 = deconv[[1]]$residuals
  resids2 = deconv[[2]]$residuals
  if(dim(het.sub)[2] == 1){ # In case we're only looking at one site, make sure diff.ests is 2 dimensional
    resids1    = as.data.frame( resids1, ncol=1 ); colnames(resids1) = colnames(het.sub)
    resids2    = as.data.frame( resids2, ncol=1 ); colnames(resids2) = colnames(het.sub)
  }
  resids = rbind( resids1, resids2 )
  
  se1 = ls.diag(deconv[[ 1 ]])$std.err           # std err of Cell-type specific expression estimates group 1
  se2 = ls.diag(deconv[[ 2 ]])$std.err           # std err of Cell-type specific expression estimates group 2

  if (nonNeg) { # Force cell-type specific expression estimates to be non-negative.
    cell.expr.ests.1[cell.expr.ests.1 < 0] =  0.000000001
    cell.expr.ests.2[cell.expr.ests.2 < 0] =  0.000000001
  }

  # Get group-wise difference estimates:
  diff.ests  = array( dim = c( numcell, numgene ))
  diff.ests = cell.expr.ests.2 - cell.expr.ests.1   # Difference ests

  if ( stdz ) {
    se = ((n1 * se1^2 + n2 * se2^2)/(n1 + n2))^(1/2)
    s0.r = apply(se, 1, quantile, 0.5, na.rm = TRUE)
    for (i in 1:numcell) {
      # diff.ests[i ] = diff.ests[i ]/(se[i] + s0.r[i])
      diff.ests[i, ] = diff.ests[i, ]/(se[i, ] + s0.r[i])
  } }

  if ( medCntr ) {
    for (i in 1:numcell) {
      # diff.ests[i] = diff.ests[i] - median(diff.ests[i])
      diff.ests[i, ] = diff.ests[i, ] - median(diff.ests[i, ])
  } }

  if(dim(het.sub)[2] == 1){ # In case we're only looking at one site, make sure diff.ests is 2 dimensional
    diff.ests = as.data.frame( diff.ests, ncol=1 )
  }
  
  colnames( diff.ests ) = colnames(het.sub); rownames( diff.ests ) = colnames( cell.props )
  rownames( resids    ) = rownames(het.sub)

  return( list( diff.ests, resids, deconv, se1, se2 ) )
} # end do.dual.decon
###############################################################################
###############################################################################
