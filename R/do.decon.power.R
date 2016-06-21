#' do.decon.power
#'
#' This function is called by the lrcde function.  Not meant to be called directly by user.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param decon.list  List output by 'do.dual.decon' (or 'do.single.decon') function.  Containst lm fit object with regression results per group.
#' @return power.list  Results of power analysis.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details
#' This function is the main analysis part of the LRCDE package.  Here is where critical difference between control and case groups is estimated.
#' @export
#' @examples
#' do.decon.power( decon.list )
do.decon.power <- function(  decon.list, groups  )
{

  ###########################################################################
  unique.groups = unique( groups )
  n = group.wise.sample.size( groups )
  n.control = n[1]
  n.case    = n[2]
  ###########################################################################

  ###############################################################################
  diffs = decon.list[[1]]

  # Extract se for coefficients estimates (cell specific expressions) to compute CI on coefs:
  se1 = decon.list[[4]]
  se2 = decon.list[[5]]

  # Use se to get tail of 95% CI: (These are used to compute 'diff.critical' below.)
  if(alternative == "two.sided"){
    se1.tail = se1 * qt( 0.975, n.control - 1 , lower.tail=TRUE )
    se2.tail = se2 * qt( 0.975, n.case    - 1 , lower.tail=TRUE )
  } else {
    se1.tail = se1 * qt( 0.950, n.control - 1 , lower.tail=TRUE )
    se2.tail = se2 * qt( 0.950, n.case   - 1 , lower.tail=TRUE )
  }


  # Cell type specific expression estimates per group:
  deconv = decon.list[[3]]
  base.expr = deconv[[ 1 ]]$coefficients   # Cell-type specific expression estimates from group 1
  case.expr = deconv[[ 2 ]]$coefficients   # Cell-type specific expression estimates from group 2

  if(nonNeg) {
    base.expr[ base.expr < 0 ] = 0.000000001
    case.expr[ case.expr < 0 ] = 0.000000001
  }

  # case.dist =   ( base.expr + se1.tail ) - ( base.expr + abs( diffs ) )
  case.dist =   ( se1.tail ) - ( abs( diffs ) )
  power     = pt( case.dist, n.case - 1 , lower.tail=FALSE )

  ###########################################################################
  diff.critical = abs(diffs) - ( se1.tail + se2.tail )

  # If tail of 95% CI on base is greater than difference estimate,...
  # ...then greater than 0.5 power not possible:
  do.not.trust.these =  ( se1.tail > abs(diffs) )

  colnames(do.not.trust.these) = colnames(diffs)
  rownames(do.not.trust.these) = rownames(diffs)
  rownames(power) = rownames(diffs)
  rownames(base.expr) = rownames(diffs); rownames(case.expr)  = rownames(diffs)
  rownames(se1.tail)  = rownames(diffs); colnames(se1.tail)   = colnames(diffs)
  #############################################################################

  #############################################################################
  # Package and return power.list object:
  power.list = list( power, diffs,  diff.critical, base.expr, case.expr, do.not.trust.these, deconv , se1.tail )
  return( power.list )
  #############################################################################
}
