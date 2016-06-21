#' het.from.synthetic
#'
#' Called from sim_driver_script.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param cell.props  Should be samples by cell types (rows by columns).  The relative cell proportions per sample.
#' @param cell.expr Matrix of simulated cell type-specific expressions (returned from sim.cell.expr()).
#' @param resids The simulated residuals returned from get.resids.synthetic().
#' @param groups Vector of 1's and 2's indicating group membership (controls and cases).
#' @return het.obs The simulated samples by genomic sites heterogeneous observations matrix with simulated residuals applied.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' het.from.synthetic( cell.props, simulated.cell.expression, simulated.residuals, groups )

het.from.synthetic <- function( cell.props, cell.expr, resids, groups )
{

  ###########################################################################
  n = group.wise.sample.size( groups )
  n.controls = n[1]
  n.cases    = n[2]
  ###########################################################################
  
  n.cells  = dim( cell.props )[2]
  n.plus.1 = n.cells + 1
  P        = 2 * n.cells
  
  # Divvy up cell.expr into controls and cases:
  cell.expr.control  = cell.expr[ 1:n.cells  , ]
  cell.expr.cases    = cell.expr[ n.plus.1:P , ]

  # Divvy up cell.props into controls and cases:
  cell.props.control  = cell.props[ 1:n.controls  , ]
  cntl.plus.1 = n.controls + 1; N = n.controls + n.cases
  cell.props.cases    = cell.props[ cntl.plus.1:N , ]

  # Make "perfect fitted values" (these "fall on the regression line"):
  het.fit.control  = cell.props.control %*% cell.expr.control
  het.fit.cases    = cell.props.cases   %*% cell.expr.cases
  het.fitted       = rbind( het.fit.control, het.fit.cases )

  # Add resids:
  het.obs = het.fitted + resids

  # Return "fitted" values:
  return( het.obs )
}

###############################################################################
###############################################################################
