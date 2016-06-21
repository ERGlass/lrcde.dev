#' sim.cell.expr
#'
#' Called from sim_driver_script.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param num.cells Number of cell types to simulate.
#' @param num.genes Number of genomic sites to simulate.
#' @param base.expr Control group expression level to simulate at the cell type-specific level.
#' @param diff.expr Absolute difference between control group and cases group to simulate at the cell type-specific level.
#' @param prop.2.fold Proportion of simulate genomic sites to "fold" or alter between controls and cases.
#' @param cell.p Number of cell type to alter between controls and cases for simulation purposes.
#' @return List containing group-wise matrix of simulated cell type-specific expressions (rbinded as controls and cases.) and "truth" vector of 0's and 1's indicating which genomic sites have been altered for cell type p.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' sim.cell.expr(num.cells.2.simulate, num.genes.2.simulate, base.expression, difference.2.simulate, proportion.of.sim.genes.2.fold, cell.2.modify )

sim.cell.expr <- function(  num.cells
                            , num.genes   = 5000
                            , base.expr   = 5
                            , diff.expr   = 1
                            , prop.2.fold = .5
                            , cell.p )
{
  
  # Create two matrices: Controls and Cases: n.cells by n.genes:
  cells.control = matrix( rep( base.expr, num.genes * num.cells )
                          , nrow = num.cells, ncol=num.genes )
  cells.cases   = cells.control
  
  num2fold      = prop.2.fold * num.genes
  genes2fold    = 1:num2fold
  
  length.folded = length( genes2fold )
  length.non    = num.genes - length.folded
  truth         = c( rep(1, length.folded), rep(0, length.non ) )
  
  # Apply differential expression to cell p in cases:
  cells.cases[ cell.p , genes2fold ] =  cells.cases[ cell.p , genes2fold ] + diff.expr
  
  # combine controls and cases:
  cell.expr = rbind( cells.control, cells.cases )
  
  return( list( cell.expr, truth ) )
}
###############################################################################
#  cells.expr = sim.cell.expr( 5, 20, 200, 200, .5 ) # testing 1, 2, 3...
###############################################################################
