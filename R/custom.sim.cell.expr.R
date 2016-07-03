#' custom.sim.cell.expr
#'
#' Called from a user wrapper script (see example 'lrcde.wrapper.sims.R').
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param num.cells Number of cell types to simulate.
#' @param base.expr.vec Control group expression level to simulate at the cell type-specific level.
#' @param diff.expr.vec Absolute difference between control group and cases group to simulate at the cell type-specific level.
#' @param cell.p Number of cell type to alter between controls and cases for simulation purposes.
#' @return cell.expr Cell expressions matrix.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' custom.sim.cell.expr( num.cells.2.simulate, base.expression.vector, difference.2.simulate.vector, cell.2.modify )

custom.sim.cell.expr <- function(  num.cells=5
                                  , base.expr.vec
                                  , diff.2.model.vec
                                  , cell.p
                                  , len.mse.vec)
{
  # Create two matrices: Controls and Cases: n.cells by n.genes:
  base.reps =  rep( base.expr.vec, times = ( length( diff.2.model.vec) * len.mse.vec ) )
  cells.control = matrix( rep( base.reps, num.cells ), nrow = num.cells, byrow=T )
  cells.cases   = cells.control

  # Apply differential expression to cell.p in cases:
  diff.reps = length( base.expr.vec)
  diff.rep.vec = rep( diff.2.model.vec, each = diff.reps, times = len.mse.vec )
  
  for(p in cell.p) {# print(p)
    cells.cases[ p ,  ] =  cells.control[ p ,  ] + diff.rep.vec
  }

  # combine controls and cases:
  cell.expr = rbind( cells.control, cells.cases )

  return( cell.expr )
}
