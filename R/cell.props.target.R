
#' cell.props.target
#'
#' Creates synthetic cell proportions which targets both a specific cell proportion standard deviation for cell 1 (the cell used for simulation purposes) and condition number over entire cell proportion matrix.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param n.samps Sample size for cell proportions matrix to create.
#' @param target.sd Standard deviation to create across samples for cell of interest.
#' @param target.cd Condition number for whole cell proportions matrix.
#' @return cell.props An n.samps by cell types matrix of simulated relative cell proportions.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details 
#' NOTE: DO not attempt to specify a condition number less than 15.  Function will crash.
#' @export
#' @examples
#' cell.props.target(  n.samps, target.sd, target.cd )
#' 
cell.props.target <- function( n.cells, n.samps, target.sd, target.cd  ) {
  #######################################################################
  # Create 1st cell with desired SD:
  start.prop    = 0.01
  tol = 0.001
  end.prop = 0.001
  #########################    
  while(TRUE){
    cell.1 = seq( start.prop, end.prop, length.out=n.samps )
    sd.1 = sd(cell.1)
    diff = target.sd - sd.1
    #     cat("diff: ", diff, ", tol: ", tol)
    if( abs(diff) <= tol ) break;
    if( diff > 0 ) {
      end.prop = end.prop + 0.0001    # increase target sd
    } else {
      end.prop = end.prop - 0.0011    # decrease target sd
    } # if diff ...
  } # while
  #######################################################################
  
  # Make more cells with noise:  
  tot.rem = 1 - cell.1
  n.cells.minus.1 = ( n.cells - 1 )
  # per.cell = tot.rem / n.cells.minus.1
  
  # New cells:
  scale.tune = 0.2
  cd.diff = 100000 # Start with arbitrarily huge difference
  cd.tol = 0.01 * target.cd
  
  while( TRUE ){ # Target a condtion number for whole matrix:
    
    if(cd.diff > 0){ # If tmp kappa is too small, decrease SD of remaining cells:
      scale.tune = scale.tune * .9999
    } else {          # Else, if tmp kappa is too big, increase SD of remaining cells:
      scale.tune = scale.tune * 1.001
    }
    
    more.cells = matrix(, ncol=n.cells.minus.1, nrow=n.samps)
    more.cells[,1] = cell.1
    for(r in 2:(n.cells.minus.1)) {
      more.cells[,r] = runif( n.samps , 0 , tot.rem * scale.tune )
      tot.rem = 1 - apply( more.cells[,1:r], 1, sum )
    }
    
    # Create LAST cell:
    sum.mor = apply( more.cells, 1, sum )
    last.cell = 1-sum.mor             # This forces VIF to Inf (row sums to 1)
    cell.props.tmp = cbind( more.cells, last.cell )
    
    # Get condition number and diff from target condtion number:
    cc.cc = t( cell.props.tmp ) %*% cell.props.tmp
    cd.tmp = kappa( cc.cc, exact=TRUE )
    cd.diff = target.cd - cd.tmp
    
    # cat("tmp kappa: ", cd.tmp, ", diff: ", cd.diff, "\n")  
    if ( abs(cd.diff) < cd.tol ) break;
  } # While
  colnames( cell.props.tmp ) = NULL  
#######################################################################
#######################################################################
return( cell.props.tmp )
} # End cell.props.target

