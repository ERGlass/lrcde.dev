#' custom.resids.synthetic - internal
#'
#' Used in 'lrcde.permutations_for_auc_variability' and 'using_lrcde'
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param mse2model.vec the vector of mean squared errors (MSEs) to model. Default - 0.05
#' @param groups  A vector of 1's and 2's indicating group membership per sample. 
#' Should align with the order of samples in het.sub and cell.props matrices. Required.
#' @param diff.2.model.vec vector of cell type-specific differential expressions to create.
#' Devault - c(0.1, 0.5, 1)
#' @param base.expr.vec vector of base (control group) cell type specific expressions to model.
#' Default - 10
#' @param adjuster a scaling factor used on the resulting standard deviation of residuals. Default - 1
#' @param n.cells number of cell types. Required
#' @return a matrix of synthetic normal residuals based on standard deviations of control 
#' and cases residuals from resids.gene.j. `length(group)` rows X `n.cells` columns
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' \dontrun{
#' resids <- custom.resids.synthetic(mse2model.vec = c(0.05), groups = c(rep(1, 15), rep(2, 15)), 
#'                                   diff.2.model.vec = c(0.1, 0.5, 1), base.expr.vec = 10, adjuster = 1, n.cells = 3)
#' }                                   

custom.resids.synthetic <- function(mse2model.vec = 0.05, groups, diff.2.model.vec = c(0.1, 0.5, 1), base.expr.vec = 10, adjuster = 1, n.cells )
{
  ###########################################################################
  n = group.wise.sample.size( groups );
  n.controls = n[1];
  n.cases    = n[2];
  n.samps = n.controls + n.cases
  ###########################################################################
  

  # Convert mses to sds:
  sd.2.model = (sqrt( mse2model.vec) / ( n.samps-n.cells ) ) * adjuster;

  basic.resids.2.model = matrix( rep(0, n.controls*length(mse2model.vec)), nrow = n.controls, ncol=length(mse2model.vec));

  reps = length(base.expr.vec) * length(diff.2.model.vec)
  # Add synthetic resids one gene at a time:
  basic.resids.2.model = matrix( rep(0, n.controls*length(mse2model.vec)), nrow = n.controls, ncol=length(mse2model.vec))
  for( j in 1:length(mse2model.vec)) {
    basic.resids.2.model[ , j ] = rnorm( n.controls, 0, sd.2.model[j])     # New column of resids
  }

  # Create zero matrix:
  num.genes.desired  = length(base.expr.vec) * length(diff.2.model.vec) * length(mse2model.vec);
  syn.resids.control = matrix( rep(0, n.controls * num.genes.desired ) , nrow = n.controls, ncol = num.genes.desired );
  add2 = 0;
  for( i in 1:length( mse2model.vec ) ) {
    for( j in 1:reps) {
      syn.resids.control[ , j + add2 ] = basic.resids.2.model[ , i ];
    }
    add2 = add2 + reps;
  }

  # Stack for identical resids on cases and controls:
  syn.resids = rbind( syn.resids.control, syn.resids.control )
  return( syn.resids )
}
