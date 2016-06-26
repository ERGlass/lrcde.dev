#' custom.resids.synthetic
#'
#' Called from user lrcde wrapper script (see example 'lrcde.wrapper.sims.R').
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param resids.gene.j The residuals for both study groups for the genomic site j of interest.
#' @param groups  A vector of 1's and 2's indicating group membership per sample.  Should align with samples in het.sub and cell.props.
#' @param num.genes.desires Number of genomic sites to simulate based on genomic site j.
#' @return syn.resids  Synthetic normal residuals based on standard deviations of control and cases residuals from resids.gene.j.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' custom.resids.synthetic(het.obs, cell.props, groups, n.perms=1000, file.ext="my.results", cell.vector.to.analyze, proj.root="Project.directory", output.dir="output.dir" )
custom.resids.synthetic <- function( mse2model.vec, groups, diff.2.model.vec, base.expr.vec, adjuster=1 )
{
  ###########################################################################
  n = group.wise.sample.size( groups );
  n.controls = n[1];
  n.cases    = n[2];
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
