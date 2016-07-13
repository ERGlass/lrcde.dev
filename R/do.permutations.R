#' do.permutations
#'
#' Called from sim_driver_script.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param het  Should be samples by genomic site (rows by columns).  The samples by genomic measures heterogeneous observations matrix.
#' @param cell.props  Should be samples by cell types (rows by columns).  The relative cell proportions per sample.
#' @param groups  A vector of 1's and 2's indicating group membership per sample.  Should align with samples in het.sub and cell.props.
#' @param nperms Integer value indicating number of group label permutations to perform when simulating each site's parameters.
#' @param medCntr Boolean indicating whether to median center differential expression estimates.
#' @param stdz Boolean indicating whether to 'standardize' differential expression estimates.
#' @param nonNeg Boolean indicating whether to force cell type-specific expression estiamtes to be non-negative.
#' @return An 3-dimensional array of cell type-specific differential expression estimates for each permutation by each cell type by each genomic site.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' do.permutations( het.obs, cell.props, groups, n.perms=1000,  medCntr, stdz, nonNeg )
do.permutations <- function ( het, cell.props, groups, nperms=100, medCntr, stdz, nonNeg  ) 
{
  numgene = ncol( het )
  numcell = ncol( cell.props )
  div.diff.perm = array(dim = c( nperms, numcell, numgene ) )
  perm = list()
  unique.groups = unique( groups)
  n = vector()
  for( i in unique.groups) { n[i] = length(groups[groups == i]) }
  n1 = n[1]; n2=n[2];
  ########################################
  for (i in 1:nperms) {
    perm.groups = sample( groups )  # Permute group labels (whole point of doing permutations here).
    
    for(curset in unique.groups) {
      perm[[ curset ]] = lm( het[ perm.groups==curset, ] ~ 0 + cell.props[ perm.groups==curset, ] ) # through the origin
    }

    cell.expr.ests.1 = perm[[ 1 ]]$coefficients   # Cell-type specific expression estimates from group 1
    cell.expr.ests.2 = perm[[ 2 ]]$coefficients   # Cell-type specific expression estimates from group 2

    se1 = ls.diag(perm[[ 1 ]])$std.err           # std err of Cell-type specific expression estimates group 1
    se2 = ls.diag(perm[[ 2 ]])$std.err           # std err of Cell-type specific expression estimates group 2

    if (nonNeg) {
      cell.expr.ests.1[cell.expr.ests.1 < 0] = 0.000000001
      cell.expr.ests.2[cell.expr.ests.2 < 0] = 0.000000001
    }
    
      div.diff.perm[ i, , ] = cell.expr.ests.2 - cell.expr.ests.1   # or diff ests

    if ( stdz ) {
      se = ((n1 * se1^2 + n2 * se2^2)/(n1 + n2))^(1/2)
      s0.r = apply(se, 1, quantile, 0.5, na.rm = TRUE)
      for (i in 1:numcell) {
        div.diff.perm[i, , ] = div.diff.perm[i, , ]/(se[i, ] + s0.r[i])
      } }
    if ( medCntr ) {
      for (i in 1:numcell) {
        div.diff.perm[i, , ] = div.diff.perm[i, , ] - median(div.diff.perm[i, , ])
      } }

  } # end permutations
  
  ########################################
  return( div.diff.perm ) # Return all the permutation fold.x estimates
} #
