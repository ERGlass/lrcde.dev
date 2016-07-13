#' get.permuted.p.values
#'
#' Called from sim_driver_script
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param diff.perms.func  Differential expression estimates from permutations (a 3D array: permutations by cell types by genomic sites).
#' @param diffs.func  Differential expression estimates from group label non-permuted observations.
#' @param alt.is  Defaults to 'greater' to test for up-regulation of a genomic site.  All differences are modled by absolute value, so default is 'greater'.
#' @param cell.p  Integer of cell type being tested for differences between truth (diffs.func) and permuted differences (diff.perms.func).
#' @return sig.perms.func  Vector of permutation p.values for each genomic site.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' get.permuted.p.values( 3.D.array.of.permuted.difference.estimates, actual.difference.estiamtes, alt.is='greater', cell.p=1 )

get.permuted.p.values <- function( diff.perms.func, diffs.func, alt.is, cell.p )
{
  # sig.perms.func <- array(dim = c(numcell, numgene)); # Create location to store empirical p-values per gene
  numgene   = dim( diffs.func )[2]
  num.perms = dim( diff.perms.func )[1]
  
  sig.perms.func <- rep(0, length=numgene); # Create location to store empirical p-values per gene
  #  for(curcell in 1:numcell)  {
  
  if(alt.is == "two.sided"){
    for (curgene in 1:numgene) {
      # sig.perms.func[ curgene ] <- (sum( abs( diff.perms.func[ , cell.p, curgene ]) > abs( diffs.func[ cell.p, curgene] ))) / num.perms
      sig.perms.func[ curgene ] <- ( 1 + sum( abs( diff.perms.func[ , cell.p, curgene ] ) >= abs( diffs.func[ cell.p, curgene] ))) / ( 1 + num.perms )
    } #
  } # two sided
  
  if(alt.is == "greater"){
    for (curgene in 1:numgene) {
      sig.perms.func[ curgene ] <- ( 1 +  sum( ( diff.perms.func[ , cell.p, curgene ]) >= ( diffs.func[ cell.p, curgene] )) ) / ( 1 + num.perms )
    } #
  } # greater
  
  if(alt.is == "less"){
    for (curgene in 1:numgene) {
      sig.perms.func[ curgene ] <- ( 1 + sum( ( diff.perms.func[ , cell.p, curgene ]) <= ( diffs.func[ cell.p, curgene] )) ) / ( 1 + num.perms )
    } #
  } # less
  
  return( sig.perms.func )
} # get.permuted.p.values
###############################################################################
###############################################################################
