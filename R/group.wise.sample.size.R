#' group.wise.sample.size
#'
#' Utility function to keep from having to pass group-wise sample sizes between functions.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/ERGmisc}
#' @param groups  A vector of 1's and 2's indicating group membership per sample.  Should align with samples in het.sub and cell.props.
#' @return n  A vector containing control and cases sample sizes.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @export
#' @examples
#' group.wise.sample.size( groups )

group.wise.sample.size <- function( groups ){
  unique.groups = unique( groups )
  n = vector()
  for( i in unique.groups) { n[i] = length( groups[ groups == i ] ) }
  # n.controls = n[1]
  # n.cases    = n[2]
  return(n)
}