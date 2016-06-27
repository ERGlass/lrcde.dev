#' group.wise.sample.size
#'
#' Utility function to avoid passing group-wise sample sizes between functions.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param groups A vector of identifiers indicating group membership per sample.
#' Should have the same length as the number of rows in het.sub and cell.props.
#' @return n  A vector containing group-wise sample sizes.
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @examples
#' \dontrun{
#' group.wise.sample.size(group = rep(1:3, each = 5)) # Three groups, 5 samples each
#' }

group.wise.sample.size <- function(groups) {
  unique.groups <- unique(groups)
  n <- vector(mode = "numeric", length = length(unique.groups))
  for (i in unique.groups) {
    n[i] <- length(groups[groups == i])
  }
  return(n)
}