#' lrcde.permutations
#'
#' Call this function to run entire functionality of LRCDE.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param het.sub  Should be samples by genomic site (rows by columns).  The samples by genomic measures heterogeneous observations matrix.
#' @param cell.props  Should be samples by cell types (rows by columns).  The relative cell proportions per sample.
#' @param groups  A vector of 1's and 2's indicating group membership per sample.  Should align with samples in het.sub and cell.props.
#' @param medCntr Boolean indicating whether to mean center differential expression estimates.
#' @param stdz Boolean indicatin whether to scale differential expression estimates with their pooled adjusted standard deviation
#' @param nonNeg Boolean indicating whether to force cell type-specific estimates to be non-negative (TRUE) or not (FALSE).
#' @param method Only "dual" is implemented in this version. This should be one of "dual" (csSAM method), or "ridge".  Default is "dual".  Specifies which type of regression deconvolution to perform.
#' @param direction Should be one of "two.sided", "up", or "down".  Which direction to test for cell type-specific expression changes.
#' @param n.perms Number of permutations of group labels to perform
#' @param cell.p The simulation differentiated target cell.
#' @param truth The truth indicator vector of which simulated "genes" in the target cell are differentiated between cases and controls.
#' @return List containing data.frame (total.frame) of analysis results and a list of parameter values supplied to lrcde function (arg.used).
#' @export

lrcde.permutations <- function(   het.sub
                                , cell.props
                                , groups
                                , medCntr     = FALSE
                                , stdz        = FALSE
                                , nonNeg      = TRUE
                                , method      = "dual"
                                , direction   = "two.sided"
                                , n.perms    = 100
                                , cell.p     = 1
                                , truth
                                ) {
  
  require( pROC  )  # Hands us area under curve as a scalar.
  
  ###########################################################################
  unique.groups <- unique(groups)
  n             <- group.wise.sample.size(groups)
  n.control     <- n[1]
  n.case        <- n[2]
  n.cells       <- dim(cell.props)[2]
  ###########################################################################
  
  ###########################################################################
  # Do group-wise regressions (two regressions per genomic site):
  if (method == "dual") {
    decon.list <- do.dual.decon( het.sub, cell.props, groups
                                 , medCntr=medCntr, stdz=stdz, nonNeg=nonNeg)
  }
  ###########################################################################
  
  ###########################################################################
  # Do permutations:
  fold.diff.perms = do.permutations( het.sub, cell.props, groups=groups
                                     , nperms=n.perms
                                     , medCntr=medCntr, stdz=stdz
                                     , nonNeg=nonNeg   )
  ###########################################################################
  
  ###########################################################################
  # Get permuted 'exact' p-values:
  # Hard coded to 'greater' since we take ABS of estimated differences
  p.perms <- get.permuted.p.values( fold.diff.perms, decon.list[[1]]
                                    , alt.is=direction, cell.p )
  ###########################################################################

  ###########################################################################
  # Uses the 'pROC' package for 'roc' function: >require(pROC)
  roc1 <- pROC::roc( truth, p.perms, ci=TRUE )
  auroc <- round( roc1$auc, 3 )
  ###########################################################################
  
  ###############################################################################
  # Cell proportions (cell.props) statistics:
  cell.props.control  <- cell.props[groups == 1, ]
  cell.props.case     <- cell.props[groups == 2, ]
  cell.SDs            <- apply(cell.props, 2, sd)
  
  # Condition numbers for the cell proportion matrixes
  # (cases and controls are identical):
  kappa.control       <- kappa(t(cell.props.control) %*% cell.props.control
                               , exact = TRUE)
  ###############################################################################

  return.list = list( p.perms , auroc, kappa.control  )
  return( return.list )
  #############################################################################

} # End LRCDE
