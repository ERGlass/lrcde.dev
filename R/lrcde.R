#' lrcde
#'
#' Cell type-specific differential expression detection given heterogeneous gene expression matrix,
#' cell type-specific cell proportions, and a vector of group assignment (e.g., case-control study).
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param het.sub the samples (rows) by genes (columns) heterogeneous gene expression matrix.
#' Should contain non-median-centered, non-standardized, positive values. log2-transformation recommended. Required
#' @param cell.props the cell proportion matrix, cell types (rows) by samples (columns). 
#' The proportion should be relative, e.g., sum up to ~1 per sample. Required
#' @param groups a vector of 1's and 2's indicating group assignment. 
#' Should correspond to the sample order in het.sub and cell.props. Required
#' @param output.file file name to output the comma-separated results. Default: LRCDE_power.analysis.csv.
#' @param VERBOSE boolean, whether to output progess to console. Default is TRUE.
#' @param method cpecifies which type of regression deconvolution to perform. Only "dual" method (deconvolution of each group separately) 
#' is implemented in the current version (default). Future methods may include "ridge" and/or "single".  
#' @param direction the type of test to perform. Should be one of "two.sided" (default), "up", or "down".
#' @return A list with three elements. The first contains a data.frame of analysis results.
#' The second contains a list of parameters supplied to the lrcde function. 
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details
#' Each row of the results data.frame contains cell type-specific differential expression analysis statistics.
#' Each column contains:
#' \describe{
#'  \item{site}{gene index, each gene is tested for differential expression in each cell type}
#'  \item{base}{cell type-specific gene expression estimate in _group 1_ (e.g., controls)}
#'  \item{case}{cell type-specific gene expression estimate in _group 2_ (e.g., cases)}
#'  \item{diff.est}{cell type-specific gene expression difference estimate. Larger differences are of interest, given sufficient power}
#'  \item{mse.control/case}{cell type-specific mean squared error gene expression estimatefor group 1/2, respectively}
#'  \item{cell}{cell type-specific index}
#'  \item{cell.sd}{cell type-specific standard deviation across samples}
#'  \item{kappa.1/2}{group-specific condition number for the cell proportion matrixes, for group 1/2, respectively}
#'  \item{t.crit}{t-statistics cutoff to reject the null hypothesis of no cell type-specific differential expression}
#'  \item{t.stat/p.val.t}{observed t-statistics/p-value, respectively, for the cell type-specific differential expression analysis. 
#'  Not corrected for multiple testing. ToDo: use '1 - p' for negative diff.est}
#'  \item{se.1/se.2}{standard error of the cell type-specific gene expression estimates for group 1/2, respectively}
#'  \item{t.power}{power of the cell type-specific differential analysis. Non-siginficant power defaults to p.val.t, power at p.val.t < 0.05 indicates potentially significant cell type-specific differences. Larger power (> 0.8) corresponds to more significant differences.}
#' }
#' @export
#' @examples
#' \dontrun{
#' library("lrcde")
#' library(CellMix)
#' library(GEOquery)
#' 
#' gedDataInfo() # Data sets in CellMix to select from
#' mix <- ExpressionMix("GSE19830") # Example dataset with gene expression form mixture of tissues
#' p.sub <- subset(pData(mix), select = c(Liver, Brain, Lung, characteristics_ch1))
#' 
#' het <- t(exprs(mix))  # Heterogeneous expression matrix
#' cell.props <- t(coef(mix)) # Cell proportions matrix
#' 
#' apply(cell.props, 2, sd) # SD of cell proportions should not be 0
#' # Cell proportion standard deviation > 0.06 across samples should
#' # produce reasonable differential detection power.  
#' apply(cell.props, 1, sum)  # These sum perfectly to 1.
#' 
#' # Create (artificially) two groups, Brain and Lung vs. Liver
#' p.sub.1 <- p.sub[(p.sub$Brain + p.sub$Lung) > (p.sub$Liver), ]
#' # Heterogeneous matrixes for each group
#' het.1 <- het[rownames(het) %in% rownames(p.sub.1), ]
#' het.2 <- het[!(rownames(het) %in% rownames(p.sub.1)), ]
#' het2use <- rbind(het.1, het.2) # Combined matrix
#' # Cell proportions for each group
#' cells.1 <- cell.props[rownames(cell.props) %in% rownames(p.sub.1), ]
#' cells.2 <- cell.props[!(rownames(cell.props) %in% rownames(p.sub.1)), ]
#' cell.props <- rbind(cells.1, cells.2) # Combined cell proportions
#' # Group vector
#' groups <- c(rep(1, dim(het.1)[1]), rep(2, dim(het.2)[1]))
#' 
#' # Now apply the lrcde function in order to deconvolve cell
#' # type-specific expression and perform power analysis on differential
#' # expression detection:
#' 
#' power.analysis.df <- lrcde(het2use, 
#'                            cell.props, 
#'                            groups, 
#'                            output.file = "LRCDE_power_analysis_results.csv", 
#'                            VERBOSE = TRUE, 
#'                            method = "dual", 
#'                            direction = "two.sided")
#' }

lrcde <- function(het.sub,
                  cell.props,
                  groups,
                  output.file = "LRCDE_power_analysis.csv",
                  VERBOSE     = TRUE,
                  method      = "dual",
                  direction   = "two.sided") {
# Default parameters kept for historical reason
  medCntr     = FALSE
  stdz        = FALSE
  nonNeg      = TRUE
  
  options(stringsAsFactors = FALSE)
  
  ###########################################################################
    # TEST cell.proportions here
  ###########################################################################
  
  ###########################################################################
  # Alpha significance level for t-tests:
  alpha = 0.05
  ###########################################################################
  
  ###########################################################################
  unique.groups <- unique(groups)
  n             <- group.wise.sample.size(groups)
  n.control     <- n[1]
  n.case        <- n[2]
  n.cells       <- dim(cell.props)[2]
  ###############################################################################

  ###############################################################################
  ###########################################################################
  # Do group-wise regressions (two regressions per genomic site):
  if (method == "dual") {
    decon.list <- do.dual.decon( het.sub, cell.props, groups, medCntr, stdz, nonNeg   )
  }
  
  # Returned by do.dual.decon function:
  # decon.list =  list( fold.diff.ests, resids, deconv, se1, se2  ) )
  ###############################################################################

  ###############################################################################
    # NOTE: For "dual" regression, t-sample t-test is performed "by hand" below:
  # Standard two-sample t-test (Welch's test), where assumption is se1 not equal se2:
  # Balanced or unbalanced samples, and variances not necessarily equal between groups:
  se1 = decon.list[[4]]
  se2 = decon.list[[5]]  # Group-wise standard errors
  diff.ests = decon.list[[1]]                   # Differential expression estimate
  # Welche's "pooled" standard error:
  se.welches =  sqrt(((se1^2)/n.control) + ((se2^2)/n.case))
  # T-statistic:
  t.stats = ( diff.ests )  /  se.welches
  # d.f.equal = ( n.control + n.case - 2 )            # Pooled degrees of freedom
  # NOTE: Welches degrees of freedom is reduced from pooled degrees of freedom when standard errors are unequal:
  d.f.numerator = (((se1^2)/n.control) + ((se2^2)/n.case))^2
  d.f.denom    = ((se1^2)/n.control)^2 / (n.control-1)  +  ((se2^2)/n.case)^2 / (n.case-1)
  d.f.welches  = d.f.numerator / d.f.denom
  
  # Alpha = 0.05 (total):
  if(direction=="two.sided"){
    alpha.half = alpha/2
  } else {
    alpha.half = alpha
  }
  
  t.crit = qt( 1-alpha.half, d.f.welches )
  
  p.val.t   =  pt( t.stats, d.f.welches, lower.tail=FALSE )
  # OUTPUT: t.stats, p.val.t
  t.stats.out = as.vector( t(t.stats))
  p.val.t.out = as.vector( t(p.val.t))
  se1.out     = as.vector( t( se1 ) )
  se2.out     = as.vector( t( se2 ) )
  se.welches.out = as.vector( t(se.welches))
  t.crit.out  = as.vector( t(t.crit) )
  ###############################################################################

  ###############################################################################
  # Power from t-statistic:
  if(direction=="two.sided"){
    power.test = 1 - pt( t.crit, d.f.welches, t.stats ) + pt( (-1 * t.crit), d.f.welches, t.stats )
  } else { # one sided
    power.test = 1 - pt( t.crit, d.f.welches, abs(t.stats) )
  }
  power.t.out = as.vector(t(power.test))
  ###############################################################################

  diffs.t  = as.data.frame( t( diff.ests ) )

  ###############################################################################
  # Cell proportions (cell.props) statistics:
  cell.props.control  <- cell.props[ groups == 1, ]
  cell.props.case     <- cell.props[ groups == 2, ]
  cell.SDs            <- apply( cell.props, 2, sd )
  # Condition numbers for the cell proportion matrixes
  kappa.control       <- kappa( t( cell.props.control) %*% cell.props.control , exact = TRUE )
  kappa.case          <- kappa( t( cell.props.case   ) %*% cell.props.case    , exact = TRUE )
  kappa.all           <- kappa( t( cell.props        ) %*% cell.props         , exact = TRUE )
  ###############################################################################

  # ###########################################################################
  # MSE of obs:
  resids = decon.list[[2]]
  resids.control   = resids[ groups == 1, , drop=FALSE ]
  resids.case      = resids[ groups == 2, , drop=FALSE ]
  resids.control.2 = resids.control ^ 2
  resids.case.2    = resids.case    ^ 2
  sse.control      = apply( resids.control.2 , 2, sum )
  sse.case         = apply( resids.case.2    , 2, sum )
  mse.control      = sse.control / ( n.control - n.cells )
  mse.case         = sse.case    / ( n.case    - n.cells )

  mse.control.frame = data.frame( as.character(colnames(het.sub)), mse.control, stringsAsFactors = F )
  mse.case.frame    = data.frame( as.character(colnames(het.sub)), mse.case   , stringsAsFactors = F )
  colnames(mse.control.frame) = c( "site", "mse.control" )
  colnames(mse.case.frame)    = c( "site", "mse.case"    )
  ###########################################################################

  #############################################################################
  # Assemble data.frame for output:

  deconv = decon.list[[3]]
  cell.expr.ests.1 = deconv[[ 1 ]]$coefficients   # Cell-type specific expression estimates from group 1
  cell.expr.ests.2 = deconv[[ 2 ]]$coefficients   # Cell-type specific expression estimates from group 2
  
  base.t   = as.data.frame(t( cell.expr.ests.1 ))
  case.t   = as.data.frame(t( cell.expr.ests.2 ))

  colnames(base.t) = colnames(diffs.t)
  colnames(case.t) = colnames(diffs.t)
  
  diffs.t$site  = rownames( diffs.t  )
  base.t$site   = rownames( base.t   )
  case.t$site   = rownames( case.t   )

  cell.names = colnames(cell.props)
  numcells   = length(cell.names)

  total.frame = data.frame()
  for( p in 1:numcells){
    tmp.frame = data.frame()
    cell.name = cell.names[p]


    text2parse = paste0("diffs.tmp = subset(diffs.t, select=c(", cell.name ,", site ))")
    eval(parse(text=text2parse));
    colnames(diffs.tmp) = c("diff.est", "site")


    text2parse = paste0("base.tmp = subset(base.t, select=c( site , ", cell.name ,"  ))")
    eval(parse(text=text2parse));
    colnames(base.tmp) = c( "site", "base")
    
    text2parse = paste0("case.tmp = subset(case.t, select=c(", cell.name ,", site ))")
    eval(parse(text=text2parse));
    colnames(case.tmp) = c("case", "site")
    
    
    tmp.frame = dplyr::left_join( base.tmp   , case.tmp          , by="site" )
    tmp.frame = dplyr::left_join( tmp.frame  , diffs.tmp         , by="site" )
    tmp.frame = dplyr::left_join( tmp.frame  , mse.control.frame , by="site" )
    tmp.frame = dplyr::left_join( tmp.frame  , mse.case.frame    , by="site" )
    tmp.frame$cell    = cell.name

    cat("cell name: ", cell.name, "\n")

    tmp.frame$cell.sd = cell.SDs[p]
    tmp.frame$kappa.1 = kappa.control
    tmp.frame$kappa.2 = kappa.case

    total.frame = rbind(total.frame, tmp.frame)
  }
  #############################################################################

  #############################################################################
  # Output t-statistic results:
  for(j in 1:length(p.val.t.out)){
    if( p.val.t.out[j] > alpha.half ){
      power.t.out[j] = p.val.t.out[j]
    }
  }
  
  old.names = colnames(total.frame)
  total.frame = cbind(total.frame, t.crit.out, t.stats.out, p.val.t.out, se1.out, se2.out, se.welches.out, power.t.out )
  colnames(total.frame) = c(old.names, "t.crit", "t.stat", "p.val.t" , "se1", "se2", "se.p", "t.power")
  #############################################################################

  #############################################################################
  # Write total.frame to .CSV file
  unlink(output.file) # Delete, if file exists
  write.csv(total.frame, file = output.file, row.names = FALSE)
  #############################################################################

  #############################################################################
  # Package the return list object:
  args.used = list(output.file,
                   nonNeg,
                   method,
                   direction)

  return.list = list( total.frame, args.used, decon.list )
  return( return.list )
  #############################################################################
} # End LRCDE
