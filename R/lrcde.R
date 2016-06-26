#' lrcde
#'
#' Call this function to run entire functionality of LRCDE.
#' @author Edmund R Glass, \email{Edmund.Glass@@gmail.com}, Mikhail G Dozmorov, \email{Mikhail.Dozmorov@@vcuhealth.org}
#' @references \url{https://github.com/ERGlass/lrcde.dev}
#' @param het.sub  Should be samples by genomic site (rows by columns).  The samples by genomic measures heterogeneous observations matrix.
#' @param cell.props  Should be samples by cell types (rows by columns).  The relative cell proportions per sample.
#' @param groups  A vector of 1's and 2's indicating group membership per sample.  Should align with samples in het.sub and cell.props.
#' @param output.file File or path and file to output.  If indicated output directory (if path indicated) does not exist, a warning is issued and program execution halts.  Default behavior is to write output to LRCDE_power.analysis.csv in the current working directory.
#' @param FEEDBACK Boolean indicating whether to output progess indication to console.  Default is TRUE.
#' @param medCntr Boolean indicating whether to mean center differential expression estimates.
#' @param stdz Boolean indicatin whether to scale differential expression estimates with their pooled adjusted standard deviation
#' @param nonNeg Boolean indicating whether to force cell type-specific estimates to be non-negative (TRUE) or not (FALSE).
#' @param method Only "dual" is implemented in this version. This should be one of "single", "dual" (csSAM method), or "ridge".  Default is "dual".  Specifies which type of regression deconvolution to perform.
#' @param direction Should be one of "two.sided", "up", or "down".  Which direction to test for cell type-specific expression changes.
#' @return List containing data.frame (total.frame) of analysis results and a list of parameter values supplied to lrcde function (arg.used).
#' @keywords Deconvolution cell type-specific differential expression detection power analysis
#' @details
#' The lrcde function is meant to be called directly by user.  It is the entry point for the cell type-specific differential expression power analysis.
#' This is mainly a wrapper script for the deconvolution step and subsequent call to power analysis function (do.decon.power).
#' Default behavior is to write output to LRCDE_power.analysis.csv in the current working directory.
#' @export
#' @examples
#' # Load lrcde library:
#' library( "lrcde" )
#'
#' # User chosen working directory - Where your output .csv file will be found:
#' setwd(getwd()); getwd()
#' 
#' # Load data to work with:
#' # Select an ExpressionSet object from the CellMix package to work with:
#' library(CellMix)
#' library(GEOquery)
#' 
#' # Data sets in CellMix to select from:
#' gedDataInfo()
#' mix <- ExpressionMix("GSE19830")
#' p.sub = subset(pData(mix), select=c( Liver, Brain, Lung,  characteristics_ch1))
#' head( p.sub )
#' unique( p.sub$characteristics_ch1 )
#' het <- t(exprs(mix))            # Heterogeneous expression matrix
#' cell.props <- t(coef(mix))      # cell proportions matrix
#' dim(het)
#' dim(cell.props)
#' apply(cell.props, 2, sd)
#' # Cell proportion standard deviation > 0.06 across samples should produce reasonable differential detection power.
#' # Brain     Liver      Lung 
#' # 0.2613150 0.2859059 0.2837341
#' apply(cell.props, 1, sum)            # These sum perfectly to 1.  NOT suitable for single-step deconvolution.
#' 
#' # Create (artificially) two groups
#' p.sub.1 = p.sub[ (p.sub$Brain + p.sub$Lung) >  ( p.sub$Liver )     ,]
#' dim(p.sub.1)
#' g.1.names = rownames(p.sub.1)
#' het.1 = het[ rownames(het) %in% g.1.names, ]
#' dim(het.1)
#' het.2 = het[ !(rownames(het) %in% g.1.names), ]
#' dim(het.2)
#' cells.1 = cell.props[ rownames(cell.props)%in%g.1.names,  ]
#' cells.2 = cell.props[ !(rownames(cell.props)%in%g.1.names),  ]
#' all.equal(rownames(het.1), rownames(cells.1))
#' all.equal(rownames(het.2), rownames(cells.2))
#' het2use = rbind(het.1, het.2)
#' cell.props = rbind(cells.1, cells.2)
#' groups = c(rep(1, dim(het.1)[1]), rep(2, dim(het.2)[1]) )
#' table(groups)
#' 
#' # Now apply the lrcde function in order to deconvolve cell type-specific expression
#' # and perform power analysis on differential expression detection:
#' 
#' power.analysis.df = lrcde( het2use, cell.props, groups
#'                            , output.file = "LRCDE_power_analysis_results"
#'                            , FEEDBACK = TRUE, medCntr = FALSE, stdz = FALSE, nonNeg = TRUE
#'                            , method = "dual", direction = "two.sided")
#' 
#' 
#' 
lrcde <- function(  het.sub
                  , cell.props
                  , groups
                  , output.file="LRCDE_power_analysis"
                  , FEEDBACK  = TRUE
                  , medCntr   = FALSE
                  , stdz      = FALSE
                  , nonNeg    = TRUE
                  , method    = "dual"
                  , direction = "two.sided"
                  )
  {

  ###########################################################################
  unique.groups = unique( groups )
  n = group.wise.sample.size( groups )
  n.control = n[1]
  n.case     = n[2]
  n.cells = dim(cell.props)[2]
  ###########################################################################

  ###############################################################################
  # Initial checks and warnings:

    # Check for existence of indicated output file and directory and break nicely if file exists or directory does not exist.
  #   test.4.file = paste0( output.file, ".csv" )
  #   if(file.exists(test.4.file)){ # Do not clobber an existing file.
  #     cat("Warning: indicated output file: '", output.file, "' exists.  You may want to rename indicated output file.")
  #   } else

  {
  # OK to write file:
  # if( !dir.exists( dirname(output.file) )) { # dir.exists function not working in Windows 7
  #   cat(" Output directory indicated: '", dirname(output.file), "' does not exist.  Please create desired output directory first, then re-run lrcde.\n")
  # } else  { # output.dir exists... carry on:
  ###############################################################################

  ###############################################################################
    # Do the actual deconvolution step:

    # Do regressions (one regression per genomic site):
    # if( method == "single" ) {
    #   decon.list   = do.single.decon(  het.sub
    #                                   , cell.props
    #                                   , groups
    #                                   , medCntr, stdz, nonNeg
    #                                   )
    # } # single

    # Do group-wise regressions (two regressions per genomic site):
    if( method == "dual" ) {
      decon.list   = do.dual.decon(   het.sub
                                    , cell.props
                                    , groups
                                    , medCntr, stdz, nonNeg
                                    )
    } # dual

    # # NOT IMPEMENTED:
    # # Do group-wise regressions (two regressions per genomic site):
    # if( method == "ridge" ) {
    #   decon.list =   do.ridge.decon(  het.sub
    #                                  , cell.props
    #                                  , groups
    #                                  , medCntr, stdz, nonNeg
    #                                 )
    # } # ridge

    ###############################################################################

  # Returned by any of the above do.*.decon functions:
  # decon.list =  list( fold.diff.ests, resids, deconv, se1, se2  ) )

  ###############################################################################
  # Cell proportions (cell.props) statistics:
  cell.props.control = cell.props[ groups == 1, ]
  cell.props.case    = cell.props[ groups == 2, ]
  cell.SDs = apply( cell.props, 2, sd )
  kappa.control = kappa( t( cell.props.control ) %*% cell.props.control, exact=TRUE )
  kappa.case    = kappa( t( cell.props.case    ) %*% cell.props.case   , exact=TRUE )
  kappa.all     = kappa( t( cell.props    ) %*% cell.props   , exact=TRUE )
  ###############################################################################

  # ###########################################################################
  # # Do Shapiro-Wilk test on heterogeneous obs
  # # ...and report p-value in final output:
  sw.obs   = apply( het.sub, 2, shapiro.test )
  sw.obs.p = sw.obs[[1]]$p.value
  ###########################################################################

  # ###########################################################################
  # # Extract R^2 for gene.j from original regression (deconvolution):
  # # Controls:
  # fit.controls          = deconv[[ 2 ]]
  # fit.summs.controls    = summary( fit.controls )
  # fit.gene.j.controls   = fit.summs.controls[[ gene.j ]]
  # R.2.adj.control       = fit.gene.j.controls$adj.r.squared
  #
  # # Cases:
  # fit.cases             = deconv[[ 2 ]]
  # fit.summs             = summary( fit.cases )
  # fit.gene.j            = fit.summs[[ gene.j ]]
  # R.2.adj.case          = fit.gene.j$adj.r.squared
  # ###########################################################################

  # ###########################################################################
  # MSE of obs:
  resids = decon.list[[2]]
  resids.control   = resids[ groups == 1, , drop=FALSE]
  resids.case      = resids[ groups == 2, , drop=FALSE]
  resids.control.2 = resids.control ^ 2
  resids.case.2    = resids.case    ^ 2
  sse.control      = apply( resids.control.2 , 2, sum )
  sse.case         = apply( resids.case.2    , 2, sum )
  mse.control      = sse.control / ( n.control - n.cells )
  mse.case         = sse.case    / ( n.case    - n.cells )

  mse.control.frame = data.frame( as.character(names(mse.control)), mse.control, stringsAsFactors = F  )
  mse.case.frame    = data.frame( as.character(names(mse.case)), mse.case, stringsAsFactors = F)
  colnames(mse.control.frame) = c("site", "mse.control")
  colnames(mse.case.frame)    = c("site", "mse.case")
  ###########################################################################

  ###############################################################################
  # Do power analysis:
  power.list = do.decon.power( decon.list, groups, direction, nonNeg )
  # power.list = list( power, diffs,  diff.critical, base.expr, case.expr, do.not.trust.these, deconv, tail.95.1 )
  ###############################################################################

  #############################################################################
  # Assemble data.frame for output:
  # power           = power.list[[1]]
  # diffs           = power.list[[2]]
  # diff.critical   = power.list[[3]]
  # base.expr       = power.list[[4]]
  # case.expr       = power.list[[5]]
  # do.not.trust    = power.list[[6]]
  # deconv          = power.list[[7]]
  # tail.95.1       = power.list[[8]]

  power.t  = as.data.frame(t( power.list[[1]] ))
  diffs.t  = as.data.frame(t( power.list[[2]] ))
  crit.t   = as.data.frame(t( power.list[[3]] ))
  base.t   = as.data.frame(t( power.list[[4]] ))
  case.t   = as.data.frame(t( power.list[[5]] ))
  # do.not.t = as.data.frame(t( ! power.list[[6]] )) # Inverting boolean here to be intuitive: trust=TRUE means power > .5 likely.
  tail.1   = as.data.frame(t( power.list[[8]] ))

  # if(dim(het.sub)[2] > 1){
  #   power.t = t(power.t)
  #   diffs.t = t(diffs.t)
  #   crit.t  = t(crit.t)
  #   base.t  = t(base.t)
  #   case.t  = t(case.t)
  #   tail.1  = t(tail.1)
  # }
  
  power.t$site  = rownames( power.t  )
  diffs.t$site  = rownames( diffs.t  )
  crit.t$site   = rownames( crit.t   )
  base.t$site   = rownames( base.t   )
  case.t$site   = rownames( case.t   )
  # do.not.t$site = rownames( do.not.t )
  tail.1$site   = rownames( tail.1   )

  cell.names = colnames(cell.props)
  numcells   = length(cell.names)

  total.frame = data.frame()
  for( p in 1:numcells){
    tmp.frame = data.frame()
    cell.name = cell.names[p]

    text2parse = paste0("power.tmp = subset(power.t, select=c( site ,", cell.name ," ))" )
    eval(parse(text=text2parse));
    colnames(power.tmp) = c( "site", "power")
    text2parse = paste0("diffs.tmp = subset(diffs.t, select=c(", cell.name ,", site ))")
    eval(parse(text=text2parse));
    colnames(diffs.tmp) = c("diff.est", "site")

    text2parse = paste0("tail.tmp = subset(tail.1, select=c(", cell.name ,", site ))")
    eval(parse(text=text2parse));
    colnames(tail.tmp) = c("95.tail.1", "site")

    text2parse = paste0("crit.tmp = subset(crit.t, select=c(", cell.name ,", site ))")
    eval(parse(text=text2parse));
    colnames(crit.tmp) = c("crit.dif", "site")
    text2parse = paste0("base.tmp = subset(base.t, select=c( site , ", cell.name ,"  ))")
    eval(parse(text=text2parse));
    colnames(base.tmp) = c( "site", "base")
    text2parse = paste0("case.tmp = subset(case.t, select=c(", cell.name ,", site ))")
    eval(parse(text=text2parse));
    colnames(case.tmp) = c("case", "site")
    # text2parse = paste0("do.not.tmp = subset(do.not.t, select=c(", cell.name ,", site ))")
    # eval(parse(text=text2parse));
    # colnames(do.not.tmp) = c("trust", "site")

    tmp.frame = left_join( base.tmp   , case.tmp          , by="site" )
    tmp.frame = left_join( tmp.frame  , diffs.tmp         , by="site" )
    tmp.frame = left_join( tmp.frame  , tail.tmp          , by="site" )
    tmp.frame = left_join( tmp.frame  , crit.tmp          , by="site" )
    tmp.frame = left_join( tmp.frame  , mse.control.frame , by="site" )
    tmp.frame = left_join( tmp.frame  , mse.case.frame    , by="site" )
    tmp.frame = left_join( tmp.frame  , power.tmp         , by="site" )
    # tmp.frame = left_join( tmp.frame  , do.not.tmp        , by="site" )
    tmp.frame$cell    = cell.name

    cat("cell name: ", cell.name, "\n")

    tmp.frame$cell.sd = cell.SDs[p]
    tmp.frame$kappa.1 = kappa.control
    tmp.frame$kappa.2 = kappa.case

    total.frame = rbind(total.frame, tmp.frame)
  }
  #############################################################################

  #############################################################################
  # Write total.frame to .CSV file:  Places content in current working directory.
  file2output = paste0( output.file, ".csv"  )
  write.csv( total.frame, file=file2output , row.names = F )
  #############################################################################

  #############################################################################
  # Package the return list object:
  args.used = list(  output.file
                    , medCntr
                    , stdz
                    , nonNeg
                    , method
                    , direction
                    )

  return.list = list( total.frame, args.used )
  return( return.list )
  #############################################################################
  # } # If output.dir ! exists

  } # output file exists
} # End LRCDE
