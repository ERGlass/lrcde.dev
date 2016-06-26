#' LRCDE package for cell type-specific deconvolution
#' 
#' @description 
#' lrcde allows estimation of cell type-specific genetic signal and subsequent
#' differential expression detection.
#' 
#' @details
#' lrcde differs from other algorithms in that it does not perform permutations
#' but rather a cell type by cell type gene-by-gene power calculation for
#' each detected difference.  Differences detected with relatively high power,
#' about 0.8, likely represent actual cell type-specific differential
#' expression.  Differences detected with lower power could potentially be
#' actually detected differences, but the possiblity exists that they are not.
#' 
#' The \code{\link{lrcde}} function is what you will call directly when your
#' inputs consist of a matrix of heterogeneous measures and measures of
#' relative cell proportions for all samples.
#' 
#' @examples
#' 
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
"_PACKAGE"
#> [1] "_PACKAGE"