#' het.mtx
#' 
#' A subset of heterogeneous gene expression measures
#' @format A matrix of 150 samples x 490 genes, log2 scale data
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391}
"het.mtx"

#' het.pheno
#' 
#' An indicator vector of healthy/disease state
#' @format A factor of length 150, 72 healthy and 78 SLE subjects
"het.pheno"

#' cell.sigs
#' 
#' A subset of the LM22 cell signatures
#' @format A matrix of 22 cell types x 490 genes, log2 scale data
#' 
#' @source \url{https://cibersort.stanford.edu/}
"cell.sigs"

#' cell.props
#' 
#' Estimated cell proportions. Estimated using LM22 cell signatures
#' @format A matrix of 150 samples X 22 cell types.
#' 
"cell.props"

#' LM22
#' 
#' The LM22 cell signatures
#' @format A data frame of 547 genes x gene symbol + 22 cell types, decimal scale data
#' 
#' @source \url{https://cibersort.stanford.edu/}
"LM22"

