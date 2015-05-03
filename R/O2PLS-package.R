#' O2PLS: Two-Way Orthogonal Partial Least Squares
#'
#' This is based on work of (Trygg, Wold, 2003). 
#' Includes some misc functions, the O2PLS fit function and some cross-validation tools.
#' 
#' @section O2PLS functions:
#' The O2PLS fit is done with \code{\link{o2m}}. Cross-validation is done with \code{\link{loocv}} or \code{\link{adjR2}} which may depend on the \code{parallel} package.
#' 
#' List of functions:
#' 
#' @author 
#' Said el Bouhaddani \email{s.el_bouhaddani@@lumc.nl}
#' 
#' Maintainer: Said el Bouhaddani \email{s.el_bouhaddani@@lumc.nl}
#' @docType package
#' @name O2PLS-package
#' @keywords O2PLS-package
#' @import parallel
