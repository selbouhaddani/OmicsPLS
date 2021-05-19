#' Data integration with O2PLS: Two-Way Orthogonal Partial Least Squares
#'
#' The OmicsPLS package is an R package for penalized integration of heterogeneous omics data. 
#' The software articles are published in (el Bouhaddani et al, 2018, \doi{10.1186/s12859-018-2371-3}) and (Gu et al, 2020, \doi{10.1186/s12859-021-03958-3}).
#' OmicsPLS includes the O2PLS fit, the GO2PLS fit, cross-validation tools and some misc functions.
#' 
#' @author
#' Said el Bouhaddani (\email{s.elbouhaddani@@umcutrecht.nl}, Twitter: @@selbouhaddani),
#' Zhujie Gu, 
#' Szymon Kielbasa,
#' Geurt Jongbloed,
#' Jeanine Houwing-Duistermaat,
#' Hae-Won Uh.
#'
#' Maintainer: Said el Bouhaddani (\email{s.elbouhaddani@@umcutrecht.nl}).
#' 
#' @section Model and assumptions:
#' \strong{Note that the rows of \code{X} and \code{Y} are the subjects and columns are variables.}
#' The number of columns may be different, but the subjects should be the same in both datasets.
#' 
#' The O2PLS model (Trygg & Wold, 2003) decomposes two datasets \eqn{X} and \eqn{Y} into three parts. 
#' \itemize{
#'  \item{1.} A joint part, representing the relationship between \eqn{X} and \eqn{Y}
#'  \item{2.} An orthogonal part, representing the unrelated latent variation in \eqn{X} and \eqn{Y} separately.
#'  \item{3.} A noise part capturing all residual variation.
#' }
#' 
#' See also the corresponding paper (el Bouhaddani et al, 2018).
#' 
#' 
#' @section Fitting:
#' The O2PLS fit is done with \code{\link{o2m}}. 
#' For data \code{X} and \code{Y} you can run \code{o2m(X,Y,n,nx,ny)} for an O2PLS fit with \code{n} joint and \code{nx, ny} orthogonal components.
#' See the help page of \code{\link{o2m}} for more information on parameters.
#' There are four ways to obtain an O2PLS fit, depending on the dimensionality.
#' \itemize{
#'  \item{} For the not-too-high dimensional case, you may use \code{\link{o2m}} with default parameters. E.g. \code{o2m(X,Y,n,nx,ny)}.
#'  \item{} In case you only want the parameters, you may add \code{stripped = TRUE} to obtain a stripped version of \code{o2m} which avoids calculating and storing some matrices. E.g. \code{o2m(X,Y,n,nx,ny,stripped=TRUE)}.
#'  \item{} For high dimensional cases, defined by \code{ncol(X)>p_thresh} and \code{ncol(Y)>q_thresh}, a NIPALS approach is used which avoids storing large matrices. E.g. \code{o2m(X,Y,n,nx,ny,p_thresh=3000,q_thresh=3000)}.
#'  The thresholds are by default both at 3000 variables.
#'  \item{} If you want a stripped version in the high dimensional case, add \code{stripped = TRUE}. E.g. \code{o2m(X,Y,n,nx,ny,stripped=TRUE,p_thresh=3000,q_thresh=3000)}.
#'  \item{} For GO2PLS, add \code{sparsity = TRUE} and specify how many variables or groups to retain. E.g. \code{o2m(X,Y,n,nx,ny,sparse=TRUE,keepx, keepy)}.
#' }
#' 
#' @section Obtaining results:
#' After fitting an O2PLS model, by running e.g. \code{fit = o2m(X,Y,n,nx,ny)}, the results can be visualised.
#' Use \code{\link{plot}(fit,...)} to plot the desired loadings with/without ggplot2.
#' Use \code{\link{summary}(fit,...)} to see the relative explained variances in the joint/orthogonal parts.
#' Also plotting the joint scores \code{fit$Tt, fit$U} and orthogonal scores \code{fit$T_Yosc, fit$U_Xosc} are of help.
#' 
#' @section Cross-validating: 
#' Determining the number of components \code{n,nx,ny} is an important task. For this we have two methods.
#' See \code{citation("OmicsPLS")} for our proposed approach for determining the number of components, implemented in \code{crossval_o2m_adjR2}!
#' \itemize{
#'  \item{} Cross-validation (CV) is done with \code{\link{crossval_o2m}} and \code{\link{crossval_o2m_adjR2}}, both have built in parallelization which relies on the \code{parallel} package.
#'  Usage is something like \code{crossval_o2m(X,Y,a,ax,ay,nr_folds)} where \code{a,ax,ay} are vectors of integers. See the help pages.
#'  \code{nr_folds} is the number of folds, with \code{nr_folds = nrow(X)} for Leave-One-Out CV.
#'  \item{} For \code{crossval_o2m_adjR2} the same parameters are to be specified. This way of cross-validating is (potentially much)
#'  faster than the standard approach. It is also recommended over the standard CV.
#'  \item{} To cross-validate the number of variables to keep, use \code{\link{crossval_sparsity}}. 
#' }
#' 
#' @section S3 methods:
#' There are S3 methods implemented for a fit obtained with \code{o2m}, i.e. \code{fit <- o2m(X,Y,n,nx,ny)}
#' \itemize{
#'   \item{} Use plot(fit) to plot the loadings, see above.
#'   \item{} Use \code{\link{loadings}(fit)} to extract a matrix with loading values
#'   \item{} Use \code{\link{scores}(fit)} to extract the scores
#'   \item{} Use \code{\link{print}} and \code{\link{summary}} to print and summarize the fit object
#' }
#' 
#' @section Imputation:
#' When the data contains missing values, one should impute them prior to using O2PLS.
#' There are many sophisticated approaches available, such as MICE and MissForest, and no one approach is the best for all situations.
#' To still allow users to quickly impute missing values in their data matrix, 
#' the \code{\link{impute_matrix}} function is implemented. 
#' It relies on the \code{\link{softImpute}} function+package and imputes based on the singular value decomposition.
#' 
#' @section Misc:
#' Also some handy tools are available
#' \itemize{
#'  \item{} \code{\link{orth}(X)} is a function to obtain an orthogonalized version of a matrix or vector \code{X}.
#'  \item{} \code{\link{ssq}(X)} is a function to calculate the sum of squares (or squared Frobenius norm) of \code{X}. See also \code{\link{vnorm}} for calculating the norm of each column in \code{X}.
#'  \item{} \code{\link{mse}(x, y)} returns the mean squared difference between two matrices/vectors. 
#' }
#' 
#' @section Citation:
#' If you use the OmicsPLS R package in your research, please cite the corresponding software paper:
#' 
#' \strong{el Bouhaddani, S., Uh, H.-W., Jongbloed, G., Hayward, C., Klarić, L., Kiełbasa, S. M., & Houwing-Duistermaat, J.} (2018).
#' \emph{Integrating omics datasets with the OmicsPLS package.}
#'  BMC Bioinformatics, 19(1). \doi{10.1186/s12859-018-2371-3}
#' 
#' The bibtex entry can be obtained with command \code{citation("OmicsPLS")}.
#' Thank you!
#' 
#' The original paper proposing O2PLS is
#' 
#' \strong{Trygg, J., & Wold, S.} (2003). 
#' \emph{O2-PLS, a two-block (X-Y) latent variable regression (LVR) method with an integral OSC filter.} 
#' Journal of Chemometrics, 17(1), 53-64. \doi{10.1002/cem.775}
#' 
#' @docType package
#' @name OmicsPLS
#' @keywords OmicsPLS
#' @import parallel ggplot2 tibble magrittr softImpute
#' @importFrom graphics abline
#' @importFrom stats cov sd
#' @importFrom dplyr mutate
NULL

#' Check if matrices satisfy input conditions
#'
#' @param X Should be numeric matrix.
#' @param Y Should be numeric matrix.
#' @return NULL
#' @details This function throws an error if any of the elements is \code{NA}, \code{Inf}, \code{NaN} or \code{nrow(X)} doesn't match \code{nrow(Y)}.
#' 
#' @keywords internal
#' @export
input_checker <- function(X, Y = NULL) {
  if(!is.numeric(X)) stop("Input is not numeric, but of mode ",mode(X),"\n")
  if(!is.matrix(X)) stop("Input is not a matrix, but of class ",class(X),"\n")
  if(any(is.na(X))) stop("Input contains NA's or NaN's, consider imputing with impute_matrix","\n")
  if(any(is.infinite(X))) stop("Input contains non-finite elements, consider imputing with impute_matrix","\n")
  
  if (!is.null(Y)) {
    if(!is.numeric(Y)) stop("Input is not numeric, but of mode ",mode(Y),"\n")
    if(!is.matrix(Y)) stop("Input is not a matrix, but of class ",class(Y),"\n")
    if(any(is.na(Y))) stop("Input contains NA's or NaN's, consider imputing with impute_matrix","\n")
    if(any(is.infinite(Y))) stop("Input contains non-finite elements, consider imputing with impute_matrix","\n")
    if(nrow(X) != nrow(Y)) stop("# rows don't match: ",nrow(X)," versus ",nrow(Y),"\n")
    if(!identical(rownames(X), rownames(Y))) warning("Caution: Rownames don't match!","\n")
  }
  NULL
}

#' Impute missing values in a matrix
#'
#' @param X A matrix with missing values in some entries.
#' @param ... Further arguments for \code{softimpute}.
#' @return An imputed version of matrix \eqn{X}
#' @details This function is based on the \code{\link{softImpute}} function in its eponymous package.
#' @examples
#' X <- matrix(rnorm(20*100),20)
#' Xmis <- X
#' Xmis[sample(length(Xmis),length(Xmis)/10)] <- NA
#' anyNA(X)
#' anyNA(impute_matrix(Xmis))
#' @export
impute_matrix <- function(X, ...){
  Xnames = dimnames(X)
  if(!is.matrix(X)){
    message("X has class ",class(X),", trying to convert with as.matrix.",sep="","\n")
    X <- as.matrix(X)
    dimnames(X) <- Xnames
  }
  if(any(is.infinite(X))) X[is.infinite(X)] = NA
  if(!anyNA(X)){
    message("X doesn't contain missings. Returning original matrix.","\n")
    return(X)
  }
  
  imp <- softImpute::softImpute(X, rank.max = min(dim(X),51)-1, 
                                type="svd", ...)
  return(softImpute::complete(X, imp))
}


#' Orthogonalize a matrix
#'
#' @param X Numeric vector or matrix.
#' @param X_true (optional) A 'true' matrix/vector. Used to correct the sign of the orthonormalized X if QR is used. Only the first column is corrected.
#' @param type A character or numeric. Should be one of "QR" or "SVD".
#' @return An orthogonalized representation of \eqn{X}
#' @details Choosing type='QR' uses a QR decomposition of X to produce orthonormal columns. For type=='SVD' it uses an SVD decomposition.
#' The columns are corrected for sign.
#' @examples
#' orth(c(3,4))
#' round(crossprod(orth(matrix(rnorm(500),100,5))),4)
#' orth(matrix(1:9,3,3),type='QR')[,1] - orth(1:3); orth(matrix(1:9,3,3),type='SVD')[,1] - orth(1:3);
#' @export
orth <- function(X, X_true = NULL, type = c("QR", "SVD")) {
  X <- as.matrix(X)
  if (!is.null(X_true)) {
    X_true <- as.matrix(X_true)
    input_checker(X,X_true)
    if(ncol(X) != ncol(X_true)) stop("# columns don't match:",ncol(X),"versus",ncol(X_true),"\n")
  }else {
    input_checker(X)
  }
  
  type <- match.arg(type)
  if (type == "SVD") {
    e <- svd(X)
    e <- tcrossprod(e$u, e$v)
  } else {
    e <- qr.Q(qr(X))
  }
  if(is.null(X_true)) {
    sign_e <- sign(crossprod(e,X)) * diag(1,ncol(e))
  } else {
    sign_e <- sign(crossprod(e,X_true)) * diag(1,ncol(e))
  }
  if(any(diag(sign_e)==0)){warning("Orthogonalization made some columns orthogonal to original columns","\n")}
  
  return(e %*% sign_e)
}

#' Calculate Sum of Squares
#'
#' @param X Numeric vector or matrix.
#' @return The sum of squared elements of \eqn{X}
#' @details This is the Frobenius norm of \eqn{X}.
#' @examples
#' ssq(tcrossprod(1:5))
#' ssq(rnorm(1e5))/1e5
#' @export
ssq <- function(X) {
  return(sum(X^2))
}

#' Calculate mean squared difference
#'
#' @param x Numeric vector or matrix.
#' @param y Numeric vector or matrix. Defaults to 0.
#' @param na.rm Remove NA's?
#' @return The mean of the squared differences elementwise.
#' @details Is equal to ssq(\code{x-y})/length(c(\code{x})). If \code{x} and \code{y} are of unequal length, the invoked minus-operator will try to make the best out of it by recycling elements of the shorter object (usually you don't want that).
#' In particular if \code{x} is an N x p matrix and \code{y} an N x 1 vector, y is subtracted from each column of \code{x}, and if \code{y=0} (default) you get the mean of vec(\code{x^2})
#' @examples
#' mse(2)
#' mse(1:10,2:11) == 1
#' mse(matrix(rnorm(500),100,5),matrix(rnorm(500),100,5))
#' @export
mse <- function(x, y = 0, na.rm = FALSE)
{
  if(length(x) != length(y) && length(y) != 0) message("Comparing lengths ",length(x)," with ",length(y),"\n")
  mean((x - y)^2, na.rm = na.rm)
}

#' Root MSE of Prediction
#'
#' Calculates the Root MSE of prediction on test data. Only tested to work inside \code{\link{loocv}}.
#'
#' @param Xtst Numeric vector or matrix.
#' @param Ytst Numeric vector or matrix.
#' @param fit \code{\link{o2m}} fit (on data without \code{Xtst} and \code{Ytst}).
#' @param combi Logical. Should the symmetrized MSE be used, i.e. add both MSEs. Not yet implemented, but see \code{\link{rmsep_combi}} 
#' @details This function is the building block for \code{\link{loocv}}, as it produced the prediction error for test (left out) data.
#' @return Mean squares difference between predicted Y and true Y
#' @export
rmsep <- function(Xtst, Ytst, fit, combi = FALSE) {
  if(!inherits(fit,c("o2m","pre.o2m"))) stop("fit should be an O2PLS fit","\n")
  
  if (!is.matrix(Xtst)) Xtst <- t(Xtst)
  
  if (!is.matrix(Ytst)) Ytst <- t(Ytst)
  
  input_checker(Xtst, Ytst)
  
  Yhat <- Xtst %*% fit$W. %*% fit$B_T %*% t(fit$C.)
  # Xhat = Ytst%*%fit$C.%*%fit$B_U%*%t(fit$W.)
  
  return(mean(c(sqrt(mse(Yhat, Ytst)))))  #,sqrt(mse(Xhat,Xtst)))))
}

#' K fold CV for O2PLS
#'
#' For (a grid of) values for \code{a}, \code{nx} and \code{ny}, \code{loocv} estimates the prediction error using k-fold CV.
#'
#' @inheritParams o2m
#' @param a Vector of integers. Contains the numbers of joint components.
#' @param a2 Vector of integers. Contains the numbers of orthogonal components in \eqn{X}.
#' @param b2 Vector of integers. Contains the numbers of orthogonal components in \eqn{Y}.
#' @param fitted_model List. Deprecated. O2PLS model fit with \code{\link{o2m}}. Is used to calculate the apparent error without recalculating this fit.
#' @param func Function to fit the O2PLS model with. Only \code{\link{o2m}} and \code{\link{o2m_stripped}} are supported.
#' @param app_err Logical. Deprecated. Should the apparent error also be computed? 
#' @param kcv Integer. The value of \eqn{k}, i.e. the number of folds. Choose \eqn{N} for LOO-CV.
#' @details Note that this function can be easily parallelized (on Windows e.g. with the \code{parallel} package.).
#' @return List with two numeric vectors:
#' \item{CVerr}{Contains the k-fold CV estimated RMSEP}
#' \item{Fiterr}{Contains the apparent error}
#' @details The parameters \code{a}, \code{a2} and \code{b2} can be integers or vectors of integers. A for loop is used to loop over all combinations.
#' The resulting output is a list, which is more easy to interpret if you use \code{array(unlist(output_of_loocv$CVerr))} as in the example below.
#' The array wil have varying \code{a} along the first dimension and \code{a2} and \code{b2} along the second and third respectively.
#' Typing \code{example(loocv)} (hopefully) clarifies this function.
#' @export
loocv <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, fitted_model = NULL, func = o2m, app_err = F, kcv,
                  stripped = TRUE, p_thresh = 3000, 
                  q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  app_err = F
  fitted_model = NULL
  stopifnot(all(a == round(a)), all(a2 == round(a2)), all(b2 == round(b2)))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  # if (!is.null(fitted_model)) {
  #   app_err <- F
  #   message("apparent error calculated with provided fit","\n")
  # }
  # determine type of model
  type <- 3  #ifelse(deparse(substitute(func))=='o2m',3,ifelse(deparse(substitute(func))=='oplsm',2,1))
  
  N <- nrow(X)
  if (N != nrow(Y)) {
    stop("N not the same","\n")
  }
  mean_err <- mean_fit <- NA * 1:max(length(a), length(a2), length(b2))
  k <- 0
  
  # blocks contains the begin and endpoints of test indices to use
  blocks <- cut(seq(1:N), breaks=kcv, labels=F)
  
  # loop through chosen parameters
  for (j in a) {
    for (j2 in a2) {
      for (j3 in b2) {
        k <- k + 1
        err <- NA * 1:kcv
        folds <- sample(N)
        # loop through number of folds
        for (i in 1:kcv) {
          ii <- which(blocks==i)
          if (type == 3) {
            pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = j, nx = j2, ny = j3,
                         stripped = stripped, p_thresh = p_thresh, 
                         q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)
          }
          # if(type==2){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j,n_orth=j2)}
          # if(type==1){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j)}
          fit <- try(do.call(func, pars), silent = T)
          err[i] <- ifelse(inherits(fit, "try-error"), NA, rmsep(X[folds[ii], ], Y[folds[ii], ], 
                                                                 fit))
        }
        mean_err[k] <- mean(err)
        # calculate apparent error
        # if (app_err && is.null(fitted_model)) {
        #   if (type == 3) {
        #     pars2 <- list(X = X, Y = Y, n = j, nx = j2, ny = j3)
        #   }
        #   # if(class(fit)=='oplsm'){pars2=list(X=X,Y=Y,ncomp=j,n_orth=j2)}
        #   # if(class(fit)=='plsm'){pars2=list(X=X,Y=Y,ncomp=j)}
        #   fit2 <- try(do.call(func, pars2), F)
        #   mean_fit[k] <- ifelse(inherits(fit, "try-error"), NA, rmsep(X, Y, fit2))
        #   # print('1e loop')
        # }
        # if (!is.null(fitted_model)) {
        #   mean_fit[k] <- rmsep(X, Y, fitted_model)
        #   # print('2e loop')
        # }
      }
    }
  }
  return(list(CVerr = mean_err, Fiterr = mean_fit))
}

#' Gridwise adjusted R2 for O2PLS
#'
#' For (a grid of) values for \code{a}, \code{nx} and \code{ny}, \code{loocv} calculates the R2 of the joint part. Parallel computing is supported on Windows with package \code{parallel}.
#'
#' @inheritParams o2m
#' @param a Vector of integers. Contains the numbers of joint components.
#' @param a2 Vector of integers. Contains the numbers of orthogonal components in \eqn{X}.
#' @param b2 Vector of integers. Contains the numbers of orthogonal components in \eqn{Y}.
#' @param func Function to fit the O2PLS model with. Only \code{\link{o2m}} and \code{\link{o2m_stripped}} are supported.
#' @param parall Integer. Should a parallel cluster be set up using package \code{parallel} (Windows)? Best is to leave it to \code{FALSE}.
#' @param cl Object of class '\code{cluster}'. If parall is \code{TRUE} and \code{cl} is not \code{NULL}, calculations are parallelized over workers in cl.
#' @details The use of this function is to calculate the R2 of the joint part, while varying the number of orthogonal components. Adding more joint components will increase the R2!
#'
#' A parallelized version is built in -tested on windows-, use package \code{parallel} and set \code{parall=TRUE} to activate this. There should not be already a cluster object with the name \code{cl}.
#' In case of some error, don't forget to invoke \code{stopCluster(cl)} to end the cluster. See Task Manager (Windows) to verify that the workers are spanned/ended.
#' @return Matrix with two rows:
#' \item{adjR2X}{Contains the joint R2 in X}
#' \item{adjR2Y}{Contains the joint R2 in Y}
#' @details See \code{\link{loocv}} for more intuition.
#' @export
adjR2 <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, func = o2m, parall = F, cl = NULL, 
                  stripped = TRUE, p_thresh = 3000, 
                  q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  stopifnot(all(a == round(a)), all(a2 == round(a2)), all(b2 == round(b2)))
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_checker(X, Y)
  cl_was_null <- FALSE
  if (!parall) {
    S_apply <- function(cl = NULL, x, fun) {
      sapply(x, fun)
    }
  }
  if (parall & is.null(cl)) {
    cl_was_null <- TRUE
    S_apply <- parSapply
    cl <- makeCluster(rep("localhost", detectCores()), type = "SOCK")
    clusterExport(cl = cl, varlist = c("ssq", "o2m_stripped", "adjR2"))
  }
  if (parall & !is.null(cl)) {
    stopifnot(inherits(cl,'cluster'))
    S_apply <- parSapply
  }
  
  pars1 <- merge(merge(data.frame(a = a), data.frame(a2 = a2)), data.frame(b2 = b2))
  pars2 <- apply(pars1, 1, as.list)
  
  # cl <- makeCluster(rep( 'localhost', detectCores()),type='SOCK') clusterExport(cl=cl,
  # varlist=c('X','Y','N','pars2','ssq','o2m'))
  outp <- S_apply(cl, pars2, function(p) {
    fit <- try(func(X, Y, p$a, p$a2, p$b2, stripped = stripped, p_thresh = p_thresh, 
                    q_thresh = q_thresh, tol = tol, max_iterations = max_iterations), silent=T)
    if(inherits(fit,'try-error')) {
      if(!grepl("exceeds #columns in X or Y", fit)){
        warning(fit[1])
      }
    }
    if(!inherits(fit, 'try-error')){
      RX <- 1 - ssq(fit$H_UT)/ssq(fit$U)
      RY <- 1 - ssq(fit$H_TU)/ssq(fit$Tt)
      adjRX <- RX  #1 - (1 - RX)*(N - 1)/(N - p$a - 1)
      adjRY <- RY  #1 - (1 - RY)*(N - 1)/(N - p$a - 1)
      return(c(adjR2X = adjRX, adjR2Y = adjRY))
    } else {
      return(c(adjR2X = NA, adjR2Y = NA))
    }
  })
  if (parall & cl_was_null == TRUE) {
    stopCluster(cl)
  }
  return(outp)
}

#' Norm of a vector or columns of a matrix
#'
#' @param x Numeric vector or matrix.
#' @return (columnwise) Euclidian norm of \eqn{x}
#' @examples
#' vnorm(orth(1:5))
#' vnorm(matrix(1:9,3,3))^2 - colSums(matrix(1:9,3)^2)
#' @export
vnorm <- function(x)
{
  x <- as.matrix(x)
  return(sqrt(apply(x^2, 2, sum)))
}

#' Symmetrized root MSE of Prediction
#'
#' Calculates the symmetrized root MSE of prediction on test data. *Expected* to work in combination with \code{\link{loocv}}.
#'
#' @param Xtst Numeric vector or matrix.
#' @param Ytst Numeric vector or matrix.
#' @param fit \code{\link{o2m}} fit (on data without \code{Xtst} and \code{Ytst}).
#' @details This function is the building block for \code{\link{loocv}}, as it produced the prediction error for test (left out) data.
#'
#' This is a symmetrized version of \code{\link{rmsep}}, and is used by \code{\link{loocv}}. The predicion error of both \code{Xtst} and \code{Ytst} are calculated and summed.
#' Whether this is a good idea depends: If \eqn{X} and \eqn{Y} have similar meanings (LC-MS versus MALDI) this is a good thing to do. If the two matrices do not have similar meanings,
#' (Metabolomics versus Transcriptomics) then you may want to not sum up the two prediction errors or include weights in the sum.
#' @return Mean squares difference between predicted Y and true Y
#' @export
rmsep_combi <- function(Xtst, Ytst, fit)
{
  if(!inherits(fit,c("o2m","pre.o2m"))) stop("fit should be an O2PLS fit","\n")
  
  if (!is.matrix(Xtst)) Xtst <- t(Xtst)
  
  if (!is.matrix(Ytst)) Ytst <- t(Ytst)
  
  input_checker(Xtst, Ytst)
  
  Yhat <- (Xtst - Xtst %*% fit$W_Yosc %*% t(fit$W_Yosc)) %*% fit$W. %*% fit$B_T %*% t(fit$C.)
  Xhat <- (Ytst - Ytst %*% fit$C_Xosc %*% t(fit$C_Xosc)) %*% fit$C. %*% fit$B_U %*% t(fit$W.)
  
  return(sqrt(mse(Yhat, Ytst)) + sqrt(mse(Xhat, Xtst)))
}

#' K-fold CV based on symmetrized prediction error
#'
#' The prediction error of both \code{X~Xhat} and \code{Y~Yhat} are summed. This provides a symmetrized version of \code{\link{loocv}}.
#' @inheritParams o2m
#' @inheritParams loocv
#' @details Note that this function can be easily parallelized (on Windows e.g. with the \code{parallel} package.).
#' If there are NAs in the CVerr component, this is due to an error in the fitting.
#' @return List with two numeric vectors:
#' \item{CVerr}{Contains the k-fold CV estimated RMSEP}
#' \item{Fiterr}{Contains the apparent error}
#'
#' @export
loocv_combi <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, fitted_model = NULL, func = o2m, app_err = F, kcv,
                        stripped = TRUE, p_thresh = 3000, 
                        q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  app_err = F
  fitted_model = NULL
  stopifnot(all(a == round(a)), all(a2 == round(a2)), all(b2 == round(b2)), kcv == round(kcv[1]))
  stopifnot(is.logical(app_err), is.function(func))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  # if (!is.null(fitted_model)) {
  #   if(inherits(fitted_model,c("o2m","pre.o2m"))){stop("fitted_model should be of class 'o2m' or NULL","\n")}
  #   app_err <- F
  #   message("apparent error calculated with provided fit")
  # }
  
  # determine type of model, deprecated 
  type <- 3  #ifelse(deparse(substitute(func))=='o2m',3,ifelse(deparse(substitute(func))=='oplsm',2,1))
  
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  mean_err <- mean_fit <- NA * 1:max(length(a), length(a2), length(b2))
  k <- 0
  
  # blocks contains the begin and endpoints of test indices to use
  # Extra N ignored by 1:kcv if N/kcv is integer
  blocks <- cut(seq(1:N), breaks=kcv, labels=F)
  
  # loop through chosen parameters
  for (j in a) {
    for (j2 in a2) {
      for (j3 in b2) {
        k <- k + 1
        err <- NA * 1:kcv
        folds <- sample(N)
        # loop through number of folds
        for (i in 1:kcv) {
          ii <- which(blocks==i)
          if (type == 3) {
            pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = j, nx = j2, ny = j3, 
                         stripped = stripped, p_thresh = p_thresh, 
                         q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)
          }
          fit <- try(do.call(func, pars), silent = T)
          if(inherits(fit,'try-error')) {
            if(!grepl("exceeds #columns in X or Y", fit)){
              warning(fit[1])
            }
          }
          err[i] <- ifelse(inherits(fit, 'try-error'), 
                           NA, 
                           rmsep_combi(X[folds[ii], ], Y[folds[ii], ], fit))
        }
        mean_err[k] <- mean(err)
        # calculate apparent error
        # if (app_err && is.null(fitted_model)) {
        #   if (inherits(fit,c("o2m","pre.o2m"))) {
        #     pars2 <- list(X = X, Y = Y, n = j, nx = j2, ny = j3)
        #   }
        #   # if (class(fit) == "oplsm") {
        #   #   pars2 <- list(X = X, Y = Y, ncomp = j, n_orth = j2)
        #   # }
        #   # if (class(fit) == "plsm") {
        #   #   pars2 <- list(X = X, Y = Y, ncomp = j)
        #   # }
        #   fit2 <- try(do.call(func, pars2), F)
        #   mean_fit[k] <- ifelse(inherits(fit,"try-error"), NA, rmsep_combi(X, Y, fit2))
        #   print("1e loop")
        # }
        # if (!is.null(fitted_model)) {
        #   mean_fit[k] <- rmsep_combi(X, Y, fitted_model)
        #   print("2e loop")
        # }
      }
    }
  }
  return(list(CVerr = mean_err, Fiterr = mean_fit))
} 


# Generic Methods ---------------------------------------------------------

#' Print function for O2PLS.
#' 
#' This function is the print method for an O2PLS fit
#' 
#' @param x An O2PLS fit (an object of class o2m)
#' @param ... For consistency
#' @return NULL
#'
#'
#' @export
print.o2m <- function (x, ...) {
  # Diagnostics or something?
  # Time to end
  # Used stripped method or high dimensional method
  # ...
  n = x$flags$n #ncol(x$W.)
  nx = x$flags$nx #ifelse(vnorm(x$P_Yosc.)[1] == 0, 0, ncol(x$P_Yosc.))
  ny = x$flags$ny #ifelse(vnorm(x$P_Xosc.)[1] == 0, 0, ncol(x$P_Xosc.))
  if(x$flags$method == "SO2PLS") cat("SO2PLS fit \n")
  else if(x$flags$method == "GO2PLS") cat("GO2PLS fit \n")
  else{
    if(x$flags$stripped) cat("O2PLS fit: Stripped \n") 
    else if(x$flags$highd) cat("O2PLS fit: High dimensional \n") 
    else cat("O2PLS fit \n")
  }

  cat("with ",n," joint components  \n",sep='')
  cat("and  ",nx," orthogonal components in X \n",sep='')
  cat("and  ",ny," orthogonal components in Y \n",sep='')
  cat("Elapsed time: ",x$flags$time, " sec\n\n", sep='')
}

#' Print function for O2PLS.
#' 
#' This function is the print method for an O2PLS fit
#' 
#' @param x An O2PLS fit (an object of class o2m)
#' @param ... For consistency
#' @return NULL
#'
#'
#' @export
print.pre.o2m <- function (x, ...) {
  cat("Internal function used to fit O2PLS, use o2m to enable print, plot, etc \n")
}


#' Plot one or two loading vectors for class o2m
#' 
#' This function plots one or two loading vectors, by default with ggplot2. 
#' 
#' @param x An O2PLS fit, with class 'o2m'
#' @param loading_name character string. One of the following: 'Xjoint', 'Yjoint', 'gr_Xjoint', 'gr_Yjoint', 'Xorth' or 'Yorth'.
#' @param i Integer. First component to be plotted.
#' @param j NULL (default) or Integer. Second component to be plotted.
#' @param use_ggplot2 Logical. Default is \code{TRUE}. If \code{FALSE}, the usual plot device will be used.
#' @param label Character, either 'number' or 'colnames'. The first option prints numbers, the second prints the colnames
#' @param ... Further arguments to \code{geom_text}, such as size, col, alpha, etc.
#' 
#' @return If \code{use_ggplot2} is \code{TRUE} a ggplot2 object. Else NULL.
#' 
#' @seealso \code{\link{summary.o2m}}
#' 
#' @export
plot.o2m <- function (x, loading_name = c("Xjoint", "Yjoint", "gr_Xjoint", "gr_Yjoint", "Xorth", "Yorth"), i = 1, j = NULL, use_ggplot2=TRUE, label = c("number", "colnames"), ...)
{
  stopifnot(i == round(i), is.logical(use_ggplot2))
  
  fit <- list()
  loading_name = match.arg(loading_name)
  if((loading_name %in% c("gr_Xjoint", "gr_Yjoint")) & x$flags$method != "GO2PLS") stop("Loading plots at group level only available in GO2PLS")
  
  which_load = switch(loading_name, Xjoint = "W.", Yjoint = "C.", 
                      gr_Xjoint = "W_gr", gr_Yjoint = "C_gr", Xorth = "P_Yosc.", Yorth = "P_Xosc.")
  fit$load = as.matrix(x[which_load][[1]])
  if(ncol(fit$load) < max(i,j) )
    stop("i and j cannot exceed #components = ",ncol(fit$load))
  fit$load = fit$load[,c(i,j)]
  
  p = nrow(as.matrix(fit$load))
  if(is.null(j)){
    fit$load = cbind(1:p,fit$load)
    colnames(fit$load) = c("index",paste(loading_name,"loadings",i))
  }else{
    stopifnot(j == round(j))
    colnames(fit$load) = c(paste(loading_name,"loadings",i),paste(loading_name,"loadings",j))
  }
  
  label2 = match.arg(label)
  if(label2 == "colnames" && !is.null(rownames(x[which_load][[1]]))) {
    label = rownames(x[which_load][[1]])
  } else label = 1:p
  if(label2 == "colnames" && is.null(rownames(x[which_load][[1]]))) message("No labels found in colnames, proceeding...","\n")
  
  if (use_ggplot2) {
    plt = with(fit, {
      ggplot(data.frame(x = load[, 1], y = load[, 2]), aes(x = x, y = y, label = I(label))) + 
        geom_text(...) + 
        labs(x = colnames(load)[1], y = colnames(load)[2])
      })
    plt = plt + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
    #print(plt)
    return(plt)
  }
  else {
    with(fit, {
      plot(load[, 1], load[, 2], type = "n")
      text(load[, 1], load[, 2])
    })
    abline(v=0,h=0)
  }
}

#' Summary of an O2PLS fit
#'
#' Until now only variational summary given by the R2's is outputted
#'
#' @param object List. Should be of class \code{\link{o2m}}.
#' @param digits Integer, number of digits.
#' @param ... For compatibility
#' @return List with R2 values.
#' @examples
#' summary(o2m(scale(-2:2),scale(-2:2*4),1,0,0))
#' 
#' @seealso \code{\link{plot.o2m}}
#' 
#' @export
summary.o2m <- function(object, digits = 3, ...) {
  fit <- object
  a <- ncol(fit$W.)
  if(digits != round(digits) || digits <= 0) stop("Digits must be a positive integer","\n")
  outp <- with( fit, list(
    Comp = a,
    R2_X = R2X,
    R2_Y = R2Y,
    R2_Xjoint = R2Xcorr,
    R2_Yjoint = R2Ycorr,
    R2_Xhat = R2Xhat,
    R2_Yhat = R2Yhat,
    R2_Xpred = R2Xhat / R2Xcorr,
    R2_Ypred = R2Yhat / R2Ycorr,
    B_T = B_T.,
    B_U = B_U,
    flags = flags,
    digits = digits
  ) )
  class(outp) <- "summary.o2m"
#   Mname <- list(c(""), c("Comp", "R2X", "R2Y", "R2Xcorr", "R2Ycorr", "R2Xhat", "R2Yhat", "XRatio", "YRatio"))
#   M <- matrix(c(ifelse(perc,a/100,a), fit$R2X, fit$R2Y, fit$R2Xcorr, fit$R2Ycorr, fit$R2Xhat, fit$R2Yhat, fit$R2Xhat/fit$R2Xcorr, 
#                 fit$R2Yhat/fit$R2Ycorr), nrow = 1, dimnames = Mname)
#   return(round((1 + perc * 99) * M, 4 - perc * 2))
  outp
}

#' Prints the summary of an O2PLS fit
#'
#' Readable output is given in the form of percentages of variances explained.
#'
#' @inheritParams summary.o2m
#' @return NULL
#' @keywords internal
#' @examples
#' summary(o2m(scale(-2:2),scale(-2:2*4),1,0,0))
#' @export
print.summary.o2m <- function(x, ...){
  digits = x$digits
  method = x$flags$method
  cat(paste("\n*** Summary of the", method, "fit *** \n\n"))
  R2_names = c("Joint","Orthogonal","Noise")
  R2_Xall = with(x,{c(R2_Xjoint, R2_X - R2_Xjoint, 1 - R2_X)})
  R2_Yall = with(x,{c(R2_Yjoint, R2_Y - R2_Yjoint, 1 - R2_Y)})
  R2_dataframe = data.frame(X = R2_Xall, Y = R2_Yall, 
                            row.names = R2_names)
  names(R2_dataframe) <- c("data X", "data Y")
  with(x,{
    cat("-  Call:",deparse(x$flags$call),"\n\n")
    cat("-  Modeled variation\n")
    cat("-- Total variation:\n")
    cat("in X:",flags$ssqX,"\n")
    cat("in Y:",flags$ssqY,"\n\n")
    cat("-- Joint, Orthogonal and Noise as proportions:\n\n")
    print(round(R2_dataframe, digits))
    cat("\n")
    cat("-- Predictable variation in Y-joint part by X-joint part:\n")
    cat("Variation in T*B_T relative to U:",round(R2_Ypred,digits),"\n")
    cat("-- Predictable variation in X-joint part by Y-joint part:\n")
    cat("Variation in U*B_U relative to T:",round(R2_Xpred,digits),"\n")
    cat("\n")
    cat("-- Variances per component:\n\n")
    with(flags,{
      ssqdf = rbind(varXjoint,varYjoint)
      ssqdf = as.data.frame(ssqdf)
      row.names(ssqdf) <- c("X joint", "Y joint")
      names(ssqdf) <- paste("Comp", 1:n)
      print(round(ssqdf, digits))
      cat("\n")
      if(nx > 0){
        ssqdf = t(varXorth)
        ssqdf = as.data.frame(ssqdf)
        row.names(ssqdf) <- "X Orth"
        names(ssqdf) <- paste("Comp", 1:nx)
        print(round(ssqdf, digits))
        cat("\n")
      }
      if(ny > 0){
        ssqdf = t(varYorth)
        ssqdf = as.data.frame(ssqdf)
        row.names(ssqdf) <- "Y Orth"
        names(ssqdf) <- paste("Comp", 1:ny)
        print(round(ssqdf, digits))
        cat("\n")
      }
    })
    cat("\n")
    cat("-  Coefficient in 'U = T B_T + H_U' model:\n")
    cat("-- Diagonal elements of B_T =\n", round(diag(B_T),3),"\n\n")
    #     cat("-- Coefficient in 'T = U B_U + H_T' model:\n")
    #     cat("Diagonal elements of B_U =\n", diag(B_U))
  })
  NULL
}

#' Extract the loadings from an O2PLS fit
#'
#' This function extracts loading parameters from an O2PLS fit
#'
#' @param x Object of class \code{o2m}
#' @param ... For consistency
#' 
#' @return Loading matrix
#' @examples
#' loadings(o2m(scale(-2:2),scale(-2:2*4),1,0,0))
#' 
#' @seealso \code{\link{scores.o2m}}
#' 
#' @rdname loadings
#' @export
loadings <- function(x, ...) UseMethod("loadings")


#' @param loading_name character string. One of the following: 'Xjoint', 'Yjoint', 'gr_Xjoint', 'gr_Yjoint', 'Xorth' or 'Yorth'.
#' @param subset subset of loading vectors to be extracted.
#' @param sorted Logical. Should the rows of the loadings be sorted according to the 
#' absolute magnitute of the first column?
#' 
#' @rdname loadings
#' @export
loadings.o2m <- function(x, loading_name = c("Xjoint", "Yjoint", "gr_Xjoint", "gr_Yjoint", "Xorth", "Yorth"), 
                         subset = 0, sorted = FALSE, ...) {
  if(any(subset != abs(round(subset)))) stop("Subset must be a vector of non-negative integers","\n")
  
  loading_name = match.arg(loading_name)
  if((loading_name %in% c("gr_Xjoint", "gr_Yjoint")) & x$flags$method != "GO2PLS") stop("Loading plots at group level only available in GO2PLS")
  
  which_load = switch(loading_name, Xjoint = "W.", Yjoint = "C.", 
                      gr_Xjoint = "W_gr", gr_Yjoint = "C_gr", Xorth = "P_Yosc.", Yorth = "P_Xosc.")
  loading_matrix = x[[which_load]]
  dim_names = dimnames(loading_matrix)
  if(length(subset) == 1 && subset == 0) subset = 1:ncol(loading_matrix)
  if(max(subset) > ncol(loading_matrix)) stop("Elements in subset exceed #components","\n")
  loading_matrix = as.matrix(loading_matrix[,subset])
  dimnames(loading_matrix) <- dim_names
  
  if(sorted){
    order_load = order(loading_matrix[,1]^2, decreasing = TRUE)
    if(is.null(dim_names[[1]])) dim_names[[1]] <- order_load
    loading_matrix = loading_matrix[order_load,]
  }
  return(loading_matrix)
}

#' Extract the scores from an O2PLS fit
#'
#' This function extracts score matrices from an O2PLS fit
#'
#' @param x Object of class \code{o2m}
#' @param ... For consistency
#' 
#' @return Scores matrix
#' @examples
#' scores(o2m(scale(-2:2),scale(-2:2*4),1,0,0))
#' 
#' @seealso \code{\link{loadings.o2m}}
#' 
#' @rdname scores
#' @export
scores <- function(x, ...) UseMethod("scores")


#' @param x Object of class \code{o2m}
#' @param which_part character string. One of the following: 'Xjoint', 'Yjoint', 'Xorth' or 'Yorth'.
#' @param subset subset of scores vectors to be extracted.
#' 
#' 
#' @rdname scores
#' @export
scores.o2m <- function(x, which_part = c("Xjoint", "Yjoint", "Xorth", "Yorth"), 
                         subset = 0, ...) {
  if(any(subset != abs(round(subset)))) stop("Subset must be a vector of non-negative integers","\n")
  
  which_part = match.arg(which_part)
  which_scores = switch(which_part, Xjoint = "Tt", Yjoint = "U", Xorth = "T_Yosc", Yorth = "U_Xosc")
  scores_matrix = x[[which_scores]]
  dim_names = dimnames(scores_matrix)
  if(length(subset) == 1 && subset == 0) subset = 1:ncol(scores_matrix)
  if(max(subset) > ncol(scores_matrix)) stop("Elements in subset exceed #components","\n")
  scores_matrix = as.matrix(scores_matrix[,subset])
  dimnames(scores_matrix) <- dim_names
  
  return(scores_matrix)
}

#' Predicts X or Y
#'
#' Predicts X or Y based on new data on Y or X
#'
#' @inheritParams summary.o2m
#' @param newdata New data, which one of X or Y is specified in \code{XorY}.
#' @param XorY Character specifying whether \code{newdata} is X or Y.
#' 
#' @return Predicted Data
#' 
#' @details Prediction is done after correcting for orthogonal parts.
#' 
#' @examples
#' predict(o2m(scale(1:10), scale(1:10), 1, 0, 0), newdata = scale(1:5), XorY = "X")
#' @export
predict.o2m <- function(object, newdata, XorY = c("X","Y"), ...) {
  XorY = match.arg(XorY)
  Xnames = dimnames(newdata)
  if(!is.matrix(newdata)){
    message("newdata has class ",class(newdata),", trying to convert with as.matrix.",sep="","\n")
    newdata <- as.matrix(newdata)
    dimnames(newdata) <- Xnames
  }
  input_checker(newdata)
  switch(XorY,
         X = if(ncol(newdata) != nrow(object$W.)) stop("Number of columns mismatch!","\n"),
         Y = if(ncol(newdata) != nrow(object$C.)) stop("Number of columns mismatch!","\n"))
  
  pred = switch(XorY, 
                Y = with(object, (newdata - newdata %*% C_Xosc %*% t(P_Xosc.)) %*% C. %*% B_U %*% t(W.)), 
                X = with(object, (newdata - newdata %*% W_Yosc %*% t(P_Yosc.)) %*% W. %*% B_T. %*% t(C.)))
  
  return(pred)
}


# Penalized -------------------------------------------------

#' Check if penalization parameters satisfy input conditions
#' 
#' @param x Should be numeric matrix.
#' @param y Should be numeric matrix.
#' @param keepx Input of \code{o2m} function.
#' @param keepy Input of \code{o2m} function.
#' @param n Number of joint components.
#' 
#' @return NULL
#' @details This function throws an error if lambda is not within the range.
#' 
#' @keywords internal
#' @export
lambda_checker <- function(x,y,keepx, keepy,n) {
  if(all(is.null(keepx), is.null(keepy))) stop("Please specify 'keepx' and 'keepy'. Otherwise please set 'sparsity' to FALSE\n")
  if(any(is.null(keepx), is.null(keepy))){
    if(is.null(keepx)){
      keepx = ncol(x)
      message("'keepx' not specified, sparsity not imposed for X\n")
    }else{
      keepy = ncol(y)
      message("'keepy' not specified, sparsity not imposed for Y\n")
    }
  }
  bl_x <- !sapply(keepx, is.numeric)
  bl_y <- !sapply(keepy, is.numeric)
  if(!length(keepx) %in% c(1,n)) stop("length of 'keepx' must be equal to 1 or number of joint components","\n")
  if(!length(keepy) %in% c(1,n)) stop("length of 'keepy' must be equal to 1 or number of joint components","\n")
  if(any(c(bl_x, bl_y)))  stop("Input of keepx, keepy must be of type numeric","\n")
  if(any(c(keepx<=0, keepy<=0)))  stop("Input of keepx, keepy must be positive","\n")
  if(max(keepx) > dim(x)[2])  stop("keepx must be less then the number of column of X","\n")
  if(max(keepy) > dim(y)[2])  stop("keepx must be less then the number of column of Y","\n")
  if(length(keepx)==1){keepx <- rep(keepx,n)}
  if(length(keepy)==1){keepy <- rep(keepy,n)}
  return(list(keepx=keepx, keepy=keepy))
}

#' Check if penalization parameters for groups satisfy input conditions
#' 
#' @param groupx Vector. Input of \code{o2m} function.
#' @param groupy Vector. Input of \code{o2m} function.
#' @param keepx Input of \code{o2m} function or \code{crossval_sparsity} function.
#' @param keepy Input of \code{o2m} function or \code{crossval_sparsity} function.
#' @param n Number of joint components.
#'
#' @return NULL
#' @details This function throws an error if lambda is not within the range.
#' 
#' @keywords internal
#' @export
lambda_checker_group <- function(groupx, groupy, keepx, keepy, n) {
  if(all(is.null(keepx), is.null(keepy))) stop("Please specify 'keepx' and 'keepy'. Otherwise please set 'sparsity' to FALSE and run O2PLS","\n")
  if(any(is.null(keepx), is.null(keepy))){
    if(is.null(keepx)){
      keepx = length(unique(groupx))
      message("'keepx' not specified, sparsity not imposed in X","\n")
    }else{
      keepy = length(unique(groupy))
      message("'keepy' not specified, sparsity not imposed in Y","\n")
    }
  }
  bl_x <- !sapply(keepx, is.numeric)
  bl_y <- !sapply(keepy, is.numeric)
  if(any(c(bl_x, bl_y)))  stop("Input of keepx, keepy must be positive numbers","\n")
  if(any(c(keepx<=0, keepy<=0)))  stop("Input of keepx, keepy must be positive","\n")
  if(max(keepx) > length(unique(groupx)))  stop("keepx must not exceed the number of groups in X","\n")
  if(max(keepy) > length(unique(groupy)))  stop("keepy must not exceed the number of groups in Y","\n")
  if(length(keepx)==1){keepx <- rep(keepx,n)}
  if(length(keepy)==1){keepy <- rep(keepy,n)}
  return(list(keepx=keepx, keepy=keepy))
}

#' Check if sparsity parameters satisfy input conditions in cross-validation functions
#' 
#' @param x Should be numeric matrix.
#' @param y Should be numeric matrix.
#' @param keepx_seq Input of \code{crossval_sparsity} function.
#' @param keepy_seq Input of \code{crossval_sparsity} function.
#' 
#' @return NULL
#' @details This function throws an error if sparsity parameters are not valid.
#' 
#' @keywords internal
#' @export
cv_lambda_checker <- function(x,y,keepx_seq=NULL, keepy_seq=NULL) {
  if(all(is.null(keepx_seq), is.null(keepy_seq))) stop("Please specify at least one of 'keepx_seq' and 'keepy_seq'","\n")
  bl_x <- !sapply(keepx_seq, is.numeric)
  bl_y <- !sapply(keepy_seq, is.numeric)
  if(any(c(bl_x, bl_y)))  stop("Input of keepx_seq, keepy_seq must all be positive numbers","\n")
  if(any(c(keepx_seq<=0, keepy_seq<=0)))  stop("Input of keepx_seq, keepy_seq must all be positive","\n")
  if(max(keepx_seq) > dim(x)[2])  stop("all numbers in keepx_seq must be less then the number of column of X","\n")
  if(max(keepy_seq) > dim(y)[2])  stop("all numbers in keepy_seq must be less then the number of column of Y","\n")
}

#' Check if sparsity parameters satisfy input conditions in cross-validation functions
#' 
#' @param groupx Vector. Input of \code{crossval_sparsity} function.
#' @param groupy Vector. Input of \code{crossval_sparsity} function.
#' @param keepx_seq Input of \code{crossval_sparsity} function.
#' @param keepy_seq Input of \code{crossval_sparsity} function.
#' 
#' @return NULL
#' @details This function throws an error if sparsity parameters are not valid.
#' 
#' @keywords internal
#' @export
cv_lambda_checker_group <- function(groupx,groupy,keepx_seq=NULL, keepy_seq=NULL) {
  if(any(is.null(keepx_seq), is.null(keepy_seq))) stop("Please specify 'keepx_seq' and 'keepy_seq'","\n")
  bl_x <- !sapply(keepx_seq, is.numeric)
  bl_y <- !sapply(keepy_seq, is.numeric)
  if(any(c(bl_x, bl_y)))  stop("Input of keepx_seq, keepy_seq must all be positive numbers","\n")
  if(any(c(keepx_seq<=0, keepy_seq<=0)))  stop("Input of keepx_seq, keepy_seq must all be positive","\n")
  if(max(keepx_seq) > length(unique(groupx)))  stop("all numbers in keepx_seq must be less then the number of groups in X","\n")
  if(max(keepy_seq) > length(unique(groupy)))  stop("all numbers in keepy_seq must be less then the number of groups in Y","\n")
}

#' Soft threshholding a vector with respect to a number of variables
#'
#' @param x Numerical vector
#' @param keepx How many non-zero
#' @return Soft-thresholded vector
#'
#' @export
#' 
thresh_n <- function (x, keepx){
  nx <- length(x) - keepx
  if (nx != 0) {
    absa = abs(x)
    if (any(rank(absa, ties.method = "max") <= nx)) {
      x = ifelse(rank(absa, ties.method = "max") <= nx, 
                 0, sign(x) * (absa - max(absa[rank(absa, ties.method = "max") <= 
                                                 nx])))
    }
  }
  x
}

#' Soft threshholding a vector with respect to a number of groups
#'
#' @param w Numerical loading vector 
#' @param keep_gr How many groups to retain
#' @param index_gr List of index and size. index are the index of variables belongs to the group in the original vector, size is the group size
#' @return A list containing sparse loading vector and names of the selected groups
#'
#' @export
#' 
thresh_n_gr <- function (w, keep_gr, index_gr){
  nr <- length(index_gr)
  # each group l2 norm
  gr_norm <- sapply(1:nr, function(j){
    wj <- w[index_gr[[j]]$index] 
    normj <- norm_vec(wj)
    return(normj)
  })
  
  # find weights for each group
  # first sorted critical lambda value for each group
  size <- sapply(index_gr, function(e) e$size)
  lambda_seq <- (gr_norm/sqrt(size)) %>% sort(decreasing = T)
  if(keep_gr == nr){
    lambda <- 0
  }else{
    lambda <- lambda_seq[keep_gr+1]
  }
  coef_seq <- (gr_norm - sqrt(size)*lambda)/gr_norm
  coef_seq[coef_seq < 0] <- 0
  
  for(j in 1:nr){
    w[index_gr[[j]]$index] <- coef_seq[j] * w[index_gr[[j]]$index]
  }
  
  select_gr <- names(index_gr)[which(coef_seq > 0)]
  return(list(w = w, select_gr = select_gr))
}


#' Norm of a vector
#'
#' @param x Numerical vector
#' @return L2 norm of a vector
#' @export
norm_vec <- function(x) sqrt(
  sum(x^2))

#' Orthogonalize a sparse loading vector with regard to a matrix
#'
#' @param x sparse loading vector to be orthogonalized
#' @param W sparse loading matrix of the previous loading vectors
#' @return A sparse loading vector
#' @export
orth_vec <- function(x, W){
  # get non-zero positions in x and subset W
  W %<>% as.matrix
  l <- length(x)
  pos <- which(x!=0)
  x <- x[pos]
  W <- W[pos, ,drop = F]
  # First check if W is already 0
  # Step1: delete row i if ith row in W contain all 0
  # Step2: delete column j in W if it contains all 0
  # this makes sure W'W is invertible (expect when W ncol > nrow)
  if(all(W == 0)){
    x_orth <- x
    # print('all 0')
  }else{
    # Step1
    #print(W)
    rowi <- sapply(1:nrow(W), function(i) any(W[i, ]!=0)) %>% which
    W <- W[rowi, ,drop = F]
    x_sub <- x[rowi]  # x_sub = x_original[pos[rowi]]
    # Step2
    colj <- sapply(1:ncol(W), function(j) any(W[ ,j]!=0)) %>% which
    W <- W[,colj,drop = F]
    
    # when W ncol > nrow, W is singular
    if(ncol(W)>nrow(W)){
      x_sub <- rep(0, length(x_sub))
      #print("ncol > nrow")
    }else{
      # Orthogonal projection
      n <- nrow(W)
      I <- diag(1, n)
      #print(W)
      x_sub <- (I - W %*% solve(t(W)%*%W) %*% t(W)) %*% x_sub
    }
    
    x[rowi] <- x_sub
    x_orth <- x
  }
  x_orth <- x_orth/norm_vec(x_orth)
  final <- rep(0, l)
  final[pos] <- x_orth
  return(final)
}

