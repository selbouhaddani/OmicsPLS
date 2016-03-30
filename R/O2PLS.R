#' O2PLS: Two-Way Orthogonal Partial Least Squares
#'
#' This is based on work of (Trygg, Wold, 2003).
#' Includes the O2PLS fit, some misc functions and some cross-validation tools.
#' 
#' @author
#' Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}),
#' Jeanine Houwing-Duistermaat (\email{J.J.Houwing@@lumc.nl}),
#' Geurt Jongbloed (\email{G.Jongbloed@@tudelft.nl}),
#' Szymon Kielbasa (\email{S.M.Kielbasa@@lumc.nl}),
#' Hae-Won Uh (\email{H.Uh@@lumc.nl}).
#'
#' Maintainer: Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}).
#' 
#' @section Model and assumptions:
#' \strong{Note that the rows of \code{X} and \code{Y} are the subjects and columns are variables.}
#' The number of columns may be different, but the subjects should be the same in both datasets.
#' 
#' The O2PLS model (Trygg \& Wold, 2003) decomposes two datasets \eqn{X} and \eqn{Y} into three parts. 
#' \itemize{
#'  \item{1.} A joint part, representing the relationship between \eqn{X} and \eqn{Y}
#'  \item{2.} An orthogonal part, representing the unrelated latent variation in \eqn{X} and \eqn{Y} separately.
#'  \item{3.} A noise part capturing all residual variation.
#' }
#' 
#' See also the corresponding paper for interpretation (el Bouhaddani et al, 2016).
#' 
#' 
#' @section Fitting:
#' The O2PLS fit is done with \code{\link{o2m}}. 
#' For data \code{X} and \code{Y} you can run \code{o2m(X,Y,n,nx,ny)} for an O2PLS fit with \code{n} joint and \code{nx, ny} orthogonal components.
#' See the help page of \code{\link{o2m}} for more information on parameters.
#' There are four ways to obtain an O2PLS fit, depending on the dimensionality.
#' \itemize{
#'  \item{} For the not-too-high dimensional case, you may use \code{\link{o2m}} with default parameters. E.g. \code{o2m(X,Y,n,nx,ny)}.
#'  \item{} In case you don't want the fancy output, but only the parameters, you may add \code{stripped = TRUE} to obtain a stripped version of \code{o2m} which avoids calculating and storing some matrices. E.g. \code{o2m(X,Y,n,nx,ny,stripped=TRUE)}.
#'  \item{} For high dimensional cases defined by \code{ncol(X)>p_thresh} and \code{ncol(Y)>q_thresh} a Power-Method approach is used which avoids storing large matrices. E.g. \code{o2m(X,Y,n,nx,ny,p_thresh=3000,q_thresh=3000)}.
#'  The thresholds are by default both at 3000 variables.
#'  \item{} If you want a stripped version in the high dimensional case, add \code{stripped = TRUE}E.g. \code{o2m(X,Y,n,nx,ny,stripped=TRUE,p_thresh=3000,q_thresh=3000)}.
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
#' See \code{citation("O2PLS")} for our proposed approach for determining the number of components.
#' \itemize{
#'  \item{} Cross-validation (CV) is done with \code{\link{loocv}} or \code{\link{adjR2}}, the last has built in parallelization (when you use Windows!) which relies on the \code{parallel} package.
#'  This part is not yet fully developed, and a somewhat advanced understanding of R is needed to fully exploit CV here.
#'  However a slow \code{for} loop implementation of the cross-validated prediction error is implemented in \code{loocv(X,Y,a,a2,b2,kcv)}, where \code{X, Y} are the data and \code{a,a2,b2} are integer vectors representing \code{n,nx,ny}.
#'  \code{kcv} is the number of folds, with \code{kcv = nrow(X)} for Leave-One-Out CV. 
#'  \code{\link{loocv}} was proposed to find the number of joint components.
#'  \item{} For \code{\link{adjR2}} the same holds, but this function implements the coefficient of determination of the inner relation \eqn{U = TB + H} rather than the prediction error.
#'  This was proposed to find the number of orthogonal components.
#' }
#' 
#' @section Misc:
#' Also some handy tools are available
#' \itemize{
#'  \item{} \code{\link{orth}(X)} is a function to obtain an orthogonalized version of a matrix or vector \code{X}.
#'  \item{} \code{\link{ssq}(X)} is a function to calculate the sum of squares (or squared Frobenius norm) of \code{X}. See also \code{\link{vnorm}} for calculating the norm of each column in \code{X}.
#'  \item{} \code{\link{mse}(x, y)} returns the mean squared difference between two matrices/vectors. By default \code{y=0}.
#' }
#' 
#' @section Citation:
#' If you use the R package in your research, please cite the corresponding paper:
#' 
#' \strong{Bouhaddani, S., Houwing-duistermaat, J., Jongbloed, G., Salo, P., Perola, M., & Uh, H.-W.} (2016).
#' \emph{Evaluation of O2PLS in Omics data integration.}
#' BMC Bioinformatics BMTL Supplement. doi:10.1186/s12859-015-0854-z
#' 
#' The bibtex entry can be obtained with command \code{citation("O2PLS")}.
#' Thank You!
#' 
#' The original paper proposing O2PLS is
#' 
#' \strong{Trygg, J., & Wold, S.} (2003). 
#' \emph{O2-PLS, a two-block (X-Y) latent variable regression (LVR) method with an integral OSC filter.} 
#' Journal of Chemometrics, 17(1), 53-64. http://doi.org/10.1002/cem.775
#' 
#' @docType package
#' @name O2PLS
#' @keywords O2PLS
#' @import parallel ggplot2
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
  if(!is.numeric(X)) stop("Input is not numeric, but of mode ",mode(X))
  if(!is.matrix(X)) stop("Input is not a matrix, but of class ",class(X))
  if(any(is.na(X))) stop("Input contains NA's or NaN's")
  if(!any(is.finite(X))) stop("Input contains non-finite elements")
  
  if (!is.null(Y)) {
    if(!is.numeric(Y)) stop("Input is not numeric, but of mode ",mode(Y))
    if(!is.matrix(Y)) stop("Input is not a matrix, but of class ",class(Y))
    if(any(is.na(Y))) stop("Input contains NA's or NaN's")
    if(!any(is.finite(Y))) stop("Input contains non-finite elements")
    if(nrow(X) != nrow(Y)) stop("# rows don't match: ",nrow(X)," versus ",nrow(Y))
    if(!identical(rownames(X), rownames(Y))) warning("Caution: Rownames don't match!")
  }
  NULL
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
    if(ncol(X) != ncol(X_true)) stop("# columns don't match:",ncol(X),"versus",ncol(X_true))
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
  if(any(diag(sign_e)==0)){warning("Orthogonalization made some columns orthogonal to original columns")}
  
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
  if(length(x) != length(y) && length(y) != 0) message("Comparing lengths ",length(x)," with ",length(y))
  mean((x - y)^2, na.rm = na.rm)
}



#' Root MSE of Prediction
#'
#' Calculates the Root MSE of prediction on test data. Only tested to work inside \code{\link{loocv}}.
#'
#' @param Xtst Numeric vector or matrix.
#' @param Ytst Numeric vector or matrix.
#' @param fit \code{\link{o2m}} fit (on data without \code{Xtst} and \code{Ytst}).
#' @param combi Logical. Should the symmetrized MSE be used, i.e. 
#' @details This function is the building block for \code{\link{loocv}}, as it produced the prediction error for test (left out) data.
#' @return Mean squares difference between predicted Y and true Y
#' @export
rmsep <- function(Xtst, Ytst, fit, combi = FALSE) {
  if(!inherits(fit,"o2m")) stop("fit should be an O2PLS fit")
  
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
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param a Vector of integers. Contains the numbers of joint components.
#' @param a2 Vector of integers. Contains the numbers of orthogonal components in \eqn{X}.
#' @param b2 Vector of integers. Contains the numbers of orthogonal components in \eqn{Y}.
#' @param fitted_model List. O2PLS model fit with \code{\link{o2m}}. Is used to calculate the apparent error without recalculating this fit.
#' @param func Function to fit the O2PLS model with. Only \code{\link{o2m}} and \code{\link{o2m_stripped}} are supported.
#' @param app_err Logical. Should the apparent error also be computed?
#' @param kcv Integer. The value of \eqn{k}, i.e. the number of folds. Choose \eqn{N} for LOO-CV.
#' @details Note that this function can be easily parallelized (on Windows e.g. with the \code{parallel} package.).
#' @return List with two numeric vectors:
#' \item{CVerr}{Contains the k-fold CV estimated RMSEP}
#' \item{Fiterr}{Contains the apparent error}
#' @details The parameters \code{a}, \code{a2} and \code{b2} can be integers or vectors of integers. A for loop is used to loop over all combinations.
#' The resulting output is a list, which is more easy to interpret if you use \code{array(unlist(output_of_loocv$CVerr))} as in the example below.
#' The array wil have varying \code{a} along the first dimension and \code{a2} and \code{b2} along the second and third respectively.
#' Typing \code{example(loocv)} (hopefully) clarifies this function.
#' @examples
#' result=loocv(matrix(rnorm(10*100),ncol=10),matrix(rnorm(10*100),ncol=10),a=1:3,a2=0:1,b2=0:1,func=o2m_stripped,kcv=2)
#' names_for_a=sapply(1:3,function(i){paste('a',i,sep='=')})
#' names_for_a2=sapply(0:1,function(i){paste('a2',i,sep='=')})
#' names_for_b2=sapply(0:1,function(i){paste('b2',i,sep='=')})
#' array(unlist(result$CVerr),dim=c(3,2,2),dimnames=list(names_for_a,names_for_a2,names_for_b2))
#' @export
loocv <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, fitted_model = NULL, func = o2m_stripped, app_err = F, 
                  kcv)
{
  stopifnot(all(a == round(a)), all(a2 == round(a2)), all(b2 == round(b2)))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  if (!is.null(fitted_model)) {
    app_err <- F
    warning("apparent error calculated with provided fit")
  }
  # determine type of model
  type <- 3  #ifelse(deparse(substitute(func))=='o2m',3,ifelse(deparse(substitute(func))=='oplsm',2,1))
  
  N <- nrow(X)
  if (N != nrow(Y)) {
    stop("N not the same")
  }
  mean_err <- mean_fit <- NA * 1:max(length(a), length(a2), length(b2))
  k <- 0
  
  # blocks contains the begin and endpoints of test indices to use
  blocks <- c(seq(0, N, by = floor(N/kcv)), N)
  
  # loop through chosen parameters
  for (j in a) {
    for (j2 in a2) {
      for (j3 in b2) {
        k <- k + 1
        err <- NA * 1:kcv
        folds <- sample(N)
        # loop through number of folds
        for (i in 1:kcv) {
          ii <- (blocks[i] + 1):(blocks[i + 1])
          if (type == 3) {
            pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = j, nx = j2, ny = j3)
          }
          # if(type==2){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j,n_orth=j2)}
          # if(type==1){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j)}
          fit <- try(do.call(func, pars), silent = T)
          err[i] <- ifelse(class(fit) == "try-error", NA, rmsep(X[folds[ii], ], Y[folds[ii], ], 
                                                                fit))
        }
        mean_err[k] <- mean(err)
        # calculate apparent error
        if (app_err && is.null(fitted_model)) {
          if (type == 3) {
            pars2 <- list(X = X, Y = Y, n = j, nx = j2, ny = j3)
          }
          # if(class(fit)=='oplsm'){pars2=list(X=X,Y=Y,ncomp=j,n_orth=j2)}
          # if(class(fit)=='plsm'){pars2=list(X=X,Y=Y,ncomp=j)}
          fit2 <- try(do.call(func, pars2), F)
          mean_fit[k] <- ifelse(class(fit) == "try-error", NA, rmsep(X, Y, fit2))
          # print('1e loop')
        }
        if (!is.null(fitted_model)) {
          mean_fit[k] <- rmsep(X, Y, fitted_model)
          # print('2e loop')
        }
      }
    }
  }
  return(list(CVerr = mean_err, Fiterr = mean_fit))
}

#' Gridwise adjusted R2 for O2PLS
#'
#' For (a grid of) values for \code{a}, \code{nx} and \code{ny}, \code{loocv} calculates the R2 of the joint part. Parallel computing is supported on Windows with package \code{parallel}.
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
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
#' @examples
#' result=adjR2(matrix(rnorm(10*100),ncol=10),matrix(rnorm(10*100),ncol=10),a=1:3,a2=0:1,b2=0:1,func=o2m_stripped)
#' names_for_a=sapply(1:3,function(i){paste('a',i,sep='=')})
#' names_for_a2=sapply(0:1,function(i){paste('a2',i,sep='=')})
#' names_for_b2=sapply(0:1,function(i){paste('b2',i,sep='=')})
#' array(unlist(result[1,]),dim=c(3,2,2),dimnames=list(names_for_a,names_for_a2,names_for_b2))
#' @export
adjR2 <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, func = o2m_stripped, parall = F, cl = NULL)
{
  stopifnot(all(a == round(a)), all(a2 == round(a2)), all(b2 == round(b2)))
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
    stopifnot("cluster" %in% class(cl))
    S_apply <- parSapply
  }
  
  pars1 <- merge(merge(data.frame(a = a), data.frame(a2 = a2)), data.frame(b2 = b2))
  pars2 <- apply(pars1, 1, as.list)
  
  # cl <- makeCluster(rep( 'localhost', detectCores()),type='SOCK') clusterExport(cl=cl,
  # varlist=c('X','Y','N','pars2','ssq','o2m'))
  outp <- S_apply(cl, pars2, function(p) {
    fit <- func(X, Y, p$a, p$a2, p$b2)
    RX <- 1 - ssq(fit$H_UT)/ssq(fit$U)
    RY <- 1 - ssq(fit$H_TU)/ssq(fit$Tt)
    adjRX <- RX  #1 - (1 - RX)*(N - 1)/(N - p$a - 1)
    adjRY <- RY  #1 - (1 - RY)*(N - 1)/(N - p$a - 1)
    return(c(adjR2X = adjRX, adjR2Y = adjRY))
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
  if(!inherits(fit,"o2m")) stop("fit should be an O2PLS fit")
  
  if (!is.matrix(Xtst)) Xtst <- t(Xtst)
  
  if (!is.matrix(Ytst)) Ytst <- t(Ytst)
  
  input_checker(Xtst, Ytst)

  Yhat <- Xtst %*% fit$W. %*% fit$B_T %*% t(fit$C.)
  Xhat <- Ytst %*% fit$C. %*% fit$B_U %*% t(fit$W.)
  
  return(sqrt(mse(Yhat, Ytst)) + sqrt(mse(Xhat, Xtst)))
}

#' K-fold CV based on symmetrized prediction error
#'
#' The prediction error of both \code{X~Xhat} and \code{Y~Yhat} are summed. This provides a symmetrized version of \code{\link{loocv}}.
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param a Vector or integers. Contains the #joint components.
#' @param a2 Vector or integers. Contains the number of orthogonal components in \eqn{X}.
#' @param b2 Vector or integers. Contains the number of orthogonal components in \eqn{Y}.
#' @param fitted_model List. O2PLS model fit with \code{\link{o2m}}. Is used to calculate the apparent error without recalculating this fit.
#' @param func Function to fit the O2PLS model with. Only \code{\link{o2m}} and \code{\link{o2m_stripped}} are supported.
#' @param app_err Logical. Should the apparent error also be computed?
#' @param kcv Integer. The value of \eqn{k}, i.e. the number of folds. Choose \eqn{N} for LOO-CV.
#' @details Note that this function can be easily parallelized (on Windows e.g. with the \code{parallel} package.).
#' @return List with two numeric vectors:
#' \item{CVerr}{Contains the k-fold CV estimated RMSEP}
#' \item{Fiterr}{Contains the apparent error}
#'
#' @examples
#' loocv_combi(matrix(c(-2:2)),matrix(c(-2:2*4)),1,0,0,func=o2m,kcv=5)
#' @export
loocv_combi <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, fitted_model = NULL, func = o2m_stripped, app_err = F, 
                        kcv)
{
  stopifnot(all(a == round(a)), all(a2 == round(a2)), all(b2 == round(b2)), kcv == round(kcv[1]))
  stopifnot(is.logical(app_err), is.function(func))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  if (!is.null(fitted_model)) {
    if(class(fitted_model) != 'o2m'){stop("fitted_model should be of class 'o2m' or NULL")}
    app_err <- F
    warning("apparent error calculated with provided fit")
  }
  
  # determine type of model
  type <- 3  #ifelse(deparse(substitute(func))=='o2m',3,ifelse(deparse(substitute(func))=='oplsm',2,1))
  
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  mean_err <- mean_fit <- NA * 1:max(length(a), length(a2), length(b2))
  k <- 0
  
  # blocks contains the begin and endpoints of test indices to use
  blocks <- c(seq(0, N, by = floor(N/kcv)), N)
  
  # loop through chosen parameters
  for (j in a) {
    for (j2 in a2) {
      for (j3 in b2) {
        k <- k + 1
        err <- NA * 1:kcv
        folds <- sample(N)
        # loop through number of folds
        for (i in 1:kcv) {
          ii <- (blocks[i] + 1):(blocks[i + 1])
          if (type == 3) {
            pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = j, nx = j2, ny = j3)
          }
          if (type == 2) {
            pars <- list(X = X[-i, ], Y = Y[-i, ], ncomp = j, n_orth = j2)
          }
          if (type == 1) {
            pars <- list(X = X[-i, ], Y = Y[-i, ], ncomp = j)
          }
          fit <- try(do.call(func, pars), silent = T)
          err[i] <- ifelse(class(fit) == "try-error", NA, rmsep_combi(X[folds[ii], ], Y[folds[ii], 
                                                                                        ], fit))
        }
        mean_err[k] <- mean(err)
        # calculate apparent error
        if (app_err && is.null(fitted_model)) {
          if (class(fit) == "o2m") {
            pars2 <- list(X = X, Y = Y, n = j, nx = j2, ny = j3)
          }
          if (class(fit) == "oplsm") {
            pars2 <- list(X = X, Y = Y, ncomp = j, n_orth = j2)
          }
          if (class(fit) == "plsm") {
            pars2 <- list(X = X, Y = Y, ncomp = j)
          }
          fit2 <- try(do.call(func, pars2), T)
          mean_fit[k] <- ifelse(class(fit) == "try-error", NA, rmsep_combi(X, Y, fit2))
          print("1e loop")
        }
        if (!is.null(fitted_model)) {
          mean_fit[k] <- rmsep_combi(X, Y, fitted_model)
          print("2e loop")
        }
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
  if(x$flags$stripped) cat("O2PLS fit: Stripped \n") 
    else if(x$flags$highd) cat("O2PLS fit: High dimensional \n") 
      else cat("O2PLS fit \n")
  cat("with ",n," joint components  \n",sep='')
  cat("and  ",nx," orthogonal components in X \n",sep='')
  cat("and  ",ny," orthogonal components in Y \n",sep='')
  cat("Elapsed time: ",x$flags$time, " sec", sep='')
}


#' Plot one or two loading vectors for class o2m
#' 
#' This function plots one or two loading vectors, by default with ggplot2. 
#' 
#' @param x An O2PLS fit, with class 'o2m'
#' @param loading_name character string. One of the following: 'Xjoint', 'Yjoint', 'Xorth' or 'Yorth'.
#' @param i Integer. First component to be plotted.
#' @param j NULL (default) or Integer. Second component to be plotted.
#' @param use_ggplot2 Logical. Default is \code{TRUE}. If \code{FALSE}, the usual plot device will be used.
#' @param label Character, either 'number' or 'colnames'. The first option prints numbers, the second prints the colnames
#' @param ... Further arguments to \code{geom_text}, such as size, col, alpha, etc.
#' 
#' @return If \code{use_ggplot2} is \code{TRUE} a ggplot2 object. Else NULL.
#' 
#' @export
plot.o2m <- function (x, loading_name = c("Xjoint", "Yjoint", "Xorth", "Yorth"), i = 1, j = NULL, use_ggplot2=TRUE, label = c("number", "colnames"), ...)
{
  stopifnot(i == round(i), is.logical(use_ggplot2))
  
  fit <- list()
  loading_name = match.arg(loading_name)
  which_load = switch(loading_name, Xjoint = "W.", Yjoint = "C.", Xorth = "P_Yosc.", Yorth = "P_Xosc.")
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
  
  label = match.arg(label)
  if(label == "colnames" && !is.null(rownames(x[which_load][[1]]))) {
    label = rownames(x[which_load][[1]])
    } else label = 1:p
  if(label == "colnames" && is.null(rownames(x[which_load][[1]]))) message("No labels found in colnames, proceeding...")
  
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
#' @export
summary.o2m <- function(object, digits = 3, ...) {
  fit <- object
  a <- ncol(fit$W.)
  if(digits != round(digits) || digits <= 0) stop("digits must be a positive integer")
  outp <- with( fit, list(
    Comp = a,
    R2_X = R2X,
    R2_Y = R2Y,
    R2_Xjoint = R2Xcorr,
    R2_Yjoint = R2Ycorr,
    R2_Xhat = R2Xhat,
    R2_Yhat = R2Yhat,
    B_T = B_T.,
    B_U = B_U,
    flags = flags
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
print.summary.o2m <- function(x, digits = 3, ...){
  cat("\n*** Summary of the O2PLS fit *** \n\n")
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
    cat("-- Predictable variation in Y by X:\n")
    cat("Variation in Yhat:",round(R2_Yhat,digits),"\n")
    cat("-- Predictable variation in X by Y:\n")
    cat("Variation in Xhat:",round(R2_Xhat,digits),"\n")
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
#' @rdname loadings
#' @export
loadings <- function(x, ...) UseMethod("loadings")


#' @param loading_name character string. One of the following: 'Xjoint', 'Yjoint', 'Xorth' or 'Yorth'.
#' @param subset subset of loading vectors to be extracted.
#' @param sorted Logical. Should the rows of the loadings be sorted according to the 
#' absolute magnitute of the first column?
#' 
#' @rdname loadings
#' @export
loadings.o2m <- function(x, loading_name = c("Xjoint", "Yjoint", "Xorth", "Yorth"), 
                         subset = 0, sorted = FALSE, ...) {
  if(any(subset != abs(round(subset)))) stop("subset must be a vector of non-negative integers")
  
  loading_name = match.arg(loading_name)
  which_load = switch(loading_name, Xjoint = "W.", Yjoint = "C.", Xorth = "P_Yosc.", Yorth = "P_Xosc.")
  loading_matrix = x[[which_load]]
  dim_names = dimnames(loading_matrix)
  if(length(subset) == 1 && subset == 0) subset = 1:ncol(loading_matrix)
  if(max(subset) > ncol(loading_matrix)) stop("Elements in subset exceed #components")
  loading_matrix = as.matrix(loading_matrix[,subset])
  dimnames(loading_matrix) <- dim_names
  
  if(sorted){
    order_load = order(loading_matrix[,1])
    if(is.null(dim_names[[1]])) dim_names[[1]] <- order_load
    loading_matrix = loading_matrix[order_load,]
  }
  return(loading_matrix)
}

#' Predicts X or Y
#'
#' Predicts X or Y based on new data on Y or X
#'
#' @inheritParams summary.o2m
#' @param newdata New data, which one of X or Y is specified in \code{XorY}.
#' @param XorY Character specifying \code{newdata} is X or Y.
#' 
#' @return Predicted Data
#' @examples
#' predict(o2m(scale(1:10), scale(1:10), 1, 0, 0), newdata = scale(1:5), XorY = "X")
#' @export
predict.o2m <- function(object, newdata, XorY = c("X","Y"), ...) {
  XorY = match.arg(XorY)
  switch(XorY,
         X = if(ncol(newdata) != nrow(object$W.)) stop("Number of columns mismatch!"),
         Y = if(ncol(newdata) != nrow(object$C.)) stop("Number of columns mismatch!"))
  
  pred = switch(XorY, 
                Y = with(object,newdata %*% C. %*% B_U %*% t(W.)), 
                X = with(object,newdata %*% W. %*% B_T. %*% t(C.)))
  
  return(pred)
}
