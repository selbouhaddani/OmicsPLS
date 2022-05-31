#' Perform O2PLS data integration with two-way orthogonal corrections
#'
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALE THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#'
#' @param X Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param Y Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param n Integer. Number of joint PLS components. Must be positive.
#' @param nx Integer. Number of orthogonal components in \eqn{X}. Negative values are interpreted as 0
#' @param ny Integer. Number of orthogonal components in \eqn{Y}. Negative values are interpreted as 0
#' @param stripped Logical. Use the stripped version of o2m (usually when cross-validating)?
#' @param p_thresh Integer. If \code{X} has more than \code{p_thresh} columns, a power method optimization is used, see \code{\link{o2m2}}
#' @param q_thresh Integer. If \code{Y} has more than \code{q_thresh} columns, a power method optimization is used, see \code{\link{o2m2}}
#' @param tol Double. Threshold for which the NIPALS method is deemed converged. Must be positive.
#' @param max_iterations Integer. Maximum number of iterations for the NIPALS method. 
#' @param sparse Boolean. Default value is \code{FALSE}, in which case O2PLS will be fitted. Set to \code{TRUE} for GO2PLS.
#' @param groupx Vector. Used when \code{sparse = TRUE}. A vector of strings indicating group names of each X-variable. Its length must be equal to the number of variables in \eqn{X}. The order of group names must corresponds to the order of the variables. 
#' @param groupy Vector. Used when \code{sparse = TRUE}. A vector of strings indicating group names of each Y-variable. The length must be equal to the number of variables in \eqn{Y}. The order of group names must corresponds to the order of the variables.
#' @param keepx Vector. Used when \code{sparse = TRUE}. A vector of length \code{n} indicating how many variables (or groups if \code{groupx} is provided) to keep in each of the joint component of \eqn{X}. If the input is an integer, all the components will have the same amount of variables or groups retained.
#' @param keepy Vector. Used when \code{sparse = TRUE}. A vector of length \code{n} indicating how many variables (or groups if \code{groupx} is provided) to keep in each of the joint component of \eqn{Y}. If the input is an integer, all the components will have the same amount of variables or groups retained.
#' @param max_iterations_sparsity Integer. Used when \code{sparse = TRUE}. Maximum number of iterations for the NIPALS method for GO2PLS.
#'
#' @return A list containing
#'    \item{Tt}{Joint \eqn{X} scores}
#'    \item{W.}{Joint \eqn{X} loadings}
#'    \item{U}{Joint \eqn{Y} scores}
#'    \item{C.}{Joint \eqn{Y} loadings}
#'    \item{E}{Residuals in \eqn{X}}
#'    \item{Ff}{Residuals in \eqn{Y}}
#'    \item{T_Yosc}{Orthogonal \eqn{X} scores}
#'    \item{P_Yosc.}{Orthogonal \eqn{X} loadings}
#'    \item{W_Yosc}{Orthogonal \eqn{X} weights}
#'    \item{U_Xosc}{Orthogonal \eqn{Y} scores}
#'    \item{P_Xosc.}{Orthogonal \eqn{Y} loadings}
#'    \item{C_Xosc}{Orthogonal \eqn{Y} weights}
#'    \item{B_U}{Regression coefficient in \code{Tt} ~ \code{U}}
#'    \item{B_T.}{Regression coefficient in \code{U} ~ \code{Tt}}
#'    \item{H_TU}{Residuals in \code{Tt} in \code{Tt} ~ \code{U}}
#'    \item{H_UT}{Residuals in \code{U} in \code{U} ~ \code{Tt}}
#'    \item{X_hat}{Prediction of \eqn{X} with \eqn{Y}}
#'    \item{Y_hat}{Prediction of \eqn{Y} with \eqn{X}}
#'    \item{R2X}{Variation (measured with \code{\link{ssq}}) of the modeled part in \eqn{X} (defined by joint + orthogonal variation) as proportion of variation in \eqn{X}}
#'    \item{R2Y}{Variation (measured with \code{\link{ssq}}) of the modeled part in \eqn{Y} (defined by joint + orthogonal variation) as proportion of variation in \eqn{Y}}
#'    \item{R2Xcorr}{Variation (measured with \code{\link{ssq}}) of the joint part in \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Ycorr}{Variation (measured with \code{\link{ssq}}) of the joint part in \eqn{Y} as proportion of variation in \eqn{Y}}
#'    \item{R2X_YO}{Variation (measured with \code{\link{ssq}}) of the orthogonal part in \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Y_XO}{Variation (measured with \code{\link{ssq}}) of the orthogonal part in \eqn{Y} as proportion of variation in \eqn{Y}}
#'    \item{R2Xhat}{Variation (measured with \code{\link{ssq}}) of the predicted \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Yhat}{Variation (measured with \code{\link{ssq}}) of the predicted \eqn{Y} as proportion of variation in \eqn{Y}}
#'    \item{W_gr}{Joint loadings of \eqn{X} at group level (only available when GO2PLS is used)}
#'    \item{C_gr}{Joint loadings of \eqn{Y} at group level (only available when GO2PLS is used)}
#'
#' @details If both \code{nx} and \code{ny} are zero, \code{o2m} is equivalent to PLS2 with orthonormal loadings.
#' This is a `slower' (in terms of memory) implementation of O2PLS, and is using \code{\link{svd}}, use \code{stripped=T} for a stripped version with less output.
#' If either \code{ncol(X) > p_thresh} or \code{ncol(Y) > q_thresh}, the NIPALS method is used which does not store the entire covariance matrix.
#' The squared error between iterands in the NIPALS approach can be adjusted with \code{tol}.
#' The maximum number of iterations in the NIPALS approach is tuned by \code{max_iterations}.
#'
#' @examples
#' test_X <- scale(matrix(rnorm(100*10),100,10))
#' test_Y <- scale(matrix(rnorm(100*11),100,11))
#' #  --------- Default run ------------ 
#' o2m(test_X, test_Y, 3, 2, 1)
#' #  ---------- Stripped version ------------- 
#' o2m(test_X, test_Y, 3, 2, 1, stripped = TRUE)
#' #  ---------- High dimensional version ---------- 
#' o2m(test_X, test_Y, 3, 2, 1, p_thresh = 1)
#' #  ------ High D and stripped version --------- 
#' o2m(test_X, test_Y, 3, 2, 1, stripped = TRUE, p_thresh = 1)
#' #  ------ Now with more iterations -------- 
#' o2m(test_X, test_Y, 3, 2, 1, stripped = TRUE, p_thresh = 1, max_iterations = 1e6)
#' #  ---------------------------------- 
#'
#' @seealso \code{\link{summary.o2m}}, \code{\link{plot.o2m}}, \code{\link{crossval_o2m_adjR2}}, \code{\link{crossval_sparsity}}
#'
#' @export
o2m <- function(X, Y, n, nx, ny, stripped = FALSE, 
                p_thresh = 3000, q_thresh = p_thresh, tol = 1e-10, max_iterations = 1000,
                sparse = F, groupx = NULL, groupy = NULL, keepx = NULL, keepy = NULL, 
                max_iterations_sparsity = 1000) {
  tic <- proc.time()
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  if(!is.matrix(X)){
    message("X has class ",class(X),", trying to convert with as.matrix.\n",sep="")
    X <- as.matrix(X)
    dimnames(X) <- Xnames
  }
  if(!is.matrix(Y)){
    message("Y has class ",class(Y),", trying to convert with as.matrix.\n",sep="")
    Y <- as.matrix(Y)
    dimnames(Y) <- Ynames
  }
  
  input_checker(X, Y)
  
  if(length(n)>1 | length(nx)>1 | length(ny)>1)
    stop("Number of components should be scalars, not vectors\n")
  if(ncol(X) < n + max(nx, ny) || ncol(Y) < n + max(nx, ny)) 
    stop("n + max(nx, ny) = ", n + max(nx, ny), " exceeds #columns in X or Y\n")
  if(nx != round(abs(nx)) || ny != round(abs(ny))) 
    stop("n, nx and ny should be non-negative integers\n")
  if(p_thresh != round(abs(p_thresh)) || q_thresh != round(abs(q_thresh))) 
    stop("p_thresh and q_thresh should be non-negative integers\n")
  if(max_iterations != round(abs(max_iterations)) ) 
    stop("max_iterations should be a non-negative integer\n")
  if(tol < 0) 
    stop("tol should be non-negative\n")
  if(nrow(X) < n + max(nx, ny)) 
    stop("n + max(nx, ny) = ", n + max(nx, ny), " exceed sample size N = ",nrow(X),"\n")
  if(nrow(X) == n + max(nx, ny)) 
    warning("n + max(nx, ny) = ", n + max(nx, ny)," equals sample size\n")
  if (n != round(abs(n)) || n <= 0) {
    stop("n should be a positive integer\n")
  }
  
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...\n")}
  
  if(!sparse){
    ssqX = ssq(X)
    ssqY = ssq(Y)
    
    highd = FALSE
    if ((ncol(X) > p_thresh && ncol(Y) > q_thresh)) {
      highd = TRUE
      message("Using high dimensional mode with tolerance ",tol," and max iterations ",max_iterations,"\n")
      model = o2m2(X, Y, n, nx, ny, stripped, tol, max_iterations)
      class(model) <- "o2m"
    } else if(stripped){
      model = o2m_stripped(X, Y, n, nx, ny)
      class(model) <- "o2m"
    } else {
      
      X_true <- X
      Y_true <- Y
      
      N <- nrow(X)
      p <- ncol(X)
      q <- ncol(Y)
      
      T_Yosc <- U_Xosc <- matrix(0, N, n)
      W_Yosc <- P_Yosc <- matrix(0, p, n)
      C_Xosc <- P_Xosc <- matrix(0, q, n)
      
      if (nx + ny > 0) {
        # larger principal subspace
        n2 <- n + max(nx, ny)
        
        cdw <- svd(t(Y) %*% X, nu = n2, nv = n2)
        C <- cdw$u
        W <- cdw$v
        
        Tt <- X %*% W
        
        if (nx > 0) {
          # Orthogonal components in Y
          E_XY <- X - Tt %*% t(W)
          
          udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
          W_Yosc <- udv$u
          T_Yosc <- X %*% W_Yosc
          P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
          X <- X - T_Yosc %*% t(P_Yosc)
          
          # Update T again
          Tt <- X %*% W
        }
        
        U <- Y %*% C
        
        if (ny > 0) {
          # Orthogonal components in Y
          F_XY <- Y - U %*% t(C)
          
          udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
          C_Xosc <- udv$u
          U_Xosc <- Y %*% C_Xosc
          P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
          Y <- Y - U_Xosc %*% t(P_Xosc)
          
          # Update U again
          U <- Y %*% C
        }
      }
      # Re-estimate joint part in n-dimensional subspace
      cdw <- svd(t(Y) %*% X, nu = n, nv = n)
      C <- cdw$u
      W <- cdw$v
      Tt <- X %*% W
      U <- Y %*% C
      
      # Inner relation parameters
      B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
      B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
      
      # Residuals and R2's
      E <- X_true - Tt %*% t(W) - T_Yosc %*% t(P_Yosc)
      Ff <- Y_true - U %*% t(C) - U_Xosc %*% t(P_Xosc)
      H_TU <- Tt - U %*% B_U
      H_UT <- U - Tt %*% B_T
      Y_hat <- Tt %*% B_T %*% t(C)
      X_hat <- U %*% B_U %*% t(W)
      
      R2Xcorr <- (ssq(Tt)/ssqX)
      R2Ycorr <- (ssq(U)/ssqY)
      R2X_YO <- (ssq(T_Yosc %*% t(P_Yosc))/ssqX)
      R2Y_XO <- (ssq(U_Xosc %*% t(P_Xosc))/ssqY)
      R2Xhat <- (ssq(U %*% B_U)/ssqX)
      R2Yhat <- (ssq(Tt %*% B_T)/ssqY)
      R2X <- R2Xcorr + R2X_YO
      R2Y <- R2Ycorr + R2Y_XO
      
      rownames(Tt) <- rownames(T_Yosc) <- rownames(E) <- rownames(H_TU) <- Xnames[[1]]
      rownames(U) <- rownames(U_Xosc) <- rownames(Ff) <- rownames(H_UT) <- Ynames[[1]]
      rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- colnames(E) <- Xnames[[2]]
      rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- colnames(Ff) <- Ynames[[2]]
      model <- list(Tt = Tt, W. = W, U = U, C. = C, E = E, Ff = Ff, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
                    U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
                    X_hat = X_hat, Y_hat = Y_hat, R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
                    R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat,
                    W_gr = NULL, C_gr = NULL)
      class(model) <- "o2m"
    }

    model$flags = c(list(n = n, nx = nx, ny = ny, 
                         stripped = stripped, highd = highd, 
                         ssqX = ssqX, ssqY = ssqY,
                         varXjoint = apply(model$Tt,2,ssq),
                         varYjoint = apply(model$U,2,ssq),
                         varXorth = apply(model$P_Y,2,ssq)*apply(model$T_Y,2,ssq),
                         varYorth = apply(model$P_X,2,ssq)*apply(model$U_X,2,ssq),
                         method = "O2PLS"))
  } else {
    model = so2m_group(X, Y, n, nx, ny, groupx, groupy, keepx, keepy, 
                             tol, max_iterations, max_iterations_sparsity)
  }
  toc <- proc.time() - tic
  model$flags$call = match.call()
  model$flags$time = toc[3]
  return(model)
}

#' Power method for PLS2
#' 
#' DEFUNCT!!
#' 
#' Performs power method for PLS2 loadings.
#' 
#' @inheritParams o2m
#' @seealso \code{\link{o2m}}
#' @keywords internal
#' @export
pow_o2m2 <- function(X, Y, n, tol = 1e-10, max_iterations = 1000) {
  .Defunct(new = "pow_o2m", msg = "pow_o2m2 is defunct. It was based on the traditional power method. Please use pow_o2m for a NIPALS-based implementation.\n")
  # input_checker(X, Y)
  # stopifnot(n == round(n))
  # #  message("High dimensional problem: switching to power method.\n")
  # #  message("initialize Power Method. Stopping crit: sq.err<", tol, " or ", max_iterations, " iter.\n")
  # Tt <- NULL
  # U <- NULL
  # W <- NULL
  # C <- NULL
  # for (indx in 1:n) {
  #   W0 <- svd(X,nu=0,nv=1)$v[,1]
  #   C0 <- svd(Y,nu=0,nv=1)$v[,1]
  #   for (indx2 in 1:max_iterations) {
  #     tmpp <- c(W0, C0)
  #     W0 <- orth(t(X) %*% (Y %*% t(Y)) %*% (X %*% W0))
  #     C0 <- orth(t(Y) %*% (X %*% t(X)) %*% (Y %*% C0))
  #     if (mse(tmpp, c(W0, C0)) < tol) {
  #       break
  #     }
  #   }
  #   if(ssq(W0) < 1e-10 || ssq(C0) < 1e-10){
  #     W0 <- orth(rep(1,ncol(X)))
  #     C0 <- orth(rep(1,ncol(Y)))
  #     for (indx2 in 1:max_iterations) {
  #       tmpp <- c(W0, C0)
  #       W0 <- orth(t(X) %*% (Y %*% t(Y)) %*% (X %*% W0))
  #       C0 <- orth(t(Y) %*% (X %*% t(X)) %*% (Y %*% C0))
  #       if (mse(tmpp, c(W0, C0)) < tol) {
  #         message("The initialization of the power method lied in a degenerate space\n")
  #         message("Initialization changed and power method rerun\n")
  #         break
  #       }
  #     }
  #   }
  #   message("Power Method (comp ", indx, ") stopped after ", indx2, " iterations.\n")
  #   Tt <- cbind(Tt, X %*% W0)
  #   U <- cbind(U, Y %*% C0)
  #   X <- X - (X %*% W0) %*% t(W0)
  #   Y <- Y - (Y %*% C0) %*% t(C0)
  #   W <- cbind(W, W0)
  #   C <- cbind(C, C0)
  # }
  # return(list(W = W, C = C, Tt = Tt, U = U))
}

#' NIPALS method for PLS2
#' 
#' Performs power method for PLS2 loadings.
#' 
#' @inheritParams o2m
#' @seealso \code{\link{o2m}}
#' @keywords internal
#' @export
pow_o2m <- function(X, Y, n, tol = 1e-10, max_iterations = 1000) {
  input_checker(X, Y)
  stopifnot(n == round(n))
  #  message("High dimensional problem: switching to power method.\n")
  #  message("initialize Power Method. Stopping crit: sq.err<", tol, " or ", max_iterations, " iter.\n")
  Tt <- NULL
  U <- NULL
  W <- NULL
  C <- NULL
  for (indx in 1:n) {
    Ui = rowSums(Y)
    for (indx2 in 1:max_iterations) {
      tmpp <- Ui
      Wi = crossprod(X, Ui)
      Wi = Wi / c(vnorm(Wi))
      Ti = X %*% Wi
      Ci = crossprod(Y, Ti)
      Ci = Ci / c(vnorm(Ci))
      Ui = Y %*% Ci
      if (mse(tmpp, Ui) < tol) {
        break
      }
    }
    if(ssq(Wi) < 1e-10 || ssq(Ci) < 1e-10){
      Wi <- orth(rep(1,ncol(X)))
      Ci <- orth(rep(1,ncol(Y)))
      for (indx2 in 1:max_iterations) {
        tmpp <- c(Wi, Ci)
        Wi <- orth(t(X) %*% (Y %*% t(Y)) %*% (X %*% Wi))
        Ci <- orth(t(Y) %*% (X %*% t(X)) %*% (Y %*% Ci))
        if (mse(tmpp, c(Wi, Ci)) < tol) {
          message("The initialization of the power method lied in a degenerate space\n")
          message("Initialization changed and power method rerun\n")
          break
        }
      }
    }
    message("Power Method (comp ", indx, ") stopped after ", indx2, " iterations.\n")
    Tt <- cbind(Tt, X %*% Wi)
    U <- cbind(U, Y %*% Ci)
    X <- X - (X %*% Wi) %*% t(Wi)
    Y <- Y - (Y %*% Ci) %*% t(Ci)
    W <- cbind(W, Wi)
    C <- cbind(C, Ci)
  }
  return(list(W = W, C = C, Tt = Tt, U = U))
}

#' Perform O2-PLS with two-way orthogonal corrections
#'
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALES THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#'
#' @inheritParams o2m
#' 
#' @return A list containing
#'    \item{Tt}{Joint \eqn{X} scores}
#'    \item{W.}{Joint \eqn{X} loadings}
#'    \item{U}{Joint \eqn{Y} scores}
#'    \item{C.}{Joint \eqn{Y} loadings}
#'    \item{E}{Residuals in \eqn{X}}
#'    \item{Ff}{Residuals in \eqn{Y}}
#'    \item{T_Yosc}{Orthogonal \eqn{X} scores}
#'    \item{P_Yosc.}{Orthogonal \eqn{X} loadings}
#'    \item{W_Yosc}{Orthogonal \eqn{X} weights}
#'    \item{U_Xosc}{Orthogonal \eqn{Y} scores}
#'    \item{P_Xosc.}{Orthogonal \eqn{Y} loadings}
#'    \item{C_Xosc}{Orthogonal \eqn{Y} weights}
#'    \item{B_U}{Regression coefficient in \code{Tt} ~ \code{U}}
#'    \item{B_T.}{Regression coefficient in \code{U} ~ \code{Tt}}
#'    \item{H_TU}{Residuals in \code{Tt} in \code{Tt} ~ \code{U}}
#'    \item{H_UT}{Residuals in \code{U} in \code{U} ~ \code{Tt}}
#'    \item{X_hat}{Prediction of \eqn{X} with \eqn{Y}}
#'    \item{Y_hat}{Prediction of \eqn{Y} with \eqn{X}}
#'    \item{R2X}{Variation (measured with \code{\link{ssq}}) of the modeled part in \eqn{X} (defined by joint + orthogonal variation) as proportion of variation in \eqn{X}}
#'    \item{R2Y}{Variation (measured with \code{\link{ssq}}) of the modeled part in \eqn{Y} (defined by joint + orthogonal variation) as proportion of variation in \eqn{Y}}
#'    \item{R2Xcorr}{Variation (measured with \code{\link{ssq}}) of the joint part in \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Ycorr}{Variation (measured with \code{\link{ssq}}) of the joint part in \eqn{Y} as proportion of variation in \eqn{Y}}
#'    \item{R2X_YO}{Variation (measured with \code{\link{ssq}}) of the orthogonal part in \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Y_XO}{Variation (measured with \code{\link{ssq}}) of the orthogonal part in \eqn{Y} as proportion of variation in \eqn{Y}}
#'    \item{R2Xhat}{Variation (measured with \code{\link{ssq}}) of the predicted \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Yhat}{Variation (measured with \code{\link{ssq}}) of the predicted \eqn{Y} as proportion of variation in \eqn{Y}}
#'
#' @details If both \code{nx} and \code{ny} are zero, \code{o2m2} is equivalent to PLS2 with orthonormal loadings.
#' For cross-validation purposes, consider using \code{stripped = TRUE}.
#' 
#' Note that in this function, a power-method based approach is used when the data dimensionality is larger than the sample size. E.g. for genomic data the covariance matrix might be too memory expensive.
#'
#' @examples
#' # This takes a couple of seconds on an intel i5
#' system.time(
#'              o2m2(matrix(rnorm(50*2000),50),matrix(rnorm(50*2000),50),1,0,0)
#' )
#' 
#' # This however takes 10 times as much...
#' # system.time(
#' #              o2m(matrix(rnorm(50*2000),50),matrix(rnorm(50*2000),50),1,0,0,
#' #              p_thresh = 1e4,q_thresh = 1e4)  # makes sure power method is not used
#' # )
#'
#' @seealso \code{\link{o2m}}
#' 
#' @keywords internal
#' @export
o2m2 <- function(X, Y, n, nx, ny, stripped = TRUE, tol = 1e-10, max_iterations = 100) {
  
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  X_true <- X
  Y_true <- Y
  
  N <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  T_Yosc <- U_Xosc <- matrix(0, N, n)
  W_Yosc <- P_Yosc <- matrix(0, p, n)
  C_Xosc <- P_Xosc <- matrix(0, q, n)
  
  if (nx + ny > 0) {
    # larger principal subspace
    n2 <- n + max(nx, ny)
    
    # if(N<p&N<q){ # When N is smaller than p and q
    W_C <- pow_o2m(X, Y, n2, tol, max_iterations)
    W <- W_C$W
    C <- W_C$C
    Tt <- W_C$Tt
    U <- W_C$U
    # } cdw = svd(t(Y)%*%X,nu=n2,nv=n2); C=cdw$uW=cdw$v
    
    # Tt = X%*%W;
    
    if (nx > 0) {
      # Orthogonal components in Y
      E_XY <- X - Tt %*% t(W)
      
      udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
      W_Yosc <- udv$u
      T_Yosc <- X %*% W_Yosc
      P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
      X <- X - T_Yosc %*% t(P_Yosc)
      
      # Update T again Tt = X%*%W;
    }
    
    # U = Y%*%C; # 3.2.1. 4
    
    if (ny > 0) {
      # Orthogonal components in Y
      F_XY <- Y - U %*% t(C)
      
      udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
      C_Xosc <- udv$u
      U_Xosc <- Y %*% C_Xosc
      P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
      Y <- Y - U_Xosc %*% t(P_Xosc)
      
      # Update U again U = Y%*%C;
    }
  }
  # Re-estimate joint part in n-dimensional subspace if(N<p&N<q){ # When N is smaller than p and q
  W_C <- pow_o2m(X, Y, n, tol, max_iterations)
  W <- W_C$W
  C <- W_C$C
  Tt <- W_C$Tt
  U <- W_C$U
  # } cdw = svd(t(Y)%*%X,nu=n,nv=n); # 3.2.1. 1 C=cdw$u;W=cdw$v Tt = X%*%W; # 3.2.1. 2 U = Y%*%C; #
  # 3.2.1. 4
  
  # Inner relation parameters
  B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  
  # Residuals and R2's
  if(stripped){
    E <- Ff <- X_hat <- Y_hat <- as.matrix(0)
  } else {
    E <- X_true - Tt %*% t(W) - T_Yosc %*% t(P_Yosc)
    Ff <- Y_true - U %*% t(C) - U_Xosc %*% t(P_Xosc)
    Y_hat <- Tt %*% B_T %*% t(C)
    X_hat <- U %*% B_U %*% t(W)
  }
  H_TU <- Tt - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  # R2
  R2Xcorr <- (ssq(Tt)/ssq(X_true))
  R2Ycorr <- (ssq(U)/ssq(Y_true))
  R2X_YO <- (ssq(T_Yosc %*% t(P_Yosc))/ssq(X_true))
  R2Y_XO <- (ssq(U_Xosc %*% t(P_Xosc))/ssq(Y_true))
  R2Xhat <- (ssq(U %*% B_U)/ssq(X_true))
  R2Yhat <- (ssq(Tt %*% B_T)/ssq(Y_true))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_UT) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, W. = W, U = U, C. = C, E = E, Ff = Ff, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
                U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
                X_hat = X_hat, Y_hat = Y_hat, R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
                R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  class(model) <- "pre.o2m"
  return(model)
}

#' Perform O2-PLS with two-way orthogonal corrections
#'
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALES THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#' A stripped version of O2PLS
#'
#' @inheritParams o2m
#' 
#' @return A list containing
#'    \item{Tt}{Joint \eqn{X} scores}
#'    \item{W.}{Joint \eqn{X} loadings}
#'    \item{U}{Joint \eqn{Y} scores}
#'    \item{C.}{Joint \eqn{Y} loadings}
#'    \item{P_Yosc.}{Orthogonal \eqn{X} loadings}
#'    \item{P_Xosc.}{Orthogonal \eqn{Y} loadings}
#'    \item{B_U}{Regression coefficient in \code{Tt} ~ \code{U}}
#'    \item{B_T.}{Regression coefficient in \code{U} ~ \code{Tt}}
#'    \item{H_TU}{Residuals in \code{Tt} in \code{Tt} ~ \code{U}}
#'    \item{H_UT}{Residuals in \code{U} in \code{U} ~ \code{Tt}}
#'
#' @details If both \code{nx} and \code{ny} are zero, \code{o2m} is equivalent to PLS2 with orthonormal loadings.
#' This is a stripped implementation of O2PLS, using \code{\link{svd}}. For data analysis purposes, consider using \code{\link{o2m}}.
#'
#' @seealso \code{\link{o2m}}
#' @keywords internal
#' @export
o2m_stripped <- function(X, Y, n, nx, ny) {
  
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  X_true <- X
  Y_true <- Y
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  T_Yosc <- U_Xosc <- matrix(0, N, 1)
  P_Yosc <- W_Yosc <- matrix(0, p, 1)
  P_Xosc <- C_Xosc <- matrix(0, q, 1)
  
  if (nx + ny > 0) {
    n2 <- n + max(nx, ny)
    
    cdw <- svd(t(Y) %*% X, nu = n2, nv = n2)
    C <- cdw$u
    W <- cdw$v
    rm(cdw)
    Tt <- X %*% W
    
    if (nx > 0) {
      # 3.2.1. 3
      E_XY <- X
      E_XY <- E_XY - Tt %*% t(W)
      
      udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
      rm(E_XY)
      W_Yosc <- udv$u
      T_Yosc <- X %*% W_Yosc
      P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
      X <- X - T_Yosc %*% t(P_Yosc)
      
      # Update T again (since X has changed)
      Tt <- X %*% W
    }
    
    U <- Y %*% C
    
    if (ny > 0) {
      # 3.2.1. 5
      F_XY <- Y
      F_XY <- F_XY - U %*% t(C)
      
      udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
      rm(F_XY)
      C_Xosc <- udv$u
      U_Xosc <- Y %*% C_Xosc
      P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
      Y <- Y - U_Xosc %*% t(P_Xosc)
      
      # Update U again (since Y has changed)
      U <- Y %*% C
    }
  }
  
  # repeat steps 1, 2, and 4 before step 6
  cdw <- svd(t(Y) %*% X, nu = n, nv = n)
  C <- cdw$u
  W <- cdw$v
  Tt <- X %*% W
  U <- Y %*% C
  
  # 3.2.1. 6
  B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  H_TU <- Tt - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  R2Xcorr <- ssq(Tt) / ssq(X_true)
  R2Ycorr <- ssq(U) / ssq(Y_true)
  R2X_YO <- ssq(T_Yosc %*% t(P_Yosc)) / ssq(X_true)
  R2Y_XO <- ssq(U_Xosc %*% t(P_Xosc)) / ssq(Y_true)
  R2Xhat <- (ssq(U %*% B_U) / ssq(X_true))
  R2Yhat <- (ssq(Tt %*% B_T) / ssq(Y_true))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_TU) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, U = U, W. = W, C. = C, P_Yosc. = P_Yosc, P_Xosc. = P_Xosc,
                T_Yosc = T_Yosc, U_Xosc = U_Xosc, W_Yosc = W_Yosc, C_Xosc = C_Xosc,
                B_T. = B_T, B_U = B_U, H_TU = H_TU, H_UT = H_UT, 
                R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, 
                R2Xhat = R2Xhat, R2Yhat = R2Yhat,
                W_gr = NULL, C_gr = NULL)
  class(model) <- c("pre.o2m","o2m_stripped")
  return(model)
}


#' Perform O2-PLS with two-way orthogonal corrections
#'
#' DEFUNCT!!
#' 
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALES THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#' A stripped version of O2PLS
#'
#' @inheritParams o2m
#' 
#' @return A list containing
#'    \item{Tt}{Joint \eqn{X} scores}
#'    \item{W.}{Joint \eqn{X} loadings}
#'    \item{U}{Joint \eqn{Y} scores}
#'    \item{C.}{Joint \eqn{Y} loadings}
#'    \item{P_Yosc.}{Orthogonal \eqn{X} loadings}
#'    \item{P_Xosc.}{Orthogonal \eqn{Y} loadings}
#'    \item{B_U}{Regression coefficient in \code{Tt} ~ \code{U}}
#'    \item{B_T.}{Regression coefficient in \code{U} ~ \code{Tt}}
#'    \item{H_TU}{Residuals in \code{Tt} in \code{Tt} ~ \code{U}}
#'    \item{H_UT}{Residuals in \code{U} in \code{U} ~ \code{Tt}}
#'
#' @details If both \code{nx} and \code{ny} are zero, \code{o2m} is equivalent to PLS2 with orthonormal loadings.
#' This is a stripped implementation of O2PLS, using \code{\link{svd}}. For data analysis purposes, consider using \code{\link{o2m}}.
#'
#' @seealso \code{\link{o2m}}
#' @keywords internal
#' @export
o2m_stripped2 <- function(X, Y, n, nx, ny, tol = 1e-10, max_iterations = 100) {
  .Defunct(new="o2m_stripped")
  # Xnames = dimnames(X)
  # Ynames = dimnames(Y)
  # 
  # X_true <- X
  # Y_true <- Y
  # 
  # N <- dim(X)[1]
  # p <- dim(X)[2]
  # q <- dim(Y)[2]
  # 
  # T_Yosc <- U_Xosc <- matrix(0, N, 1)
  # P_Yosc <- W_Yosc <- matrix(0, p, 1)
  # P_Xosc <- C_Xosc <- matrix(0, q, 1)
  # 
  # if (nx + ny > 0) {
  #   n2 <- n + max(nx, ny)
  #   
  #   #    if (N < p & N < q) {
  #   # When N is smaller than p and q
  #   W_C <- pow_o2m(X, Y, n2, tol, max_iterations)
  #   W <- W_C$W
  #   C <- W_C$C
  #   Tt <- W_C$Tt
  #   U <- W_C$U
  #   rm(W_C)
  #   gc()
  #   #    }
  #   # 3.2.1. 2
  #   
  #   if (nx > 0) {
  #     # 3.2.1. 3
  #     E_XY <- X - Tt %*% t(W)
  #     
  #     udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
  #     rm(E_XY)
  #     W_Yosc <- udv$u
  #     T_Yosc <- X %*% W_Yosc
  #     P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
  #     X <- X - T_Yosc %*% t(P_Yosc)
  #     
  #     # Update T again (since X has changed) Tt = X%*%W;
  #   }
  #   
  #   # U = Y%*%C; # 3.2.1. 4
  #   
  #   if (ny > 0) {
  #     # 3.2.1. 5
  #     F_XY <- Y
  #     F_XY <- F_XY - U %*% t(C)
  #     
  #     udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
  #     rm(F_XY)
  #     C_Xosc <- udv$u
  #     U_Xosc <- Y %*% C_Xosc
  #     P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
  #     Y <- Y - U_Xosc %*% t(P_Xosc)
  #     
  #     # Update U again (since Y has changed) U = Y%*%C;
  #   }
  # }
  # 
  # # repeat steps 1, 2, and 4 before step 6 When N is smaller than p and q
  # #  if (N < p & N < q) {
  # W_C <- pow_o2m(X, Y, n, tol, max_iterations)
  # W <- W_C$W
  # C <- W_C$C
  # Tt <- W_C$Tt
  # U <- W_C$U
  # rm(W_C)
  # gc()
  # #  }
  # 
  # # 3.2.1. 6
  # B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  # B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  # H_TU <- Tt - U %*% B_U
  # H_UT <- U - Tt %*% B_T
  # 
  # R2Xcorr <- ssq(Tt) / ssq(X_true)
  # R2Ycorr <- ssq(U) / ssq(Y_true)
  # R2X_YO <- ssq(T_Yosc %*% t(P_Yosc)) / ssq(X_true)
  # R2Y_XO <- ssq(U_Xosc %*% t(P_Xosc)) / ssq(Y_true)
  # R2Xhat <- (ssq(U %*% B_U) / ssq(X_true))
  # R2Yhat <- (ssq(Tt %*% B_T) / ssq(Y_true))
  # R2X <- R2Xcorr + R2X_YO
  # R2Y <- R2Ycorr + R2Y_XO
  # 
  # rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  # rownames(U) <- rownames(U_Xosc) <- rownames(H_TU) <- Ynames[[1]]
  # rownames(W) <- rownames(P_Yosc) <- Xnames[[2]]
  # rownames(C) <- rownames(P_Xosc) <- Ynames[[2]]
  # 
  # model <- list(Tt = Tt, U = U, W. = W, C. = C, P_Yosc. = P_Yosc, P_Xosc. = P_Xosc,
  #               T_Yosc = T_Yosc, U_Xosc = U_Xosc, W_Yosc = W_Yosc, C_Xosc = C_Xosc,
  #               B_T. = B_T, B_U = B_U, H_TU = H_TU, H_UT = H_UT, 
  #               R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, 
  #               R2Xhat = R2Xhat, R2Yhat = R2Yhat,
  #               W_gr = NULL, C_gr = NULL)
  # 
  # class(model) <- c("pre.o2m","o2m_stripped")
  # return(model)
}


#' Perform Group Sparse O2PLS 
#'
#' @inheritParams o2m
#'
#' @return A list containing
#'    \item{Tt}{Joint \eqn{X} scores}
#'    \item{W.}{Joint \eqn{X} loadings}
#'    \item{U}{Joint \eqn{Y} scores}
#'    \item{C.}{Joint \eqn{Y} loadings}
#'    \item{T_Yosc}{Orthogonal \eqn{X} scores}
#'    \item{P_Yosc.}{Orthogonal \eqn{X} loadings}
#'    \item{W_Yosc}{Orthogonal \eqn{X} weights}
#'    \item{U_Xosc}{Orthogonal \eqn{Y} scores}
#'    \item{P_Xosc.}{Orthogonal \eqn{Y} loadings}
#'    \item{C_Xosc}{Orthogonal \eqn{Y} weights}
#'    \item{B_U}{Regression coefficient in \code{Tt} ~ \code{U}}
#'    \item{B_T.}{Regression coefficient in \code{U} ~ \code{Tt}}
#'    \item{H_TU}{Residuals in \code{Tt} in \code{Tt} ~ \code{U}}
#'    \item{H_UT}{Residuals in \code{U} in \code{U} ~ \code{Tt}}
#'    \item{W_gr}{Joint weights of X variables at group level. They are the norms of the X-joint loadings per group}
#'    \item{C_gr}{Joint weights of Y variables at group level. They are the norms of the Y-joint loadings per group}
#' 
#' @keywords internal
#' @seealso \code{\link{o2m}}
#' @export
so2m_group <- function(X, Y, n, nx, ny, groupx=NULL, groupy=NULL, keepx=NULL, keepy=NULL, 
                       tol = 1e-10, max_iterations=1000, max_iterations_sparsity=1000){

  if(is.null(groupx) & is.null(groupy)){
    method = "SO2PLS"
    #message("Group information not provided, using SO2PLS")
    keepxy <- lambda_checker(X, Y, keepx, keepy, n)
    keepx <- keepxy$keepx
    keepy <- keepxy$keepy
  }else{
    method = "GO2PLS"
    #message("Group information provided, using GO2PLS")
    # check if only information for one dataset is provided
    if(is.null(groupx)){
      if(is.null(colnames(X))) stop("Please provide 'groupx' or colnames of X")
      groupx = colnames(X)
    }
    if(is.null(groupy)){
      if(is.null(colnames(Y))) stop("Please provide 'groupy' or colnames of Y")
      groupy = colnames(Y)
    }
    keepxy <- lambda_checker_group(groupx, groupy, keepx, keepy, n)
    keepx <- keepxy$keepx
    keepy <- keepxy$keepy
  }
  
  ssqX = ssq(X)
  ssqY = ssq(Y)
  #############################################################
  # Orthogonal filtering
  #############################################################
  # setup
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  X_true <- X
  Y_true <- Y
  
  N <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  T_Yosc <- U_Xosc <- matrix(0, N, n)
  W_Yosc <- P_Yosc <- matrix(0, p, n)
  C_Xosc <- P_Xosc <- matrix(0, q, n)
  
  # filtering
  if (nx + ny > 0) {
    # larger principal subspace
    n2 <- n + max(nx, ny)
    
    W_C <- pow_o2m(X, Y, n2, tol, max_iterations)
    W <- W_C$W
    C <- W_C$C
    Tt <- W_C$Tt
    U <- W_C$U

    # Tt = X%*%W;
    
    if (nx > 0) {
      # Orthogonal components in Y
      E_XY <- X - Tt %*% t(W)
      
      udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
      W_Yosc <- udv$u
      T_Yosc <- X %*% W_Yosc
      P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
      X <- X - T_Yosc %*% t(P_Yosc)
      
      # Update T again Tt = X%*%W;
    }
    
    # U = Y%*%C; # 3.2.1. 4
    
    if (ny > 0) {
      # Orthogonal components in Y
      F_XY <- Y - U %*% t(C)
      
      udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
      C_Xosc <- udv$u
      U_Xosc <- Y %*% C_Xosc
      P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
      Y <- Y - U_Xosc %*% t(P_Xosc)
      
      # Update U again U = Y%*%C;
    }
  } 
  #############################################################
  W <- matrix(0, dim(X)[2], n)
  C <- matrix(0, dim(Y)[2], n)
  Tt <- matrix(0, dim(X)[1], n)
  U <- matrix(0, dim(Y)[1], n)
  
  if(method == "SO2PLS"){
    W_C <- pow_o2m(X, Y, n, tol, max_iterations)
    w_ini <- W_C$W
    
    # get u,v iteratively
    for(j in 1: n){
      v <- w_ini[,j]
      for (i in 1: max_iterations_sparsity){
        v_old <- v
        u <- t(Y) %*% (X %*% v)
        u <- thresh_n(u, keepy[j])
        u <- u/norm_vec(u)
        v <- t(X) %*% (Y %*% u)
        v <- thresh_n(v, keepx[j])
        v <- v/norm_vec(v)
        if (mse(v, v_old) < tol) {
          break
        }
      }
      
      # post-orthogonalizing
      if(j>1){
        # message('W')
        v <- orth_vec(v, W[,1:j-1])
        # message('C')
        u <- orth_vec(u, C[,1:j-1])
      }
      
      W[,j] <- v
      C[,j] <- u
      Tt[,j] <- X %*% W[,j]
      U[,j] <- Y %*% C[,j]
      
      p <- t(X) %*% Tt[,j] / drop(crossprod(Tt[,j]))
      q <- t(Y) %*% U[,j] / drop(crossprod(U[,j]))
      X <- X - Tt[,j] %*% t(p)
      Y <- Y - U[,j] %*% t(q)
    }
    select_grx <- NULL
    select_gry <- NULL
  }
  
  if(method == "GO2PLS"){
    names_grx <- groupx %>% unique # names of groups
    names_gry <- groupy %>% unique  
    nr_grx <- names_grx %>% length # number of groups
    nr_gry <- names_gry %>% length
    index_grx <- lapply(1:nr_grx, function(j){
      index <- which(groupx == names_grx[j])
      size <- length(index)
      return(list(index=index, size=size))
    })
    index_gry <- lapply(1:nr_gry, function(j){
      index <- which(groupy == names_gry[j])
      size <- length(index)
      return(list(index=index, size=size))
    })
    names(index_grx) <- names_grx
    names(index_gry) <- names_gry 
    
    select_grx <- select_gry <- list()
    
    W_C <- pow_o2m(X, Y, n, tol, max_iterations)
    w_ini <- W_C$W
    
    # get u,v iteratively
    for(j in 1: n){
      if(length(keepx)==1){keepx <- rep(keepx,n)}
      if(length(keepy)==1){keepy <- rep(keepy,n)}
      v <- w_ini[,j]
      for (i in 1: max_iterations_sparsity){
        v_old <- v
        u <- t(Y) %*% (X %*% v)
        ul <- thresh_n_gr(u, keepy[j], index_gry)
        u <- ul$w
        u <- u/norm_vec(u)
        v <- t(X) %*% (Y %*% u)
        vl <- thresh_n_gr(v, keepx[j], index_grx)
        v <- vl$w
        v <- v/norm_vec(v)
        if (mse(v, v_old) < tol) {
          select_grx[[j]] <- sapply(1:length(index_grx), function(k){
            wj <- v[index_grx[[k]]$index] 
            normj <- norm_vec(wj)
            return(normj)
          })
          select_gry[[j]] <- sapply(1:length(index_gry), function(k){
            wj <- u[index_gry[[k]]$index] 
            normj <- norm_vec(wj)
            return(normj)
          })
          
          names(select_grx[[j]]) <- names(index_grx)
          names(select_gry[[j]]) <- names(index_gry)
          
          break
        }
      }
      
      # post-orthogonalizing
      if(j>1){
        # message('W')
        v <- orth_vec(v, W[,1:j-1])
        # message('C')
        u <- orth_vec(u, C[,1:j-1])
      }
      
      W[,j] <- v
      C[,j] <- u
      Tt[,j] <- X %*% W[,j]
      U[,j] <- Y %*% C[,j]
      
      p <- t(X) %*% Tt[,j] / drop(crossprod(Tt[,j]))
      q <- t(Y) %*% U[,j] / drop(crossprod(U[,j]))
      X <- X - Tt[,j] %*% t(p)
      Y <- Y - U[,j] %*% t(q)
    }
    
    select_grx <- do.call(cbind,select_grx)
    select_gry <- do.call(cbind,select_gry)
  }
  
  # Inner relation parameters
  B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  
  # Residuals and R2's
  #E <- X_true - Tt %*% t(W) - T_Yosc %*% t(P_Yosc)
  #Ff <- Y_true - U %*% t(C) - U_Xosc %*% t(P_Xosc)
  H_TU <- Tt - U %*% B_U
  H_UT <- U - Tt %*% B_T
  #Y_hat <- Tt %*% B_T %*% t(C)
  #X_hat <- U %*% B_U %*% t(W)
  
  R2Xcorr <- (ssq(Tt)/ssqX)
  R2Ycorr <- (ssq(U)/ssqY)
  R2X_YO <- (ssq(T_Yosc %*% t(P_Yosc))/ssqX)
  R2Y_XO <- (ssq(U_Xosc %*% t(P_Xosc))/ssqY)
  R2Xhat <- (ssq(U %*% B_U)/ssqX)
  R2Yhat <- (ssq(Tt %*% B_T)/ssqY)
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_UT) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, W. = W, U = U, C. = C, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
                U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
                R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
                R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat, 
                W_gr = select_grx, C_gr = select_gry)
  class(model) <- "o2m"
  model$flags = c(list(n = n, nx = nx, ny = ny, 
                       stripped = FALSE, highd = TRUE, 
                       ssqX = ssqX, ssqY = ssqY,
                       varXjoint = apply(model$Tt,2,ssq),
                       varYjoint = apply(model$U,2,ssq),
                       varXorth = apply(model$P_Y,2,ssq)*apply(model$T_Y,2,ssq),
                       varYorth = apply(model$P_X,2,ssq)*apply(model$U_X,2,ssq),
                       method = method))
  return(model)
}
