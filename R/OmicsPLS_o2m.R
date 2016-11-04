#' Perform O2-PLS with two-way orthogonal corrections
#'
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALES THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#'
#' @param X Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param Y Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param n Integer. Number of joint PLS components. Must be positive!
#' @param nx Integer. Number of orthogonal components in \eqn{X}. Negative values are interpreted as 0
#' @param ny Integer. Number of orthogonal components in \eqn{Y}. Negative values are interpreted as 0
#' @param stripped Logical. Use the stripped version of o2m (usually when cross-validating)?
#' @param p_thresh Integer. If \code{X} has more than \code{p_thresh} columns, a power method optimization is used, see \code{\link{o2m2}}
#' @param q_thresh Integer. If \code{Y} has more than \code{q_thresh} columns, a power method optimization is used, see \code{\link{o2m2}}
#' @param tol double. Threshold for power method iteration
#' @param max_iterations Integer, Maximum number of iterations for power method
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
#' @details If both \code{nx} and \code{ny} are zero, \code{o2m} is equivalent to PLS2 with orthonormal loadings.
#' This is a `slower' (in terms of memory) implementation of O2PLS, and is using \code{\link{svd}}, use \code{stripped=T} for a stripped version with less output.
#' If either \code{ncol(X) > p_thresh} or \code{ncol(Y) > q_thresh}, an alternative method is used (NIPALS) which does not store the entire covariance matrix.
#' The squared error between iterands in the NIPALS approach can be adjusted with \code{tol}.
#' The maximum number of iterations in the NIPALS approach is tuned by \code{max_iterations}.
#'
#' @examples
#' test.data <- scale(matrix(rnorm(100)))
#' hist(replicate(1000,
#'          o2m(test.data,scale(matrix(rnorm(100))),1,0,0)$B_T.
#'      ),main='No joint variation',xlab='B_T',xlim=c(0,0.6));
#' hist(replicate(1000,
#'          o2m(test.data,scale(test.data+rnorm(100))/2,1,0,0)$B_T.
#'     ),main='B_T = 0.5 \n 25% joint variation',xlab='B_T',xlim=c(0,0.6));
#' hist(replicate(1000,
#'          o2m(test.data,scale(test.data+rnorm(100,0,0.1))/2,1,0,0)$B_T.
#'     ),main='B_T = 0.5 \n 90% joint variation',xlab='B_T',xlim=c(0,0.6));
#'
#' @seealso \code{\link{ssq}}, \code{\link{summary.o2m}}, \code{\link{plot.o2m}}, \code{\link{crossval_o2m}}
#'
#' @export
o2m <- function(X, Y, n, nx, ny, stripped = FALSE, 
                p_thresh = 3000, q_thresh = p_thresh, tol = 1e-10, max_iterations = 100) {
  tic <- proc.time()
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  if(!is.matrix(X)){
    message("X has class ",class(X),", trying to convert with as.matrix.",sep="")
    X <- as.matrix(X)
    dimnames(X) <- Xnames
  }
  if(!is.matrix(Y)){
    message("Y has class ",class(Y),", trying to convert with as.matrix.",sep="")
    Y <- as.matrix(Y)
    dimnames(Y) <- Ynames
  }
  input_checker(X, Y)
  
  ssqX = ssq(X)
  ssqY = ssq(Y)
  
  if(ncol(X) < n + max(nx, ny) || ncol(Y) < n + max(nx, ny)) 
    stop("n + max(nx, ny) =", n + max(nx, ny), " exceed # columns in X or Y")
  if(nx != round(abs(nx)) || ny != round(abs(ny))) 
    stop("n, nx and ny should be non-negative integers")
  if(p_thresh != round(abs(p_thresh)) || q_thresh != round(abs(q_thresh))) 
    stop("p_thresh and q_thresh should be non-negative integers")
  if(max_iterations != round(abs(max_iterations)) ) 
    stop("max_iterations should be a non-negative integer")
  if(tol < 0) 
    stop("tol should be non-negative")
  if(nrow(X) < n + max(nx, ny)) 
    stop("n + max(nx, ny) = ", n + max(nx, ny), " exceed sample size N = ",nrow(X))
  if(nrow(X) == n + max(nx, ny)) 
    warning("n + max(nx, ny) = ", n + max(nx, ny)," equals sample size")
  if (n != round(abs(n)) || n <= 0) {
    stop("n should be a positive integer")
  }
  
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...")}
  
  highd = FALSE
  if ((ncol(X) > p_thresh && ncol(Y) > q_thresh)) {
    highd = TRUE
    message("Using Power Method with tolerance ",tol," and max iterations ",max_iterations)
    model = o2m2(X, Y, n, nx, ny, stripped, tol, max_iterations)
  } else if(stripped){
    model = o2m_stripped(X, Y, n, nx, ny)
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
    
    R2Xcorr <- (ssq(Tt %*% t(W))/ssqX)
    R2Ycorr <- (ssq(U %*% t(C))/ssqY)
    R2X_YO <- (ssq(T_Yosc %*% t(P_Yosc))/ssqX)
    R2Y_XO <- (ssq(U_Xosc %*% t(P_Xosc))/ssqY)
    R2Xhat <- 1 - (ssq(U %*% B_U %*% t(W) - X_true)/ssqX)
    R2Yhat <- 1 - (ssq(Tt %*% B_T %*% t(C) - Y_true)/ssqY)
    R2X <- R2Xcorr + R2X_YO
    R2Y <- R2Ycorr + R2Y_XO
    
    rownames(Tt) <- rownames(T_Yosc) <- rownames(E) <- rownames(H_TU) <- Xnames[[1]]
    rownames(U) <- rownames(U_Xosc) <- rownames(Ff) <- rownames(H_UT) <- Ynames[[1]]
    rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- colnames(E) <- Xnames[[2]]
    rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- colnames(Ff) <- Ynames[[2]]
    model <- list(Tt = Tt, W. = W, U = U, C. = C, E = E, Ff = Ff, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
                  U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
                  X_hat = X_hat, Y_hat = Y_hat, R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
                  R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat)
    class(model) <- "o2m"
  }
  toc <- proc.time() - tic
  model$flags = c(time = toc[3], 
                  list(n = n, nx = nx, ny = ny, 
                       stripped = stripped, highd = highd, 
                       call = match.call(), ssqX = ssqX, ssqY = ssqY,
                       varXjoint = apply(model$Tt,2,ssq),
                       varYjoint = apply(model$U,2,ssq),
                       varXorth = apply(model$P_Y,2,ssq)*apply(model$T_Y,2,ssq),
                       varYorth = apply(model$P_X,2,ssq)*apply(model$U_X,2,ssq)))
  return(model)
}

#' Power method for PLS2
#' 
#' Performs power method for PLS2 loadings.
#' 
#' @inheritParams o2m
#' 
#' @keywords internal
#' @export
pow_o2m2 <- function(X, Y, n, tol = 1e-10, max_iterations = 100) {
  input_checker(X, Y)
  stopifnot(n == round(n))
  #  message("High dimensional problem: switching to power method.\n")
  #  message("initialize Power Method. Stopping crit: sq.err<", tol, " or ", max_iterations, " iter.\n")
  Tt <- NULL
  U <- NULL
  W <- NULL
  C <- NULL
  for (indx in 1:n) {
    W0 <- svd(X,nu=0,nv=1)$v[,1]
    C0 <- svd(Y,nu=0,nv=1)$v[,1]
    for (indx2 in 1:max_iterations) {
      tmpp <- c(W0, C0)
      W0 <- orth(t(X) %*% (Y %*% t(Y)) %*% (X %*% W0))
      C0 <- orth(t(Y) %*% (X %*% t(X)) %*% (Y %*% C0))
      if (mse(tmpp, c(W0, C0)) < tol) {
        break
      }
    }
    if(ssq(W0) < 1e-10 || ssq(C0) < 1e-10){
      W0 <- orth(rep(1,ncol(X)))
      C0 <- orth(rep(1,ncol(Y)))
      for (indx2 in 1:max_iterations) {
        tmpp <- c(W0, C0)
        W0 <- orth(t(X) %*% (Y %*% t(Y)) %*% (X %*% W0))
        C0 <- orth(t(Y) %*% (X %*% t(X)) %*% (Y %*% C0))
        if (mse(tmpp, c(W0, C0)) < tol) {
          message("The initialization of the power method lied in a degenerate space\n")
          message("Initialization changed and power method rerun\n")
          break
        }
      }
    }
    message("Power Method (comp ", indx, ") stopped after ", indx2, " iterations.\n")
    Tt <- cbind(Tt, X %*% W0)
    U <- cbind(U, Y %*% C0)
    X <- X - (X %*% W0) %*% t(W0)
    Y <- Y - (Y %*% C0) %*% t(C0)
    W <- cbind(W, W0)
    C <- cbind(C, C0)
  }
  return(list(W = W, C = C, Tt = Tt, U = U))
}

#' NIPALS method for PLS2
#' 
#' Performs power method for PLS2 loadings.
#' 
#' @inheritParams o2m
#' 
#' @keywords internal
#' @export
pow_o2m <- function(X, Y, n, tol = 1e-10, max_iterations = 100) {
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
#' This is a `slower' implementation of O2PLS, and is using \code{\link{svd}}. For cross-validation purposes, consider using \code{\link{o2m_stripped}}.
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
#' @seealso \code{\link{o2m}}, \code{\link{ssq}}, \code{\link{summary.o2m}}
#' 
#' @keywords internal
#' @export
o2m2 <- function(X, Y, n, nx, ny, stripped = FALSE, tol = 1e-10, max_iterations = 100) {
  
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  if(stripped){
    return(o2m_stripped2(X, Y, n, nx, ny, tol, max_iterations))
  }
  
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
  H_TU <- Tt - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  # R2
  R2Xcorr <- (ssq(Tt %*% t(W))/ssq(X_true))
  R2Ycorr <- (ssq(U %*% t(C))/ssq(Y_true))
  R2X_YO <- (ssq(T_Yosc %*% t(P_Yosc))/ssq(X_true))
  R2Y_XO <- (ssq(U_Xosc %*% t(P_Xosc))/ssq(Y_true))
  R2Xhat <- 1 - (ssq(U %*% B_U %*% t(W) - X_true)/ssq(X_true))
  R2Yhat <- 1 - (ssq(Tt %*% B_T %*% t(C) - Y_true)/ssq(Y_true))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_UT) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, W. = W, U = U, C. = C, E = 0, Ff = 0, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
                U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
                X_hat = 0, Y_hat = 0, R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
                R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  class(model) <- "o2m"
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
#' @seealso \code{\link{ssq}}, \code{\link{o2m}}, \code{\link{loocv}}, \code{\link{adjR2}}
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
  R2Xhat <- 1 - (ssq(U %*% B_U %*% t(W) - X_true) / ssq(X_true))
  R2Yhat <- 1 - (ssq(Tt %*% B_T %*% t(C) - Y_true) / ssq(Y_true))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_TU) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, U = U, W. = W, C. = C, P_Yosc. = P_Yosc, P_Xosc. = P_Xosc,
                T_Yosc. = T_Yosc, U_Xosc. = U_Xosc, W_Yosc = W_Yosc, C_Xosc = C_Xosc,
                B_T. = B_T, B_U = B_U, H_TU = H_TU, H_UT = H_UT, 
                R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, 
                R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  class(model) <- c("o2m","o2m_stripped")
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
#' @seealso \code{\link{ssq}}, \code{\link{o2m}}, \code{\link{loocv}}, \code{\link{adjR2}}
#' @keywords internal
#' @export
o2m_stripped2 <- function(X, Y, n, nx, ny, tol = 1e-10, max_iterations = 100) {
  
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  X_true <- X
  Y_true <- Y
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  T_Yosc <- U_Xosc <- matrix(NA, N, 1)
  P_Yosc <- W_Yosc <- matrix(NA, p, 1)
  P_Xosc <- C_Xosc <- matrix(NA, q, 1)
  
  if (nx + ny > 0) {
    n2 <- n + max(nx, ny)
    
    #    if (N < p & N < q) {
    # When N is smaller than p and q
    W_C <- pow_o2m(X, Y, n2, tol, max_iterations)
    W <- W_C$W
    C <- W_C$C
    Tt <- W_C$Tt
    U <- W_C$U
    rm(W_C)
    gc()
    #    }
    # 3.2.1. 2
    
    if (nx > 0) {
      # 3.2.1. 3
      E_XY <- X - Tt %*% t(W)
      
      udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
      rm(E_XY)
      W_Yosc <- udv$u
      T_Yosc <- X %*% W_Yosc
      P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
      X <- X - T_Yosc %*% t(P_Yosc)
      
      # Update T again (since X has changed) Tt = X%*%W;
    }
    
    # U = Y%*%C; # 3.2.1. 4
    
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
      
      # Update U again (since Y has changed) U = Y%*%C;
    }
  }
  
  # repeat steps 1, 2, and 4 before step 6 When N is smaller than p and q
  #  if (N < p & N < q) {
  W_C <- pow_o2m(X, Y, n, tol, max_iterations)
  W <- W_C$W
  C <- W_C$C
  Tt <- W_C$Tt
  U <- W_C$U
  rm(W_C)
  gc()
  #  }
  
  # 3.2.1. 6
  B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  H_TU <- Tt - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  R2Xcorr <- ssq(Tt) / ssq(X_true)
  R2Ycorr <- ssq(U) / ssq(Y_true)
  R2X_YO <- ssq(T_Yosc %*% t(P_Yosc)) / ssq(X_true)
  R2Y_XO <- ssq(U_Xosc %*% t(P_Xosc)) / ssq(Y_true)
  R2Xhat <- 1 - (ssq(U %*% B_U %*% t(W) - X_true) / ssq(X_true))
  R2Yhat <- 1 - (ssq(Tt %*% B_T %*% t(C) - Y_true) / ssq(Y_true))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_TU) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, U = U, W. = W, C. = C, P_Yosc. = P_Yosc, P_Xosc. = P_Xosc,
                T_Yosc. = T_Yosc, U_Xosc. = U_Xosc, W_Yosc = W_Yosc, C_Xosc = C_Xosc,
                B_T. = B_T, B_U = B_U, H_TU = H_TU, H_UT = H_UT, 
                R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, 
                R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  
  class(model) <- c("o2m","o2m_stripped")
  return(model)
}
