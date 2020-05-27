#' Perform cross-validation to find the optimal number of groups to keep for each joint component
#'
#' @param X Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param Y Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param n Integer. Number of joint PLS components. Must be positive!
#' @param nx Integer. Number of orthogonal components in \eqn{X}. Negative values are interpreted as 0
#' @param ny Integer. Number of orthogonal components in \eqn{Y}. Negative values are interpreted as 0
#' @param lambda_kcv Integer. Number of folds of CV
#' @param groupx Character. A vecter or character indicating group names of the variables. The order of group names must correspond to the order of the vairables in X.
#' @param groupy Character. A vecter or character indicating group names of the variables. The order of group names must correspond to the order of the vairables in Y.
#' @param keepx_seq Vector. A vector indicating how many groups to keep for CV in each of the joint component of X.
#' @param keepy_gr Vector. A vector indicating how many groups to keep for CV in each of the joint component of Y.
#' @param tol double. Threshold for power method iteration
#' @param max_iterations Integer, Maximum number of iterations for power method
#' @return A list containing
#'    \item{x}{A vector with length n, giving the optimal number of groups to keep for each X-joint compoent. One standard error rule is applied}
#'    \item{y}{A vector with length n, giving the optimal number of groups to keep for each Y-joint compoent. One standard error rule is applied}
#'    \item{x_max}{A vector with length n, giving the optimal number of groups to keep for each X-joint compoent, without applying the one standard error rule}
#'    \item{y_max}{A vector with length n, giving the optimal number of groups to keep for each Y-joint compoent, without applying the one standard error rule}
#' @export
best_keepgr_srr <- function(X, Y, n, nx, ny, lambda_kcv, groupx, groupy, keepx_seq, keepy_seq, tol = 1e-10, max_iterations = 100){
  
  # Check format
  stopifnot(all(n == round(n)), lambda_kcv == round(lambda_kcv))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
  #######################
  # Filter O2
  if (nx + ny > 0) {
    # larger principal subspace
    n2 <- n + max(nx, ny)
    
    # if(N<p&N<q){ # When N is smaller than p and q
    W_C <- suppressMessages(pow_o2m(X, Y, n2, tol, max_iterations))
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
  ##########################
  # Initiating variables
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  
  # initiate
  mean_covTU <- srr_covTU <- matrix(NA, nrow = length(keepy_seq), ncol = length(keepx_seq))
  
  rownames(mean_covTU) <- rownames(srr_covTU) <- keepy_seq
  colnames(mean_covTU) <- colnames(srr_covTU) <- keepx_seq
  
  covTU <- NA * 1: lambda_kcv
  keepxy_x <- keepxy_y <- x_max <- y_max <- vector()
  
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
  
  
  for (comp in 1:n) {
    kx <- 0
    
    # Creating blocks and folds
    blocks <- cut(seq(1:N), breaks=lambda_kcv, labels=F)
    folds <- sample(N)
    
    # Loop through a grid of n_lambda * n_lambda
    for (lx in keepx_seq) {
      kx <- kx +1
      ky <- 0
      for (ly in keepy_seq) {
        ky <- ky + 1
        
        # loop through number of folds
        for (i in 1:lambda_kcv) {
          ii <- which(blocks==i)
          X_tr <- X[-folds[ii], ]
          X_tst <- X[folds[ii], ]
          Y_tr <- Y[-folds[ii], ]
          Y_tst <- Y[folds[ii], ]
          
          v <- X_tr[1,]/norm_vec(X_tr[1,])
          for (k in 1: max_iterations){
            v_old <- v
            u <- t(Y_tr) %*% (X_tr %*% v)
            ul <- thresh_n_gr(u, ly, index_gry)
            u <- ul$w
            u <- u/norm_vec(u)
            v <- t(X_tr) %*% (Y_tr %*% u)
            vl <- thresh_n_gr(v, lx, index_grx)
            v <- vl$w
            v <- v/norm_vec(v)
            if (mse(v, v_old) < tol) {
              break
            }
          }

          t_tst <- X_tst %*% v
          u_tst <- Y_tst %*% u
          covTU[i] <- drop(cov(t_tst, u_tst))
        }
        
        # Test cov
        mean_covTU[ky,kx] <- mean(covTU)
        srr_covTU[ky,kx] <- sd(covTU)/sqrt(lambda_kcv)
      }
    }
    # 1 stardard err rule
    cov_max <- max(mean_covTU)
    cov_1srr <- cov_max - srr_covTU[which.max(mean_covTU)]
    keepxy <- err_back(mean_covTU, which(mean_covTU > cov_1srr, arr.ind = T), nr_grx, nr_gry)
    
    v <- X[1,]/norm_vec(X[1,])
    for (k in 1: max_iterations){
      v_old <- v
      u <- t(Y) %*% (X %*% v)
      ul <- thresh_n_gr(u, keepxy$y, index_gry)
      u <- ul$w
      u <- u/norm_vec(u)
      v <- t(X) %*% (Y %*% u)
      vl <- thresh_n_gr(v, keepxy$x, index_grx)
      v <- vl$w
      v <- v/norm_vec(v)
      if (mse(v, v_old) < tol) {
        break
      }
    }
    
    t_tmp <- X %*% v
    u_tmp <- Y %*% u
    
    p <- t(X) %*% t_tmp / drop(crossprod(t_tmp))
    q <- t(Y) %*% u_tmp / drop(crossprod(u_tmp))
    X <- X - t_tmp %*% t(p)
    Y <- Y - u_tmp %*% t(q)
    
    keepxy_x[comp] <- keepxy$x
    keepxy_y[comp] <- keepxy$y
    y_max[comp] <- as.numeric(rownames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[1]])
    x_max[comp] <- as.numeric(colnames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[2]])
  }
  
  # Output Change here the standard
  bestsp <- list()
  bestsp$x <- keepxy_x
  bestsp$y <- keepxy_y
  # bestsp$err_tu <- mean_covTU
  # bestsp$srr <- srr_covTU
  bestsp$xmax <- x_max 
  bestsp$ymax <- y_max
  return(bestsp)
}


#' Perform cross-validation to find the optimal number of variables to keep for each joint component
#'
#' @param X Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param Y Numeric matrix. Vectors will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param n Integer. Number of joint PLS components. Must be positive!
#' @param nx Integer. Number of orthogonal components in \eqn{X}. Negative values are interpreted as 0
#' @param ny Integer. Number of orthogonal components in \eqn{Y}. Negative values are interpreted as 0
#' @param lambda_kcv Integer. Number of folds of CV
#' @param keepx_seq Vector. A vector indicating how many variables to keep for CV in each of the joint component of X.
#' @param keepy_gr Vector. A vector indicating how many variables to keep for CV in each of the joint component of Y.
#' @param tol double. Threshold for power method iteration
#' @param max_iterations Integer, Maximum number of iterations for power method
#' @return A list containing
#'    \item{x}{A vector with length n, giving the optimal number of variables to keep for each X-joint compoent. One standard error rule is applied}
#'    \item{y}{A vector with length n, giving the optimal number of variables to keep for each Y-joint compoent. One standard error rule is applied}
#'    \item{x_max}{A vector with length n, giving the optimal number of variables to keep for each X-joint compoent, without applying the one standard error rule}
#'    \item{y_max}{A vector with length n, giving the optimal number of variables to keep for each Y-joint compoent, without applying the one standard error rule}
#' @export
best_keepxy_srr <- function(X, Y, n, nx, ny, lambda_kcv, keepx_seq, keepy_seq, tol = 1e-10, max_iterations = 100){
  
  # Check format
  stopifnot(all(n == round(n)), lambda_kcv == round(lambda_kcv))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
  #######################
  # Filter O2
  if (nx + ny > 0) {
    # larger principal subspace
    n2 <- n + max(nx, ny)
    
    # if(N<p&N<q){ # When N is smaller than p and q
    W_C <- suppressMessages(pow_o2m(X, Y, n2, tol, max_iterations))
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
  ##########################
  # Initiating variables
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  
  # initiate
  mean_covTU <- srr_covTU <- matrix(NA, nrow = length(keepy_seq), ncol = length(keepx_seq))
  
  rownames(mean_covTU) <- rownames(srr_covTU) <- keepy_seq
  colnames(mean_covTU) <- colnames(srr_covTU) <- keepx_seq
  
  covTU <- NA * 1: lambda_kcv
  keepxy_x <- keepxy_y <- x_max <- y_max <- vector()
  
  for (comp in 1:n) {
    kx <- 0
    
    # Creating blocks and folds
    blocks <- cut(seq(1:N), breaks=lambda_kcv, labels=F)
    folds <- sample(N)
    
    # Loop through a grid of n_lambda * n_lambda
    for (lx in keepx_seq) {
      kx <- kx +1
      ky <- 0
      for (ly in keepy_seq) {
        ky <- ky + 1
        
        # loop through number of folds
        for (i in 1:lambda_kcv) {
          ii <- which(blocks==i)
          X_tr <- X[-folds[ii], ]
          X_tst <- X[folds[ii], ]
          Y_tr <- Y[-folds[ii], ]
          Y_tst <- Y[folds[ii], ]
          
            v <- X_tr[1,]/norm_vec(X_tr[1,])
            for (k in 1: max_iterations){
              v_old <- v
              u <- t(Y_tr) %*% (X_tr %*% v)
              u <- thresh_n(u, ly)
              u <- u/norm_vec(u)
              v <- t(X_tr) %*% (Y_tr %*% u)
              v <- thresh_n(v, lx)
              v <- v/norm_vec(v)
              if (mse(v, v_old) < tol) {
                break
              }
            }
          
          t_tst <- X_tst %*% v
          u_tst <- Y_tst %*% u
          covTU[i] <- drop(cov(t_tst, u_tst))
        }
        
        # Test cov
        mean_covTU[ky,kx] <- mean(covTU)
        srr_covTU[ky,kx] <- sd(covTU)/sqrt(lambda_kcv)
      }
    }
    # 1 stardard err rule
    cov_max <- max(mean_covTU)
    cov_1srr <- cov_max - srr_covTU[which.max(mean_covTU)]
    keepxy <- err_back(mean_covTU, which(mean_covTU > cov_1srr, arr.ind = T), dim(X)[2], dim(Y)[2])
    
    v <- X[1,]/norm_vec(X[1,])
    for (k in 1: max_iterations){
      v_old <- v
      u <- t(Y) %*% (X %*% v)
      u <- thresh_n(u, keepxy$y)
      u <- u/norm_vec(u)
      v <- t(X) %*% (Y %*% u)
      v <- thresh_n(v, keepxy$x)
      v <- v/norm_vec(v)
      if (mse(v, v_old) < tol) {
        break
      }
    }
    t_tmp <- X %*% v
    u_tmp <- Y %*% u
    
    p <- t(X) %*% t_tmp / drop(crossprod(t_tmp))
    q <- t(Y) %*% u_tmp / drop(crossprod(u_tmp))
    X <- X - t_tmp %*% t(p)
    Y <- Y - u_tmp %*% t(q)
    
    keepxy_x[comp] <- keepxy$x
    keepxy_y[comp] <- keepxy$y
    y_max[comp] <- as.numeric(rownames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[1]])
    x_max[comp] <- as.numeric(colnames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[2]])
  }
  
  # Output Change here the standard
  bestsp <- list()
  bestsp$x <- keepxy_x
  bestsp$y <- keepxy_y
  # bestsp$err_tu <- mean_covTU
  # bestsp$srr <- srr_covTU
  bestsp$x_max <- x_max 
  bestsp$y_max <- y_max
  return(bestsp)
}
