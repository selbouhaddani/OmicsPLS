#' @export
# Different sparsity for comps
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
    # 1 stardard err back
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
  bestsp$xmax <- x_max 
  bestsp$ymax <- y_max
  return(bestsp)
}

######################################################
# Same sparsity across all comps
best_keepxy_srr_old <- function(X, Y, n, nx, ny, lambda_kcv, keepx_seq, keepy_seq, tol = 1e-10, max_iterations = 100){
  
  # Check format
  stopifnot(all(n == round(n)), lambda_kcv == round(lambda_kcv))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
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
        
        pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                     stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                     sparsity_it = T, keepx = lx, keepy = ly, max_iterations_sparsity = max_iterations)
        
        fit <- suppressMessages(do.call(o2m, pars))
        
        # Test error
        if(inherits(fit, 'try-error')){
          covTU[i] <- NA
        }else{
          sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
          covTU[i] <- sum_R2$cov_tu
        }
      }
      mean_covTU[ky,kx] <- mean(covTU)
      srr_covTU[ky,kx] <- sd(covTU)/sqrt(lambda_kcv)
    }
  }
  
  
  # 1 stardard err back
  cov_max <- max(mean_covTU)
  cov_1srr <- cov_max - srr_covTU[which.max(mean_covTU)]
  keepxy <- err_back(mean_covTU, which(mean_covTU > cov_1srr, arr.ind = T), dim(X)[2], dim(Y)[2])
  
  
  # Output Change here the standard
  bestsp <- list()
  bestsp$y <- keepxy$y
  bestsp$x <- keepxy$x
  
  bestsp$err_tu <- mean_covTU
  bestsp$srr <- srr_covTU
  
  bestsp$ymax <- as.numeric(rownames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[1]])
  bestsp$xmax <- as.numeric(colnames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[2]])
  
  return(bestsp)
}


best_keepxy_errtu_srr <- function(X, Y, n, nx, ny, lambda_kcv, keepx_seq, keepy_seq, tol = 1e-10, max_iterations = 100){
  
  # Check format
  stopifnot(all(n == round(n)), lambda_kcv == round(lambda_kcv))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
  # Initiating variables
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  
  # initiate
  mean_err_tu <- srr_tu <- matrix(NA, nrow = length(keepy_seq), ncol = length(keepx_seq))
  
  rownames(mean_err_tu) <- rownames(srr_tu) <- keepy_seq
  colnames(mean_err_tu) <- colnames(srr_tu) <- keepx_seq
  
  err_t <- NA * 1: lambda_kcv
  err_u <- NA * 1: lambda_kcv
  
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
        
        pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                     stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                     sparsity_it = T, keepx = lx, keepy = ly, max_iterations_sparsity = max_iterations)
        
        fit <- suppressMessages(do.call(o2m, pars))
        
        # Test error
        if(inherits(fit, 'try-error')){
          err_t[i] <- NA
          err_u[i] <- NA
        }else{
          sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
          err_t[i] <- sum_R2$SSE_t
          err_u[i] <- sum_R2$SSE_u
        }
      }
      err_tu <- err_t + err_u
      mean_err_tu[ky,kx] <- mean(err_tu)
      srr_tu[ky,kx] <- sd(err_tu)/sqrt(lambda_kcv)
    }
  }
  
  # 1 stardard err back
  err_min <- min(mean_err_tu)
  err_1srr <- err_min + srr_tu[which.min(mean_err_tu)]
  keepxy <- err_back(mean_err_tu, which(mean_err_tu < err_1srr, arr.ind = T), dim(X)[2], dim(Y)[2])
  

  # Output Change here the standard
  bestsp <- list()
  bestsp$y <- keepxy$y
  bestsp$x <- keepxy$x
  
  bestsp$err_tu <- mean_err_tu
  bestsp$srr <- srr_tu
  
  return(bestsp)
}

# dat is matrix with numeric row/col names, index get from which(..., arr.ind = T)
err_back <- function(dat, index, p, q){
  index <- index %>% as_tibble() %>% mutate(sp = as.numeric(rownames(dat)[row])/q + as.numeric(colnames(dat)[col])/p)
  # find most sparse model
  temp <- which(index$sp == min(index$sp))
  # if draw
  final <- vector()
  for(i in 1:length(temp)){
    final[i] <- dat[index$row[temp[i]], index$col[temp[i]]]
  }
  # which.max if to maximize, which.min if to minimize
  final <- temp[which.max(final)]
  return(list(x = as.numeric(colnames(dat)[index$col[final]]), 
              y = as.numeric(rownames(dat)[index$row[final]])))
}



## Two rounds one-dimensional search (min + 1standard err)
best_keepxy_errtu_srr1d <- function(X, Y, n, nx, ny, lambda_kcv, keepx_seq, keepy_seq, tol = 1e-10, max_iterations = 100){
  
  # Check format
  stopifnot(all(n == round(n)), lambda_kcv == round(lambda_kcv))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
  # Initiating variables
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  
  # initiate
  mean_err_tu_x1 <- mean_err_tu_x2 <- vector()
  mean_err_tu_y1 <- mean_err_tu_y2 <- vector()
  srr_tu_x1 <- srr_tu_x2 <- srr_tu_y1 <- srr_tu_y2 <- vector()
  
  err_t <- NA * 1: lambda_kcv
  err_u <- NA * 1: lambda_kcv

  # Creating blocks and folds
  blocks <- cut(seq(1:N), breaks=lambda_kcv, labels=F)
  folds <- sample(N)
  
  # First find keepx
  kx <- 0
  for (lx in keepx_seq) {
    kx <- kx +1
      # loop through number of folds
      for (i in 1:lambda_kcv) {
        ii <- which(blocks==i)
        
        pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                     stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                     sparsity_it = T, keepx = lx, keepy = tail(keepy_seq, 1), max_iterations_sparsity = max_iterations)
        
        fit <- suppressMessages(do.call(o2m, pars))
        
        # Test error
        if(inherits(fit, 'try-error')){
          err_t[i] <- NA
          err_u[i] <- NA
        }else{
          sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
          err_t[i] <- sum_R2$SSE_t
          err_u[i] <- sum_R2$SSE_u
        }
      }
      mean_err_tu_x1[kx] <- mean(err_t + err_u)
      srr_tu_x1[kx] <- sd(err_t + err_u)/sqrt(lambda_kcv)
      err_min <- min(mean_err_tu_x1)
      err_1srr <- err_min + srr_tu_x1[which.min(mean_err_tu_x1)]
      keepx1 <- keepx_seq[which.min(mean_err_tu_x1)]
      #keepx1 <- keepx_seq[which(mean_err_tu_x < err_1srr)[1]]
  }
  
  # First-round find keepy
  ky <- 0
  for (ly in keepy_seq) {
    ky <- ky +1
    # loop through number of folds
    for (i in 1:lambda_kcv) {
      ii <- which(blocks==i)
      
      pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                   stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                   sparsity_it = T, keepx = keepx1, keepy = ly, max_iterations_sparsity = max_iterations)
      
      fit <- suppressMessages(do.call(o2m, pars))
      
      # Test error
      if(inherits(fit, 'try-error')){
        err_t[i] <- NA
        err_u[i] <- NA
      }else{
        sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
        err_t[i] <- sum_R2$SSE_t
        err_u[i] <- sum_R2$SSE_u
      }
    }
    mean_err_tu_y1[ky] <- mean(err_t + err_u)
    srr_tu_y1[ky] <- sd(err_t + err_u)/sqrt(lambda_kcv)
    err_min <- min(mean_err_tu_y1)
    err_1srr <- err_min + srr_tu_y1[which.min(mean_err_tu_y1)]
    keepy1 <- keepy_seq[which.min(mean_err_tu_y1)]
    #keepy1 <- keepy_seq[which(mean_err_tu_y < err_1srr)[1]]
  }
  
  # Second-round find keepx
  kx <- 0
  for (lx in keepx_seq) {
    kx <- kx +1
    # loop through number of folds
    for (i in 1:lambda_kcv) {
      ii <- which(blocks==i)
      
      pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                   stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                   sparsity_it = T, keepx = lx, keepy = keepy1, max_iterations_sparsity = max_iterations)
      
      fit <- suppressMessages(do.call(o2m, pars))
      
      # Test error
      if(inherits(fit, 'try-error')){
        err_t[i] <- NA
        err_u[i] <- NA
      }else{
        sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
        err_t[i] <- sum_R2$SSE_t
        err_u[i] <- sum_R2$SSE_u
      }
    }
    mean_err_tu_x2[kx] <- mean(err_t + err_u)
    srr_tu_x2[kx] <- sd(err_t + err_u)/sqrt(lambda_kcv)
    err_min <- min(mean_err_tu_x2)
    err_1srr <- err_min + srr_tu_x2[which.min(mean_err_tu_x2)]
    keepx2 <- keepx_seq[which(mean_err_tu_x2 < err_1srr)[1]]
  }
  
  # Second-round find keepy
  ky <- 0
  for (ly in keepy_seq) {
    ky <- ky +1
    # loop through number of folds
    for (i in 1:lambda_kcv) {
      ii <- which(blocks==i)
      
      pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                   stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                   sparsity_it = T, keepx = keepx2, keepy = ly, max_iterations_sparsity = max_iterations)
      
      fit <- suppressMessages(do.call(o2m, pars))
      
      # Test error
      if(inherits(fit, 'try-error')){
        err_t[i] <- NA
        err_u[i] <- NA
      }else{
        sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
        err_t[i] <- sum_R2$SSE_t
        err_u[i] <- sum_R2$SSE_u
      }
    }
    mean_err_tu_y2[ky] <- mean(err_t + err_u)
    srr_tu_y2[ky] <- sd(err_t + err_u)/sqrt(lambda_kcv)
    err_min <- min(mean_err_tu_y2)
    err_1srr <- err_min + srr_tu_y2[which.min(mean_err_tu_y2)]
    keepy2 <- keepy_seq[which(mean_err_tu_y2 < err_1srr)[1]]
  }

  
  # Output Change here the standard
  bestsp <- list()
  bestsp$y <- keepy2
  bestsp$x <- keepx2
  bestsp$xy1 <- c(keepx1, keepy1)
  bestsp$err_tu_x <- rbind(mean_err_tu_x1, mean_err_tu_x2)
  bestsp$err_tu_y <- rbind(mean_err_tu_y1, mean_err_tu_y2)
  bestsp$srr_tu_x <- rbind(srr_tu_x1, srr_tu_x2)
  bestsp$srr_tu_y <- rbind(srr_tu_y1, srr_tu_y2)
  
  return(bestsp)
}



best_keepxy_test <- function(X, Y, n, nx, ny, lambda_kcv, keepx_seq, keepy_seq, tol = 1e-10, max_iterations = 100){
  
  # Check format
  stopifnot(all(n == round(n)), lambda_kcv == round(lambda_kcv))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
  # Initiating variables
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  
  # initiate
  mean_err_t <- mean_err_u <- mean_err_tu <- matrix(NA, nrow = length(keepy_seq), ncol = length(keepx_seq))
  
  rownames(mean_err_t) <- rownames(mean_err_u) <- rownames(mean_err_tu) <- keepy_seq
  colnames(mean_err_t) <- colnames(mean_err_u) <- colnames(mean_err_tu) <- keepx_seq

  err_t <- NA * 1: lambda_kcv
  err_u <- NA * 1: lambda_kcv
  
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
        
        pars <- list(X = X[-folds[ii], ], Y = Y[-folds[ii], ], n = n, nx = nx, ny = ny,
                     stripped = FALSE, tol = tol, p_thresh = 1, max_iterations = max_iterations,
                     sparsity_it = T, keepx = lx, keepy = ly, max_iterations_sparsity = max_iterations)
        
        fit <- suppressMessages(do.call(o2m, pars))
        
        # Test error
        if(inherits(fit, 'try-error')){
          err_t[i] <- NA
          err_u[i] <- NA
        }else{
          sum_R2 <- sumR2_combi_test(X[folds[ii], ], Y[folds[ii], ], fit)
          err_t[i] <- sum_R2$SSE_t
          err_u[i] <- sum_R2$SSE_u
        }
      }
      mean_err_t[ky,kx] <- mean(err_t)
      mean_err_u[ky,kx] <- mean(err_u)
    }
  }
  mean_err_tu <- mean_err_t + mean_err_u
  
  # Output Change here the standard
  bestsp <- list()
  bestsp$y <- as.numeric(rownames(mean_err_tu)[which(mean_err_tu == min(mean_err_tu), arr.ind = T)[1]])
  bestsp$x <- as.numeric(colnames(mean_err_tu)[which(mean_err_tu == min(mean_err_tu), arr.ind = T)[2]])

  bestsp$err_tu <- mean_err_tu

  return(bestsp)
}



sumR2_combi_test <- function(Xtst, Ytst, fit){
  if(!inherits(fit,"o2m")) stop("fit should be an O2PLS fit")
  
  if (!is.matrix(Xtst)) Xtst <- t(Xtst)
  
  if (!is.matrix(Ytst)) Ytst <- t(Ytst)
  
  input_checker(Xtst, Ytst)
  
  # Orthogonal correction
  Xtst <- Xtst - Xtst %*% fit$W_Yosc %*% t(fit$P_Yosc)
  Ytst <- Ytst - Ytst %*% fit$C_Xosc %*% t(fit$P_Xosc)

  T_tst <- Xtst %*% fit$W.
  U_tst <- Ytst %*% fit$C.
  
  # predict scores t and u
  SSE_u <- ssq(U_tst - T_tst %*% fit$B_T)/ssq(U_tst)
  SSE_t <- ssq(T_tst - Ytst %*% fit$C. %*% fit$B_U)/ssq(Xtst %*% fit$W.)
  cov_tu <- sum(diag(cov(T_tst, U_tst)))
  
  sum_R2 <- list()
  sum_R2$SSE_u <- SSE_u
  sum_R2$SSE_t <- SSE_t
  sum_R2$cov_tu <- cov_tu
  return(sum_R2)
}
