#' Cross-validate procedure for O2PLS
#'
#' @inheritParams o2m
#' @param a Vector of positive integers. Denotes the numbers of joint components to consider.
#' @param ax Vector of non-negative integers. Denotes the numbers of X-specific components to consider.
#' @param ay Vector of non-negative integers. Denotes the numbers of Y-specific components to consider.
#' @param nr_folds Positive integer. Number of folds to consider. Note: \code{kcv=N} gives leave-one-out CV. Note that CV with less than two folds does not make sense.
#' @param nr_cores Positive integer. Number of cores to use for CV. You might want to use \code{\link[parallel:detectCores]{detectCores}()}. Defaults to 1.
#' @param seed Integer. A random seed to make the analysis reproducible.
#'
#' @details This is the standard CV approach. It minimizes the sum of the prediction errors of X and Y over a three-dimensional grid of integers.
#' Parallelization is possible on all platforms. On Windows it uses \code{\link[parallel:makePSOCKcluster]{makePSOCKcluster}}, then exports all necessary objects to the workers,
#'  and then calls \code{\link{parLapply}}. On OSX and Linux the more friendly \code{\link[parallel:mclapply]{mclapply}} is used, which uses forking (but copies your global workspace).
#'  A print method is defined, printing the minimizers and minimum in a readable way. Also the elapsed time is tracked and reported.
#'
#' @return List of class \code{"cvo2m"} with the original and sorted Prediction errors and the number of folds used.
#'
#' @examples
#' local({
#' X = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' Y = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' crossval_o2m(X, Y, a = 1:4, ax = 1:2, ay = 1:2,
#'              nr_folds = 5, nr_cores = 1)
#' })
#'
#' @export
crossval_o2m <- function(X, Y, a, ax, ay, nr_folds, nr_cores = 1,
                         stripped = TRUE, p_thresh = 3000,
                         q_thresh = p_thresh, tol = 1e-10, max_iterations = 100, seed = "off") {
  tic = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_checker(X, Y)
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...\n")}
  kcv = nr_folds
  if(ncol(X) < max(a)+max(ax,ay) | ncol(Y) < max(a)+max(ay,ay))
    message("Some combinations of # components exceed data dimensions, these combinations are not considered\n")
  if(ncol(X) < min(a)+min(ax,ay) | ncol(Y) < min(a)+min(ay,ay))
    stop("There is no valid combination of numbers of components! Please select fewer components in a, ax, ay.\n")
  if(nrow(X) < kcv) stop("There are more folds than samples, please set nr_folds <= ",nrow(X),"\n")
  stopifnot(nr_cores == abs(round(nr_cores)))
  if(nr_folds==1){stop("Cross-validation needs at least two folds, to train and test\n")}
  
  parms = data.frame(nx = ax)
  parms = merge(parms,data.frame(ny = ay))
  parms = merge(parms,data.frame(a = a))
  parms = apply(parms,1,as.list)
  cl_crossval_o2m <- NULL
  
  on.exit({if(!is.null(cl_crossval_o2m)) stopCluster(cl_crossval_o2m)})
  
  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_o2m <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_o2m, library(OmicsPLS))
    clusterExport(cl_crossval_o2m, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_o2m,parms,function(e){
      suppressMessages(loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                                   stripped = stripped, p_thresh = p_thresh, seed = seed, 
                                   q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      suppressMessages(loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                                   stripped = stripped, p_thresh = p_thresh, seed = seed,
                                   q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
    })
  }
  dnams = list(paste("ax=",ax,sep=""),paste("ay=",ay,sep=""),paste("a=",a,sep=""))
  outp1 = array(unlist(outp),dim=c(length(ax),length(ay),length(a)),
                dimnames=dnams)
  outp = aperm(outp1,order(dim(outp1),decreasing=TRUE))
  dnams = dimnames(outp)
  
  if(dim(outp)[3]==1){
    dim(outp) = dim(outp)[-3]
    dimnames(outp) = dnams[1:2]
  }
  
  outpt = list(Original=outp1,Sorted=outp,kcv=kcv)
  class(outpt) <- "cvo2m"
  toc = proc.time() - tic
  outpt$time = round(toc[3],2)
  return(outpt)
}

#' Adjusted Cross-validate procedure for O2PLS
#'
#' Combines CV with R2 optimization
#'
#' @inheritParams crossval_o2m
#'
#' @details This is an alternative way of cross-validating. It is proposed in \code{citation(OmicsPLS)}.
#' This approach is (much) faster than the standard \code{crossval_o2m} approach and works fine even with two folds.
#' For each element in \code{n} it looks for nx and ny that maximize the \eqn{R^2} between T and U in the O2PLS model.
#' This approach often yields similar integer as the standard approach.
#' We however suggest to use the standard approach to minimize the prediction error around the found integers.
#'
#' @return data.frame with four columns: MSE, n, nx and ny. Each row corresponds to an element in \code{a}.
#' @examples
#' local({
#' X = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' Y = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' crossval_o2m_adjR2(X, Y, a = 1:4, ax = 1:2, ay = 1:2,
#'              nr_folds = 5, nr_cores = 1)
#' })
#' @export
crossval_o2m_adjR2 <- function(X, Y, a, ax, ay, nr_folds, nr_cores = 1,
                               stripped = TRUE, p_thresh = 3000, seed = "off", 
                               q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  tic = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_checker(X, Y)
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...\n")}
  kcv = nr_folds
  if(ncol(X) < max(a)+max(ax,ay) | ncol(Y) < max(a)+max(ay,ay))
    message("Some combinations of # components exceed data dimensions, these combinations are not considered\n")
  if(ncol(X) < min(a)+min(ax,ay) | ncol(Y) < min(a)+min(ay,ay))
    stop("There is no valid combination of numbers of components! Please select fewer components in a, ax, ay.\n")
  if(nrow(X) < kcv) stop("There are more folds than samples, please set nr_folds <= ",nrow(X),"\n")
  stopifnot(nr_cores == abs(round(nr_cores)))
  if(nr_folds==1){stop("Cross-validation needs at least two folds, to train and test\n")}
  
  cl_crossval_o2m <- NULL
  on.exit({if(!is.null(cl_crossval_o2m)) stopCluster(cl_crossval_o2m)})
  
  parms = data.frame(a = a)
  parms = apply(parms,1,as.list)
  
  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_o2m <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_o2m, library(OmicsPLS))
    clusterExport(cl_crossval_o2m, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_o2m,parms,function(e){
      parms = data.frame(nx = ax)
      parms = merge(parms,data.frame(ny = ay))
      parms = apply(parms,1,as.list)
      R2grid = matrix(colMeans(suppressMessages(adjR2(Y, X, e$a, ax, ay,
                                                      stripped = stripped, p_thresh = p_thresh,
                                                      q_thresh = q_thresh, tol = tol, max_iterations = max_iterations))),
                      nrow = length(ay), byrow=TRUE)
      nxny = which(R2grid == max(R2grid), arr.ind = TRUE)[1,]
      a_mse = suppressMessages(loocv_combi(X,Y,e$a,ax[nxny[2]],ay[nxny[1]],app_err=F,func=o2m,kcv=kcv,
                                           stripped = stripped, p_thresh = p_thresh, seed = seed, 
                                           q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
      c(a_mse, e$a, ax[nxny[2]],ay[nxny[1]])
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      parms = data.frame(nx = ax)
      parms = merge(parms,data.frame(ny = ay))
      parms = apply(parms,1,as.list)
      R2grid = matrix(colMeans(suppressMessages(adjR2(Y, X, e$a, ax, ay,
                                                      stripped = stripped, p_thresh = p_thresh,
                                                      q_thresh = q_thresh, tol = tol, max_iterations = max_iterations))),
                      nrow = length(ay), byrow=TRUE)
      #R2grid[which(is.na(R2grid))] = -999
      nxny = which(R2grid == max(R2grid,na.rm = TRUE), arr.ind = TRUE)[1,]
      a_mse = suppressMessages(loocv_combi(X,Y,e$a,ax[nxny[2]],ay[nxny[1]],app_err=F,func=o2m,kcv=kcv,
                                           stripped = stripped, p_thresh = p_thresh, seed = seed, 
                                           q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
      c(a_mse, e$a, ax[nxny[2]],ay[nxny[1]])
    })
  }
  outp2 = matrix(unlist(outp), nrow = length(a), byrow = T)
  outp2 <- as.data.frame(outp2)
  names(outp2) <- c("MSE", "n", "nx", "ny")
  message("Minimum is at n = ", outp2[,2][which.min(outp2[,1])], sep = ' ')
  message("Elapsed time: ", round((proc.time() - tic)[3],2), " sec")
  return(outp2)
}

#' Cross-validate procedure for O2PLS
#'
#' @param x List of class \code{"cvo2m"}, produced by \code{\link{crossval_o2m}}.
#' @param include_matrix Logical. Should the 3-d array with Prediction errors also be printed.
#' @param ... For consistency.
#'
#' @return NULL
#' @export
print.cvo2m <- function(x,include_matrix=FALSE,...) {
  wmCV = which(min(x$Or,na.rm = TRUE)==x$Or,arr.ind = TRUE,useNames = FALSE)
  dnams = dimnames(x$Or)
  dnams1 = dnams[[1]][wmCV[1]]
  dnams2 = dnams[[2]][wmCV[2]]
  dnams3 = dnams[[3]][wmCV[3]]
  
  cat("*******************\n")
  cat("Elapsed time: ",x$time, " sec", '\n', sep='')
  cat("*******\n")
  cat("Minimal ",x$kcv,"-CV error is at ",dnams1," ",dnams2," ",dnams3," ","\n",sep="")
  cat("*******\n")
  cat("Minimum MSE is",min(x$Sor,na.rm = TRUE),"\n")
  if(include_matrix){
    cat("*******\n")
    cat("Simplified CV matrix is \n")
    print(x$Sorted)
  }
  cat("*******************\n")
}


###### Penalized part ######
#' Perform cross-validation to find the optimal number of variables/groups to keep for each joint component
#'
#' @inheritParams o2m
#' @param nr_folds Integer. Number of folds of CV
#' @param keepx_seq Numeric vector. A vector indicating how many variables/groups to keep for CV in each of the joint component of X. Sparsity of each joint component will be selected sequentially.
#' @param keepy_seq Numeric vector. A vector indicating how many variables/groups to keep for CV in each of the joint component of Y. Sparsity of each joint component will be selected sequentially.
#' @return A list containing
#'    \item{x_1sd}{A vector with length n, giving the optimal number of variables/groups to keep for each X-joint compoent. One standard error rule is applied}
#'    \item{y_1sd}{A vector with length n, giving the optimal number of variables/groups to keep for each Y-joint compoent. One standard error rule is applied}
#'    \item{x}{A vector with length n, giving the optimal number of variables/groups to keep for each X-joint compoent, without applying the one standard error rule}
#'    \item{y}{A vector with length n, giving the optimal number of variables/groups to keep for each Y-joint compoent, without applying the one standard error rule}
#' @export
crossval_sparsity <- function(X, Y, n, nx, ny, nr_folds, keepx_seq=NULL, keepy_seq=NULL, groupx=NULL, groupy=NULL, tol = 1e-10, max_iterations = 100){
  
  if(is.null(groupx) & is.null(groupy)){
    method = "SO2PLS"
    message("Group information not provided, CV for number of variables to keep\n")
    cv_lambda_checker(X, Y, keepx_seq, keepy_seq)
    if(is.null(keepx_seq)) keepx_seq <- ncol(X)
    if(is.null(keepy_seq)) keepy_seq <- ncol(Y)
  }else{
    method = "GO2PLS"
    message("Group information provided, CV for number of groups to keep\n")
    # check if only information for one dataset is provided
    if(is.null(groupx))  groupx = colnames(X)
    if(is.null(groupy))  groupy = colnames(Y)
    cv_lambda_checker_group(groupx, groupy, keepx_seq, keepy_seq)
    if(is.null(keepx_seq)) keepx_seq <- ncol(X)
    if(is.null(keepy_seq)) keepy_seq <- ncol(Y)
  }
  
  # Check format
  stopifnot(all(n == round(n)), nr_folds == round(nr_folds))
  X = as.matrix(X)
  Y = as.matrix(Y)
  input_checker(X, Y)
  
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
  # Initiating variables
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  
  # initiate
  mean_covTU <- srr_covTU <- matrix(NA, nrow = length(keepy_seq), ncol = length(keepx_seq))
  
  rownames(mean_covTU) <- rownames(srr_covTU) <- keepy_seq
  colnames(mean_covTU) <- colnames(srr_covTU) <- keepx_seq
  
  covTU <- NA * 1:nr_folds
  keepxy_x <- keepxy_y <- x_max <- y_max <- vector()
  
  if(method == "SO2PLS"){
    for (comp in 1:n) {
      kx <- 0
      
      # Creating blocks and folds
      blocks <- cut(seq(1:N), breaks=nr_folds, labels=F)
      folds <- sample(N)
      
      # Loop through a grid of n_lambda * n_lambda
      for (lx in keepx_seq) {
        kx <- kx +1
        ky <- 0
        for (ly in keepy_seq) {
          ky <- ky + 1
          
          # loop through number of folds
          for (i in 1:nr_folds) {
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
            #covTU[i] <- drop(cov(t_tst, u_tst))
            covTU[i] <- abs(drop(cov(t_tst, u_tst)))
          }
          
          # Test cov
          mean_covTU[ky,kx] <- mean(covTU)
          srr_covTU[ky,kx] <- sd(covTU)/sqrt(nr_folds)
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
      #solve(t(0))
      keepxy_x[comp] <- keepxy$x
      keepxy_y[comp] <- keepxy$y
      y_max[comp] <- as.numeric(rownames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[1]])
      x_max[comp] <- as.numeric(colnames(mean_covTU)[which(mean_covTU == max(mean_covTU), arr.ind = T)[2]])
    }
  }else{
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
      blocks <- cut(seq(1:N), breaks=nr_folds, labels=F)
      folds <- sample(N)
      
      # Loop through a grid of n_lambda * n_lambda
      for (lx in keepx_seq) {
        kx <- kx +1
        ky <- 0
        for (ly in keepy_seq) {
          ky <- ky + 1
          
          # loop through number of folds
          for (i in 1:nr_folds) {
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
            covTU[i] <- abs(drop(cov(t_tst, u_tst)))
          }
          
          # Test cov
          mean_covTU[ky,kx] <- mean(covTU)
          srr_covTU[ky,kx] <- sd(covTU)/sqrt(nr_folds)
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
  }
  
  # Output Change here the standard
  bestsp <- list()
  bestsp$x_1sd <- keepxy_x
  bestsp$y_1sd <- keepxy_y
  #bestsp$err_tu <- mean_covTU
  #bestsp$srr <- srr_covTU
  bestsp$x <- x_max
  bestsp$y <- y_max
  return(list(Best = unlist(bestsp), Covs = mean_covTU, SEcov = srr_covTU))
}


#' Internal function for crossval_sparsity
#'
#' @param dat Matrix with numeric row/col names
#' @param index Get from which(..., arr.ind = T)
#' @param p Number of variables in X
#' @param q Number of variables in Y
#'
#' @details This function finds the most sparse combination of keepx and keepy (min(keepx/p + keepy/q)) that yields cov(T,U) within 1 std error of the largest cov(T,U). Note that it's possible that the resulting keepx or keepy is larger than the orignal when p >> q or p << q.
#' @keywords internal
#' @export
err_back <- function(dat, index, p, q){
  index <- index %>% tibble::as_tibble() %>% dplyr::mutate(sp = as.numeric(rownames(dat)[row])/q + as.numeric(colnames(dat)[col])/p)
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
