#' Cross-validate procedure for O2PLS
#'
#' @inheritParams o2m
#' @param a Vector of positive integers. Denotes the numbers of joint components to consider. 
#' @param ax Vector of non-negative integers. Denotes the numbers of X-specific components to consider.
#' @param ay Vector of non-negative integers. Denotes the numbers of Y-specific components to consider.
#' @param nr_folds Positive integer. Number of folds to consider. Note: \code{kcv=N} gives leave-one-out CV. Note that CV with less than two folds does not make sense.
#' @param nr_cores Positive integer. Number of cores to use for CV. You might want to use \code{\link{detectCores}()}. Defaults to 1.
#' 
#' @details This is the standard CV approach. It minimizes the sum of the prediction errors of X and Y over a three-dimensional grid of integers.
#' Parallelization is possible on all platforms. On Windows it uses \code{\link{makePSOCKcluster}}, then exports all necessary objects to the workers,
#'  and then calls \code{\link{parLapply}}. On OSX and Linux the more friendly \code{\link{mclapply}} is used, which uses forking.
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
                         q_thresh = p_thresh, tol = 1e-10, max_iterations = 100) {
  tic = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_checker(X, Y)
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...\n")}
  kcv = nr_folds
  if(ncol(X) < max(a)+max(ax,ay) | ncol(Y) < max(a)+max(ay,ay)) 
    warning("Some combinations of # components exceed data dimensions, these combinations are not considered\n")
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
                  stripped = stripped, p_thresh = p_thresh, 
                  q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      suppressMessages(loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                  stripped = stripped, p_thresh = p_thresh, 
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
                               stripped = TRUE, p_thresh = 3000, 
                               q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  tic = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_checker(X, Y)
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...\n")}
  kcv = nr_folds
  if(ncol(X) < max(a)+max(ax,ay) | ncol(Y) < max(a)+max(ay,ay)) 
    warning("Some combinations of # components exceed data dimensions, these combinations are not considered\n")
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
                          stripped = stripped, p_thresh = p_thresh, 
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
                          stripped = stripped, p_thresh = p_thresh, 
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