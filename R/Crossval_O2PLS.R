#' Cross-validate procedure for O2PLS
#'
#' @inheritParams o2m
#' @param a Vector of positive integers. Denotes the numbers of joint components to consider.
#' @param a2 Vector of non-negative integers. Denotes the numbers of X-specific components to consider.
#' @param b2 Vector of non-negative integers. Denotes the numbers of Y-specific components to consider.
#' @param kcv Positive integer. Number of folds to consider. Note: \code{kcv=N} gives leave-one-out CV.
#' @param nr_cores Positive integer. Number of cores to use for CV. You might want to use \code{\link{detectCores}()}.
#' 
#' @return List of class \code{"cvo2m"} with the original and sorted Prediction errors and the number of folds used.
#' 
#' @examples
#' local({
#' X = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' Y = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' crossval_o2m(X, Y, a = 1:4, a2 = 1:2, b2 = 1:2, 
#'              kcv = 5, nr_cores = 1)
#' })
#' 
#' @export
crossval_o2m <- function(X, Y, a = 1, a2 = 1, b2 = 1, kcv, nr_cores = 1,
                         stripped = TRUE, p_thresh = 3000, 
                         q_thresh = p_thresh, tol = 1e-10, max_iterations = 100) {
  stopifnot(ncol(X) > max(a)+max(a2) , ncol(Y) > max(a)+max(b2) , nrow(X) >= kcv)
  stopifnot(nr_cores == abs(round(nr_cores)))
  parms = data.frame(nx = a2)
  parms = merge(parms,data.frame(ny = b2))
  parms = merge(parms,data.frame(a = a))
  parms = apply(parms,1,as.list)
  cl_crossval_o2m <- NULL
  
  on.exit({if(!is.null(cl_crossval_o2m)) stopCluster(cl_crossval_o2m)})
  
  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_o2m <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_o2m, library(O2PLS))
    clusterExport(cl_crossval_o2m, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_o2m,parms,function(e){
      loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                  stripped = stripped, p_thresh = p_thresh, 
                  q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]]
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                  stripped = stripped, p_thresh = p_thresh, 
                  q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]]
    })
  }
  dnams = list(paste("a2=",a2,sep=""),paste("b2=",b2,sep=""),paste("a=",a,sep=""))
  outp1 = array(unlist(outp),dim=c(length(a2),length(b2),length(a)),
                dimnames=dnams)
  outp = aperm(outp1,order(dim(outp1),decreasing=TRUE))
  dnams = dimnames(outp)
  
  if(dim(outp)[3]==1){
    dim(outp) = dim(outp)[-3]
    dimnames(outp) = dnams[1:2]
  }
  
  outpt = list(Original=outp1,Sorted=outp,kcv=kcv)
  class(outpt) <- "cvo2m"
  return(outpt)
}

#' Adjusted Cross-validate procedure for O2PLS
#'
#' Combines CV with R2 optimization
#' 
#' @inheritParams crossval_o2m
#' 
#' @return List of class \code{"cvo2m"} with the original and sorted Prediction errors and the number of folds used.
#' @examples 
#' local({
#' X = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' Y = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' crossval_o2m_adjR2(X, Y, a = 1:4, a2 = 1:2, b2 = 1:2, 
#'              kcv = 5, nr_cores = 1)
#' })
#' @export
crossval_o2m_adjR2 <- function(X, Y, a = 1, a2 = 1, b2 = 1, kcv, nr_cores = 1,
                               stripped = TRUE, p_thresh = 3000, 
                               q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  stopifnot(ncol(X) > max(a)+max(a2) , ncol(Y) > max(a)+max(b2) , nrow(X) >= kcv)
  stopifnot(nr_cores == abs(round(nr_cores)))
  cl_crossval_o2m <- NULL
  on.exit({if(!is.null(cl_crossval_o2m)) stopCluster(cl_crossval_o2m)})
  
  parms = data.frame(a = a)
  parms = apply(parms,1,as.list)
  
  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_o2m <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_o2m, library(O2PLS))
    clusterExport(cl_crossval_o2m, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_o2m,parms,function(e){
      parms = data.frame(nx = a2)
      parms = merge(parms,data.frame(ny = b2))
      parms = apply(parms,1,as.list)
      R2grid = matrix(colMeans(suppressMessages(adjR2(Y, X, e$a, a2, b2,
                                                      stripped = stripped, p_thresh = p_thresh, 
                                                      q_thresh = q_thresh, tol = tol, max_iterations = max_iterations))), 
                      nrow = length(b2), byrow=TRUE)
      nxny = which(R2grid == max(R2grid), arr.ind = TRUE)[1,]
      a_mse = loocv_combi(X,Y,e$a,a2[nxny[2]],b2[nxny[1]],app_err=app_err,func=o2m,kcv=kcv,
                          stripped = stripped, p_thresh = p_thresh, 
                          q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]]
      c(a_mse, e$a, a2[nxny[2]],b2[nxny[1]])
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      parms = data.frame(nx = a2)
      parms = merge(parms,data.frame(ny = b2))
      parms = apply(parms,1,as.list)
      R2grid = matrix(colMeans(suppressMessages(adjR2(Y, X, e$a, a2, b2,
                                     stripped = stripped, p_thresh = p_thresh, 
                                     q_thresh = q_thresh, tol = tol, max_iterations = max_iterations))), 
                      nrow = length(b2), byrow=TRUE)
      nxny = which(R2grid == max(R2grid), arr.ind = TRUE)[1,]
      a_mse = loocv_combi(X,Y,e$a,a2[nxny[2]],b2[nxny[1]],app_err=app_err,func=o2m,kcv=kcv,
                          stripped = stripped, p_thresh = p_thresh, 
                          q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]]
      c(a_mse, e$a, a2[nxny[2]],b2[nxny[1]])
    })
  }
  outp2 = matrix(unlist(outp), nrow = length(a), byrow = T)
  outp2 <- as.data.frame(outp2)
  names(outp2) <- c("MSE", "n", "nx", "ny")
  message("minimum is at n = ", outp2[,2][which.min(outp2[,1])])
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
  wmCV = which(min(x$Or)==x$Or,TRUE,FALSE)
  dnams = dimnames(x$Or)
  dnams1 = dnams[[1]][wmCV[1]]
  dnams2 = dnams[[2]][wmCV[2]]
  dnams3 = dnams[[3]][wmCV[3]]
  
  cat("*******************\n")
  cat("Minimal ",x$kcv,"-CV error is at ",dnams1," ",dnams2," ",dnams3," ","\n",sep="")
  cat("*******\n")
  cat("Minimum is",min(x$Sor),"\n")
  if(include_matrix){
    cat("*******\n")
    cat("Simplified CV matrix is \n")
    print(x$Sorted)
  }
  cat("*******************\n")
}