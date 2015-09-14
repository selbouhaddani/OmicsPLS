#' O2PLS: Two-Way Orthogonal Partial Least Squares
#' 
#' This is based on work of (Trygg, Wold, 2003). 
#' Includes the O2PLS fit, some misc functions and some cross-validation tools.
#' @author 
#' Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}), 
#' Jeanine Houwing-Duistermaat (\email{J.J.Houwing@@lumc.nl}), 
#' Geurt Jongbloed (\email{G.Jongbloed@@tudelft.nl}), 
#' Hae-Won Uh (\email{H.Uh@@lumc.nl}).
#' 
#' Maintainer: Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}).
#' 
#' @section Functions:
#' The O2PLS fit is done with \code{\link{o2m}}.
#' Cross-validation is done with \code{\link{loocv}} or \code{\link{adjR2}}, the last has built in parallelization (when you use Windows!) which relies on the \code{parallel} package.
#' 
#' List of all functions:\itemize{
#' \item{\code{\link{adjR2}}}
#' \item{\code{\link{loocv}}}
#' \item{\code{\link{loocv_combi}}}
#' \item{\code{\link{mse}}}
#' \item{\code{\link{orth}}}
#' \item{\code{\link{o2m}}}
#' \item{\code{\link{o2m_stripped}}}
#' \item{\code{\link{rmsep}}}
#' \item{\code{\link{rmsep_combi}}}
#' \item{\code{\link{ssq}}}
#' \item{\code{\link{summary.o2m}}}
#' \item{\code{\link{vnorm}}}
#' }
#' 
#' @docType package
#' @name O2PLS
#' @keywords O2PLS
#' @import parallel
NULL

#' Orthogonalize a matrix
#' 
#' @param X Numeric vector or matrix.
#' @param X_true (optional) A "true" matrix/vector. Used to correct the sign of the orthonormalized X if QR is used. Only the first column is corrected.
#' @param type A character or numeric. type="QR" and type=1 are equivalent.
#' @return An orthogonalized representation of \eqn{X}
#' @details Choosing type="QR" (or type=1) uses a QR decomposition of X to produce orthonormal columns. For type=="SVD" (or type not 1) it uses an SVD decomposition.
#' @examples
#' orth(c(3,4))
#' round(crossprod(orth(matrix(rnorm(500),100,5))),4)
#' orth(matrix(1:9,3,3),type="QR")[,1] - orth(1:3); orth(matrix(1:9,3,3),type="SVD")[,1] - orth(1:3);
#' @export
orth<-function(X,X_true=NULL,type=c("QR","SVD"))
{
  #input is some matrix X
  #output is a matrix with orthonormal columns (via SVD(X))
  if(is.character(type)) {type = match.arg(type)} else {type=type[1]}
  if(type!="QR" && type!=1){
    e=svd(X); 
    return(tcrossprod(e$u,e$v))
  }
  e = qr.Q(qr(X))
  sign_e = ifelse(!is.null(X_true),c(sign(crossprod(e[,1],as.matrix(X_true)[,1]))),c(sign(crossprod(e[,1],as.matrix(X)[,1]))))
  return(sign_e*e)
}

#' Calculate Sum of Squares
#' 
#' @param X Numeric vector or matrix.
#' @return The sum of squared elements of \eqn{X}
#' @details This is the Frobenius norm of \eqn{X}.
#' @examples
#' ssq(1:5)
#' ssq(rnorm(1e5))/1e5
#' @export
ssq<-function(X)
{
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
mse <- function(x,y=0,na.rm=TRUE)
  #Input: 2 matrices x,y
  #Output: mean squared distance between 2 matrices (number)
{
  if (length(x) != length(y) && length(y)!=1) {
    warning("unequal length:result may not be sensible")
  }
  mean((x - y)^2,na.rm=na.rm)
}

#' Perform O2-PLS with two-way orthogonal corrections
#' 
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALES THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#' 
#' @param X Numeric matrix. Other types will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param Y Numeric matrix. Other types will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param n Integer. Number of joint PLS components. Must be positive!
#' @param nx Integer. Number of orthogonal components in \eqn{X}. Negative values are interpreted as 0
#' @param ny Integer. Number of orthogonal components in \eqn{Y}. Negative values are interpreted as 0
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
#'    @details If both \code{nx} and \code{ny} are zero, \code{o2m} is equivalent to PLS2 with orthonormal loadings.
#'    This is a `slower' implementation of O2PLS, and is using \code{\link{svd}}. For cross-validation purposes, consider using \code{\link{o2m_stripped}}.
#'    
#'    @examples
#'      test.data=matrix(rnorm(100))
#'      hist(replicate(1000,
#'                    o2m(test.data,matrix(rnorm(100)),1,0,0)$B_T.
#'                  ),main="No joint variation",xlab="B_T",xlim=c(0,1.5));
#'      hist(replicate(1000,
#'                    o2m(test.data,test.data+rnorm(100),1,0,0)$B_T.
#'                  ),main="B_T=1; 25% joint variation",xlab="B_T",xlim=c(0,1.5));
#'      hist(replicate(1000,
#'                    o2m(test.data,test.data+rnorm(100,0,0.1),1,0,0)$B_T.
#'                  ),main="B_T=1; 90% joint variation",xlab="B_T",xlim=c(0,1.5));
#'    @seealso \code{\link{ssq}}, \code{\link{summary.o2m}}, \code{\link{o2m_stripped}}
#' @export
o2m<-function(X,Y,n,nx,ny)
{
  X = as.matrix(X)
  Y = as.matrix(Y)
  if(n<=0){stop("#joint components must be >0")}
  X_true = X
  Y_true = Y
  
  T_Yosc = U_Xosc = matrix(0,length(X[,1]),1)
  P_Yosc = W_Yosc = matrix(0,length(X[1,]),1)
  P_Xosc = C_Xosc = matrix(0,length(Y[1,]),1)
  
  N=dim(X)[1]
  p=dim(X)[2]; q=dim(Y)[2]
  
  if(nx*ny>0){
    # larger principal subspace
    n2=n+max(nx,ny)
    
    cdw = svd(t(Y)%*%X,nu=n2,nv=n2); 
    C=cdw$u;W=cdw$v
    
    Tt = X%*%W;                    
    
    if(nx > 0){
      # Orthogonal components in Y
      E_XY = X - Tt%*%t(W);
      
      udv = svd(t(E_XY)%*%Tt,nu=nx,nv=0);
      W_Yosc = udv$u ; s = udv$d 
      T_Yosc = X%*%W_Yosc;
      P_Yosc = t(solve(t(T_Yosc)%*%T_Yosc)%*%t(T_Yosc)%*%X);
      X = X - T_Yosc%*%t(P_Yosc);
      
      # Update T again
      Tt = X%*%W;
    }
    
    U = Y%*%C;                   # 3.2.1. 4
    
    if(ny > 0){
      # Orthogonal components in Y
      F_XY = Y - U%*%t(C);
      
      udv = svd(t(F_XY)%*%U,nu=ny,nv=0);
      C_Xosc = udv$u ; s = udv$d 
      U_Xosc = Y%*%C_Xosc;
      P_Xosc = t(solve(t(U_Xosc)%*%U_Xosc)%*%t(U_Xosc)%*%Y);
      Y = Y - U_Xosc%*%t(P_Xosc);
      
      # Update U again
      U = Y%*%C;
    }
  }
  # Re-estimate joint part in n-dimensional subspace
  cdw = svd(t(Y)%*%X,nu=n,nv=n);    # 3.2.1. 1
  C=cdw$u;W=cdw$v
  Tt = X%*%W;                    # 3.2.1. 2
  U = Y%*%C;                    # 3.2.1. 4
  
  # Inner relation parameters
  B_U = solve(t(U)%*%U)%*%t(U)%*%Tt;
  B_T = solve(t(Tt)%*%Tt)%*%t(Tt)%*%U;
  
  #Residuals and R2's
  E = X_true - Tt%*%t(W) - T_Yosc%*%t(P_Yosc);
  Ff = Y_true - U%*%t(C) - U_Xosc%*%t(P_Xosc);
  H_TU = Tt - U%*%B_U;
  H_UT = U - Tt%*%B_T;
  Y_hat = Tt%*%B_T%*%t(C);
  X_hat = U%*%B_U%*%t(W);
  
  R2Xcorr = (ssq(Tt%*%t(W))/ssq(X_true))
  R2Ycorr = (ssq(U%*%t(C))/ssq(Y_true))
  R2X_YO = (ssq(T_Yosc%*%t(P_Yosc))/ssq(X_true))
  R2Y_XO = (ssq(U_Xosc%*%t(P_Xosc))/ssq(Y_true))
  R2Xhat = 1 - (ssq(U%*%B_U%*%t(W) - X_true)/ssq(X_true))
  R2Yhat = 1 - (ssq(Tt%*%B_T%*%t(C) - Y_true)/ssq(Y_true))
  R2X = R2Xcorr + R2X_YO
  R2Y = R2Ycorr + R2Y_XO
  
  model=list(
    Tt=Tt,W.=W,U=U,C.=C,E=E,Ff=Ff,T_Yosc=T_Yosc,P_Yosc.=P_Yosc,W_Yosc=W_Yosc,U_Xosc=U_Xosc,P_Xosc.=P_Xosc,
    C_Xosc=C_Xosc,B_U=B_U,B_T.=B_T,H_TU=H_TU,H_UT=H_UT,X_hat=X_hat,Y_hat=Y_hat,R2X=R2X,R2Y=R2Y,
    R2Xcorr=R2Xcorr,R2Ycorr=R2Ycorr,R2X_YO=R2X_YO,R2Y_XO=R2Y_XO,R2Xhat=R2Xhat,R2Yhat=R2Yhat
  )
  class(model)='o2m'
  return(model)
}

#' Summary of an O2PLS fit
#' 
#' Until now only variational summary given by the R2's is outputted
#' 
#' @param fit List. Contains the R2's as produced by \code{\link{o2m}}.
#' @return Matrix with R2 values given in percentage in two decimals.
#' @examples
#' summary.o2m(o2m(matrix(-2:2),matrix(-2:2*4),1,0,0))
#' @export
summary.o2m<-function(fit)
{
  a=length(fit$W.[1,])
  Mname=list(c(""),c("Comp","R2X","R2Y","R2Xcorr",'R2Ycorr','R2Xhat',"R2Yhat","XRatio","YRatio"))
  M=matrix(c(a/100,fit$R2X,fit$R2Y,fit$R2Xcorr,fit$R2Ycorr,fit$R2Xhat,fit$R2Yhat,fit$R2Xhat/fit$R2Xcorr,fit$R2Yhat/fit$R2Ycorr),nrow=1,dimnames=Mname)
  return(round(100*M,2))
}

#' Root MSE of Prediction
#' 
#' Calculates the Root MSE of prediction on test data. Only tested to work inside \code{\link{loocv}}.
#' 
#' @param Xtst Numeric vector or matrix.
#' @param Ytst Numeric vector or matrix.
#' @param fit \code{\link{o2m}} fit (on data without \code{Xtst} and \code{Ytst}).
#' @details This function is the building block for \code{\link{loocv}}, as it produced the prediction error for test (left out) data.
#' @return Mean squares difference between predicted Y and true Y
#' @export
rmsep <- function(Xtst,Ytst,fit)
{
  n1=n2=T
  if(!is.matrix(Xtst)){Xtst=t(Xtst);n1=F} # If Xtst is a row-vector
  if(!is.matrix(Ytst)){Ytst=t(Ytst);n2=F} # If Xtst is a row-vector
  
  Yhat = Xtst%*%fit$W.%*%fit$B_T%*%t(fit$C.)
  Xhat = Ytst%*%fit$C.%*%fit$B_U%*%t(fit$W.)
  
  return(mean(c(sqrt(mse(Yhat,Ytst)))))#,sqrt(mse(Xhat,Xtst)))))
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
#' names_for_a=sapply(1:3,function(i){paste("a",i,sep="=")})
#' names_for_a2=sapply(0:1,function(i){paste("a2",i,sep="=")})
#' names_for_b2=sapply(0:1,function(i){paste("b2",i,sep="=")})
#' array(unlist(result$CVerr),dim=c(3,2,2),dimnames=list(names_for_a,names_for_a2,names_for_b2))
#' @export
loocv <- function(X,Y,a=1:2,a2=1,b2=1,fitted_model=NULL,func=o2m_stripped,app_err=F,kcv)
  #Input: Data X,Y; model to fit; loop through nr of components; 
  # calculate apparent error?; nr of folds (loo:kcv=N)
  #Output: several MSE's per chosen nr of component
{
  if(!is.null(fitted_model)){app_err=F;warning("apparent error calculated with provided fit")}
  p=dim(X)[2]
  # determine type of model
  type=3#ifelse(deparse(substitute(func))=="o2m",3,ifelse(deparse(substitute(func))=="oplsm",2,1))
  
  N = length(X[,1]);if(N!=length(Y[,1])){stop('N not the same')}
  mean_err=mean_fit=NA*1:max(length(a),length(a2),length(b2))
  k=0
  
  #blocks contains the begin and endpoints of test indices to use
  blocks = c(seq(0,N,by=floor(N/kcv)),N)
  
  #loop through chosen parameters
  for(j in a){for(j2 in a2){for(j3 in b2){
    k=k+1
    err=NA*1:kcv
    folds=sample(N)
    #loop through number of folds
    for(i in 1:kcv){
      ii = (blocks[i]+1):(blocks[i+1])
      if(type==3){pars=list(X=X[-folds[ii],],Y=Y[-folds[ii],],n=j,nx=j2,ny=j3)}
      if(type==2){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j,n_orth=j2)}
      if(type==1){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j)}
      fit=try(do.call(func,pars),silent=T)
      err[i] = ifelse(class(fit)=="try-error",NA,rmsep(X[folds[ii],],Y[folds[ii],],fit))
    }
    mean_err[k]=mean(err)
    #calculate apparent error
    if(app_err && is.null(fitted_model)){
      if(type==3){pars2=list(X=X,Y=Y,n=j,nx=j2,ny=j3)}
#      if(class(fit)=='oplsm'){pars2=list(X=X,Y=Y,ncomp=j,n_orth=j2)}
#      if(class(fit)=='plsm'){pars2=list(X=X,Y=Y,ncomp=j)}
      fit2=try(do.call(func,pars2),F)
      mean_fit[k]=ifelse(class(fit)=="try-error",NA,rmsep(X,Y,fit2))
#      print('1e loop')
    }
    if(!is.null(fitted_model))
    {
      mean_fit[k]=rmsep(X,Y,fitted_model)
#      print('2e loop')
    }
  }}}
  return(list(CVerr=mean_err,Fiterr=mean_fit))
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
#' @param cl Object of class "\code{cluster}". If parall is \code{TRUE} and \code{cl} is not \code{NULL}, calculations are parallelized over workers in cl.
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
#' names_for_a=sapply(1:3,function(i){paste("a",i,sep="=")})
#' names_for_a2=sapply(0:1,function(i){paste("a2",i,sep="=")})
#' names_for_b2=sapply(0:1,function(i){paste("b2",i,sep="=")})
#' array(unlist(result[1,]),dim=c(3,2,2),dimnames=list(names_for_a,names_for_a2,names_for_b2))
#' @export
adjR2 <- function(X,Y,a=1:2,a2=1,b2=1,func=o2m_stripped,parall=F,cl=NULL)
  #Input: Data X,Y; model to fit; loop through nr of components; 
  # calculate apparent error?; nr of folds (loo:kcv=N)
  #Output: several MSE's per chosen nr of component
{
  X;Y;
  cl_was_null=FALSE
  if(!parall){S_apply=function(cl=NULL,x,fun){sapply(x,fun)}}
  if(parall&is.null(cl)){
    cl_was_null=TRUE
    S_apply=parSapply
    cl <- makeCluster(rep('localhost', detectCores()),type='SOCK')
    clusterExport(cl=cl, varlist=c("ssq",'o2m_stripped',"adjR2"))
  }
  if(parall&!is.null(cl)){S_apply=parSapply}
  
  pars1 = merge(merge(data.frame(a = a),data.frame(a2=a2)),data.frame(b2=b2))
  pars2 = apply(pars1, 1, as.list)
  N = dim(X)[1]
  #cl <- makeCluster(rep( 'localhost', detectCores()),type='SOCK')
  #clusterExport(cl=cl, varlist=c("X","Y","N","pars2","ssq",'o2m'))
  outp=S_apply(cl,pars2,function(p){
    fit=func(X,Y,p$a,p$a2,p$b2);
    RX = 1-ssq(fit$H_UT)/ssq(fit$U)
    RY = 1-ssq(fit$H_TU)/ssq(fit$Tt)
    adjRX = RX#1 - (1 - RX)*(N - 1)/(N - p$a - 1)
    adjRY = RY#1 - (1 - RY)*(N - 1)/(N - p$a - 1)
    return(c(adjR2X = adjRX, adjR2Y = adjRY))
  })
  if(parall&cl_was_null==TRUE){stopCluster(cl)}
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
  #Input: matrix x
  #Output: 2-norm of X per column ( ssq(X) = sum(vnorm(X)) )
{
  x = as.matrix(x)
  return(sqrt(apply(x^2,2,sum)))
}
#' Perform O2-PLS with two-way orthogonal corrections
#' 
#' NOTE THAT THIS FUNCTION DOES NOT CENTER NOR SCALES THE MATRICES! Any normalization you will have to do yourself. It is best practice to at least center the variables though.
#' A stripped version of O2PLS
#' 
#' @param X Numeric matrix. Other types will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param Y Numeric matrix. Other types will be coerced to matrix with \code{as.matrix} (if this is possible)
#' @param n Integer. Number of joint PLS components. Must be positive!
#' @param nx Integer. Number of orthogonal components in \eqn{X}. Negative values are interpreted as 0
#' @param ny Integer. Number of orthogonal components in \eqn{Y}. Negative values are interpreted as 0
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
#'    @details If both \code{nx} and \code{ny} are zero, \code{o2m} is equivalent to PLS2 with orthonormal loadings.
#'    This is a stripped implementation of O2PLS, using \code{\link{svd}}. For data analysis purposes, consider using \code{\link{o2m}}.
#'    
#'    @seealso \code{\link{ssq}}, \code{\link{o2m}}, \code{\link{loocv}}, \code{\link{adjR2}}
#' @export
o2m_stripped<-function(X,Y,n,nx,ny)
{
  X = as.matrix(X)
  Y = as.matrix(Y)
  if(n<=0){stop("#joint components must be >0")}
  
  X_true = X
  Y_true = Y
  
  N=dim(X)[1]
  p=dim(X)[2]; q=dim(Y)[2]
  
  T_Yosc = U_Xosc = matrix(NA,N,1)
  P_Yosc = W_Yosc = matrix(NA,p,1)
  P_Xosc = C_Xosc = matrix(NA,q,1)
  
  if(nx+ny>0){
    n2=n+max(nx,ny)
    
    cdw = svd(t(Y)%*%X,nu=n2,nv=n2); # 3.2.1. 1
    C=cdw$u;W=cdw$v
    rm(cdw)
    Tt = X%*%W;                    # 3.2.1. 2
    
    if(nx > 0){
      # 3.2.1. 3
      E_XY = X 
      E_XY = E_XY - Tt%*%t(W);
      
      udv = svd(t(E_XY)%*%Tt,nu=nx,nv=0);
      rm(E_XY)
      W_Yosc = udv$u 
      T_Yosc = X%*%W_Yosc;
      P_Yosc = t(solve(t(T_Yosc)%*%T_Yosc)%*%t(T_Yosc)%*%X);
      X = X - T_Yosc%*%t(P_Yosc);
      
      # Update T again (since X has changed)
      Tt = X%*%W;
    }
    
    U = Y%*%C;                   # 3.2.1. 4
    
    if(ny > 0){
      # 3.2.1. 5
      F_XY = Y 
      F_XY = F_XY - U%*%t(C);
      
      udv = svd(t(F_XY)%*%U,nu=ny,nv=0);
      rm(F_XY)
      C_Xosc = udv$u ; s = udv$d 
      U_Xosc = Y%*%C_Xosc;
      P_Xosc = t(solve(t(U_Xosc)%*%U_Xosc)%*%t(U_Xosc)%*%Y);
      Y = Y - U_Xosc%*%t(P_Xosc);
      
      # Update U again (since Y has changed) 
      U = Y%*%C;
    }
  }
  
  # repeat steps 1, 2, and 4 before step 6
  cdw = svd(t(Y)%*%X,nu=n,nv=n);    # 3.2.1. 1
  C=cdw$u;W=cdw$v
  Tt = X%*%W;                    # 3.2.1. 2
  U = Y%*%C;                    # 3.2.1. 4
  
  # 3.2.1. 6
  B_U = solve(t(U)%*%U)%*%t(U)%*%Tt;
  B_T = solve(t(Tt)%*%Tt)%*%t(Tt)%*%U;
  H_TU = Tt - U%*%B_U;
  H_UT = U - Tt%*%B_T;
  
  model=list(Tt=Tt,U=U,W.=W,C.=C,P_Yosc.=P_Yosc,P_Xosc.=P_Xosc,B_T.=B_T,B_U=B_U,H_TU=H_TU,H_UT=H_UT)
  class(model)='o2m'
  return(model)
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
rmsep_combi <- function(Xtst,Ytst,fit)
  #Input: test matrices (now: row i of data X,Y); used model
  #Output: MSE of Yhat from model, wrt Ytest
{
  n1=n2=T
  if(!is.matrix(Xtst)){Xtst=t(Xtst);n1=F} # If Xtst is a row-vector
  if(!is.matrix(Ytst)){Ytst=t(Ytst);n2=F} # If Xtst is a row-vector
  
  if(class(fit)=='plsm')
    # PLS fit is quickly done
  {
    Yhat=cbind(1,Xtst)%*%fit$ori
  }
  
  if(class(fit)=='oplsm')
    # OPLS corrects Xtst
  {
    if(n1){Xtst=oscr(Xtst,Ytst,n_orth=n_orth)$Xcorr}
    Yhat=cbind(1,Xtst)%*%fit$ori
  }
  
  if(class(fit)=='o2m')
    # O2PLS we should correct both Xtst and Ytst?
  {
    if(n1)
    {
      #To=Xtst%*%fit$W_Yosc
      #Po=t(solve(t(To)%*%To)%*%t(To)%*%Xtst)
      Xtst=Xtst# - tcrossprod(To,Po)
    }
    Yhat = Xtst%*%fit$W.%*%fit$B_T%*%t(fit$C.)
    Xhat = Ytst%*%fit$C.%*%fit$B_U%*%t(fit$W.)
  }
  return( sqrt(mse(Yhat,Ytst)) + sqrt(mse(Xhat,Xtst)) )
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
loocv_combi <- function(X,Y,a=1:2,a2=1,b2=1,fitted_model=NULL,func=o2m_stripped,app_err=F,kcv)
  #Input: Data X,Y; model to fit; loop through nr of components; 
  # calculate apparent error?; nr of folds (loo:kcv=N)
  #Output: several MSE's per chosen nr of component
{
  if(!is.null(fitted_model)){app_err=F;warning("apparent error calculated with provided fit")}
  p=dim(X)[2]
  # determine type of model
  type=3#ifelse(deparse(substitute(func))=="o2m",3,ifelse(deparse(substitute(func))=="oplsm",2,1))
  
  N = length(X[,1]);if(N!=length(Y[,1])){stop('N not the same')}
  mean_err=mean_fit=NA*1:max(length(a),length(a2),length(b2))
  k=0
  
  #blocks contains the begin and endpoints of test indices to use
  blocks = c(seq(0,N,by=floor(N/kcv)),N)
  
  #loop through chosen parameters
  for(j in a){for(j2 in a2){for(j3 in b2){
    k=k+1
    err=NA*1:kcv
    folds=sample(N)
    #loop through number of folds
    for(i in 1:kcv){
      ii = (blocks[i]+1):(blocks[i+1])
      if(type==3){pars=list(X=X[-folds[ii],],Y=Y[-folds[ii],],n=j,nx=j2,ny=j3)}
      if(type==2){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j,n_orth=j2)}
      if(type==1){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j)}
      fit=try(do.call(func,pars),silent=T)
      err[i] = ifelse(class(fit)=="try-error",NA,rmsep_combi(X[folds[ii],],Y[folds[ii],],fit))
    }
    mean_err[k]=mean(err)
    #calculate apparent error
    if(app_err && is.null(fitted_model)){
      if(class(fit)=='o2m'){pars2=list(X=X,Y=Y,n=j,nx=j2,ny=j3)}
      if(class(fit)=='oplsm'){pars2=list(X=X,Y=Y,ncomp=j,n_orth=j2)}
      if(class(fit)=='plsm'){pars2=list(X=X,Y=Y,ncomp=j)}
      fit2=try(do.call(func,pars2),T)
      mean_fit[k]=ifelse(class(fit)=="try-error",NA,rmsep_combi(X,Y,fit2))
      print('1e loop')
    }
    if(!is.null(fitted_model))
    {
      mean_fit[k]=rmsep_combi(X,Y,fitted_model)
      print('2e loop')
    }
  }}}
  return(list(CVerr=mean_err,Fiterr=mean_fit))
}