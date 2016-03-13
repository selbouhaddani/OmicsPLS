crossval_o2m <- function(X, Y, a = 1:2, a2 = 1, b2 = 1, kcv, app_err = F) {
  stopifnot(ncol(X) > max(a)+max(a2) , ncol(Y) > max(a)+max(b2) , nrow(X) >= kcv)
  
  parms = data.frame(nx = a2)
  parms = merge(parms,data.frame(ny = b2))
  parms = merge(parms,data.frame(a = a))
  parms = apply(parms,1,as.list)
  outp=mclapply(mc.cores=detectCores(),parms,function(e){
    loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=app_err,func=o2m_stripped,kcv=kcv)[[1]]
  })
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
print.cvo2m <- function(x,include_matrix=FALSE) {
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