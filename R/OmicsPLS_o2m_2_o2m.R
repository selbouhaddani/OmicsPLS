o2m_to_o2m <- function(fit){
  LVs <- with(fit, list(T=Tt, U=U, To=T_Yosc, Uo=U_Xosc) %>% as.data.frame %>% as_tibble)
  pars <- with(fit, list(W=W., C=C., Wo=P_Yosc., Co=P_Xosc.))
  R2s <- with(fit, c(
    R2X = R2X,
    R2Y = R2Y,
    R2Xj = R2Xcorr,
    R2Yj = R2Ycorr,
    R2Xhat = R2Xhat,
    R2Yhat = R2Yhat
  ))
  outp <- list(parameters = pars,
               latent_vars = LVs,
               meta_data = list(
                 loglikelihood = NA,
                 explained_vars = R2s,
                 comps = with(pars, c(r=ncol(W), rx=ncol(Wo),ry=ncol(Co))),
                 time = fit$flags$time,
                 call = fit$flags$call,
                 convergence = NA,
                 ssqs = c(X=fit$flags$ssqX, Y=fit$flags$ssqY),
                 stripped = fit$flags$stripped,
                 highd = fit$flags$highd
               ))
  class(outp) <- "o2m"
  return(outp)
}

print.o2m <- function(x, ...){
  cat('\n')
  cat("Call:", x$meta_data$call %>% deparse)
  cat('\n')
  
  cat('\n')
  cat(paste0("This is a", ifelse(sum(x$meta_data$comps[2:3])==0, " ", "n O2"), "PLS fit"))
  cat('\n')
  
  #  cat('\n')
  cat(paste0("With r=",x$meta$comps[1],", rx=",x$meta$comps[2]," and ry=", x$meta$comps[3]," components"))
  cat('\n')
  
  #  cat('\n')
  cat(paste0("Elapsed time: ", signif(x$meta$time,3), " sec"))
  cat('\n')
  
  if(x$meta$stripped | x$meta$highd)
  cat(paste0("Note: ", 
             ifelse(x$meta$stripped, ifelse(x$meta$highd, "Stripped and High Dim", "Stripped"), 
                    "High Dim"), " mode active"))
  cat('\n')
  
}

#
summary.o2m <- function(object, ...){
  outp <- "\n Not implemented yet \n"
  class(outp) <- "summary.po2m"
}

#
print.summary.o2m <- function(x, ...){
  cat(x)
}

