#' @export
spo2m <- function(X, Y, n, nx, ny, tol = 1e-10, max_iterations = 100, lambda_kcv = 2, 
                keepx = NULL, keepy = NULL) {
#  nn <- n+max(nx,ny)
  fitini <- o2m(X,Y,n,0,0,p_thresh=1,sparsity_it = T, keepx = keepx, keepy = keepy)
  fito2m <- o2m(X,Y,n,0,0,p_thresh=1)
  Ee <- X - fito2m$Tt %*% t(fito2m$W.)
  Ff <- Y - fito2m$U %*% t(fito2m$C.)
  
  Wosc <- svd(t(Ee)%*%fitini$Tt, nu = nx, nv = 0)$u
  Tosc <- X %*% Wosc
  Posc <- t(solve(t(Tosc) %*% Tosc) %*% t(Tosc) %*% X)
  
  Cosc <- svd(t(Ff)%*%fitini$U, nu = ny, nv = 0)$u
  Uosc <- Y %*% Cosc
  Qosc <- t(solve(t(Uosc) %*% Uosc) %*% t(Uosc) %*% Y)
  
  x <- X - Tosc %*% t(Posc)
  y <- Y - Uosc %*% t(Qosc)
  
  # update loadings
  # fit <- o2m(x,y,n,0,0,p_thresh=1,sparsity_it = T, keepx = keepx, keepy = keepy)
  
  Tt <- x %*% fitini$W.
  U <- y %*% fitini$C.
  
  B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  
  model <- list(Tt = Tt, W. = fitini$W., U = U, C. = fitini$C., T_Yosc = Tosc, P_Yosc. = Posc, W_Yosc = Wosc, 
                U_Xosc = Uosc, P_Xosc. = Qosc, C_Xosc = Cosc, B_U = B_U, B_T. = B_T)
  return(model)
}