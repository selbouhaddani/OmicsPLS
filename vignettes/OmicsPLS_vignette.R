## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=7.5, fig.height=4.5, fig.path='Figs/',
                      echo=TRUE, warning=TRUE, message=TRUE)

## -----------------------------------------------------------------------------
set.seed(564785412L)
X = rnorm(100) %*% t(c(rep(1,5), rep(0,45))/sqrt(5)) + # Component 1 = joint
  rnorm(100) %*% t(c(rep(0,45), rep(1,5))/sqrt(5)) # Component 2 = specific
Y = X[,c(6:25, 1:5, 26:45)] # Reorder columns of X and leave out last 5
X = X + matrix(rnorm(100*50), nrow=100) # add noise
Y = Y + matrix(rnorm(100*45), nrow=100) # add noise

X = scale(X, scale=F)
Y = scale(Y, scale=F)

## -----------------------------------------------------------------------------
try(
  gplots::heatmap.2(cor(X,Y), Rowv=F,Colv=F, col=gplots::bluered(100),
                    symm = TRUE, trace="none", dendrogram="none"),
  silent = TRUE)

## -----------------------------------------------------------------------------
library(OmicsPLS)
set.seed(1221L)
crossval_o2m_adjR2(X, Y, 1:3, 0:3, 0:3, nr_folds = 2)
crossval_o2m(X, Y, 1:3, 0:3, 0:3, nr_folds = 10)

## -----------------------------------------------------------------------------
fit0 = o2m(X, Y, 1, 1, 0)
fit0
summary(fit0)

## -----------------------------------------------------------------------------
plot(fit0)
plot(fit0, "Yj")

## -----------------------------------------------------------------------------
plot(fit0, "Xo")

