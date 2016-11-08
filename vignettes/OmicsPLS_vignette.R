## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=8, fig.path='Figs/',
                      echo=TRUE, warning=TRUE, message=TRUE)

## ------------------------------------------------------------------------
set.seed(564785412L)
X = rnorm(100) %*% t(c(rep(1,5), rep(0,45))/sqrt(5)) + # Component 1 = joint
  rnorm(100) %*% t(c(rep(0,45), rep(1,5))/sqrt(5)) # Component 2 = specific
Y = X[,c(6:25, 1:5, 26:45)] # Reorder columns of X and leave out last 5
X = X + matrix(rnorm(100*50), nrow=100) # add noise
Y = Y + matrix(rnorm(100*45), nrow=100) # add noise

X = scale(X, scale=F)
Y = scale(Y, scale=F)

## ------------------------------------------------------------------------
try(
  gplots::heatmap.2(cor(X,Y), Rowv=F,Colv=F, col=gplots::bluered(100)),
  silent = TRUE)

## ------------------------------------------------------------------------
library(OmicsPLS)
set.seed(1221L)
crossval_o2m_adjR2(X, Y, 1:3, 0:3, 0:3, nr_folds = 2)
crossval_o2m(X, Y, 1:3, 0:3, 0:3, nr_folds = 10)

## ------------------------------------------------------------------------
fit0 = o2m(X, Y, 1, 1, 0)
fit0
summary(fit0)

## ------------------------------------------------------------------------
plot(fit0)
plot(fit0, "Yj")

## ------------------------------------------------------------------------
plot(fit0, "Xo")

## ----Function: filter genes----------------------------------------------
filter_rna <- function(rna=rna, prop = 0.75){
  #first, calculate the maximum of gene expression per each gene and take the quantiles, we are interested in top 25% of the gene expressions
  maxGE <- apply(rna, 2, max)
  propGEmax <- quantile(maxGE, prop)
  #next, take the IQR of each gene, to capture the variability and check the top 25% genes according to the size of IQR
  IQRGE <- apply(rna, 2, IQR, na.rm=TRUE) 
  propGEIQR <- quantile(IQRGE, prop)
  #selected genes/probes is the intersection of the two previous sets
  filter2 <- (intersect(which(maxGE> propGEmax), which(IQRGE> propGEIQR)))
  return(filter2)
}

## ----Load RNA data-------------------------------------------------------
set.seed(31*12*2016)
rna0 <- data.table::fread("C:/Users/selbouhaddani/MEGANZ/Downloads/O2PLS_BiB/test.tab",header=T, sep='\t')
rna1 <- t(as.data.frame.matrix(rna0[-1,-1,with=F]))
rna2 <- matrix(as.numeric(rna1), nrow = nrow(rna1))
dimnames(rna2) <- list(colnames(rna0)[-1],unlist(rna0[-1,1,with=F]))
rna2 <- rna2[order(row.names(rna2)), ] # Order rows according to the participant ID
rna3 <- rna2[,filter_rna(rna2)]
rm(rna0)
rm(rna1)

## ----Load metabolite data------------------------------------------------
metab0 <- read.table("C:/Users/selbouhaddani/MEGANZ/Downloads/O2PLS_BiB/metabonomic_data.txt",header=T, sep='\t')
metab1 <- t(metab0[,-1])
colnames(metab1) <- metab0$Metabolite

## ----Visualize missingness-----------------------------------------------
VIM::aggr(metab1, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
     labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

## ----Remove complete missings--------------------------------------------
NAs_in_metab1 <- which(apply(metab1, 1, function(e) sum(is.na(e))/length(e))==1)
metab2 <- metab1[-NAs_in_metab1,]
rna4 <- rna3[-NAs_in_metab1,]

## ----Impute metabolites, cache=TRUE--------------------------------------
metab2.imp <- missForest::missForest(metab2, verbose = T)
X2 <- scale(metab2.imp$ximp, scale=F)
X1 <- scale(rna4, scale = F)

## ----Heatmap of correlations---------------------------------------------
try(
  gplots::heatmap.2(cor(metab1,use = 'pair'), Rowv=F, Colv=F, 
                    breaks=seq(-1,1,length.out = 25), col=gplots::bluered),
  silent = TRUE)
try(
  gplots::heatmap.2(cor(X2,use = 'pair'), Rowv=F, Colv=F, 
                    breaks=seq(-1,1,length.out = 25), col=gplots::bluered),
  silent = TRUE)

## ----Eigenvalues plot----------------------------------------------------
svd(X1, 0, 0)$d[1:10]^2 / sum(X1^2)
svd(X2, 0, 0)$d[1:10]^2 / sum(X2^2)
svd(crossprod(X1,X2),0,0)$d[1:20]

## ----Boxplots------------------------------------------------------------
boxplot(X1[,filter_rna(X1, .95)])
boxplot(X2)

## ----Cross-Validation----------------------------------------------------
library(OmicsPLS)
set.seed(1+1+2016)
#CV1 <- crossval_o2m_adjR2(X1, X2, 1:3, c(0,1,5,10), c(0,1,5,10), nr_folds = 2, nr_cores = 4)
#CV2 <- crossval_o2m(X1, X2, 1:2, 0:2, 9:11, nr_folds = 10, nr_cores = 4)
#CV1
#CV2

## ----Fit O2PLS-----------------------------------------------------------
fit = o2m(X1, X2, 2, 1, 10)
fit

## ----Summary of the fit--------------------------------------------------
summary(fit)

## ----Loadings plot-------------------------------------------------------
require(magrittr)
require(ggplot2)
require(gridExtra)
# Color names
name_col = 1 + sapply( #First sapply loops over column names
  X = colnames(X2),
  FUN = function(arg){
    crossprod(
      c(1, 1, 3), # Weights to be used as categories
      sapply(c("VLDL", "LDL", "HDL"), # metabolite classes
             function(arg2){grepl(arg2, arg)} # compare class of metabolites
      )
    )
    }
  )

alpX2 <- loadings(fit, "Yjoint", 1:2) %>%  # Retreive loadings
  abs %>% # Absolute loading values for positive weights
  rowSums %>% # Sum over the components
  sqrt + (name_col>1) # Take square root

######### Plot loadings with OmicsPLS plot method ###
p_X2 <- plot(fit, loading_name="Yj", i=1, j=2, label="c", # Plot the loadings
             alpha=0) + # set points to be 100% transparant
##################### Add all layers ###
  theme_bw() + 
  coord_fixed(ratio = 1, xlim=c(-.2,.2),ylim=c(-.2,.2)) +
  geom_point( # Set color and size
    aes(col=factor(name_col, levels = 4:1), size = I(1+(name_col>1)), shape = 
          factor(name_col, levels = 4:1)),show.legend = T) +
  theme(legend.position="right") + 
  scale_color_discrete(name="Metabolite\nGroup", 
                       labels=c("LDL","VLDL","HDL","Other")) +
  guides(size=F) + scale_shape_discrete(name="Metabolite\nGroup", 
                                        labels=c("LDL","VLDL","HDL","Other")) +
labs(title = "Metabolite joint loadings", 
     x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'), 
        legend.title=element_text(face='bold')) 

alpX1 <- loadings(fit, "Xjoint", 1:2) %>% raise_to_power(2) %>% rowSums
alpX1[-(order(alpX1,decreasing=T)[1:10])] = 0
alpX1 <- sign(alpX1)
######### Plot loadings with OmicsPLS plot method ###
p_X1 <- plot(fit, loading_name="Xj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
             alpha = alpX1, 
             aes(label = stringr::str_sub(colnames(X1), start = 6)),
             hjust = rep(c(0, 1), length.out = ncol(X1))) + 
##################### Add all layers ###
  theme_bw() + 
  coord_fixed(.8, c(-.15,.15),c(-.15,.15)) +
  geom_point(alpha = 0.5+0.5*alpX1, col = 'grey') + 
  labs(title = "Transcript joint loadings", 
       x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'))

## Finally plot both plots in one figure.
grid.arrange(p_X2, p_X1, ncol=2)

## ----Heatmap decomposition-----------------------------------------------
library(gplots)
hm.2 <- function(obj){
  try(
    heatmap.2(obj,Rowv=F,symm=T,col=colorRampPalette(c("blue","white","red"))(25),
            dendrogram='none',trace='none',symkey=TRUE,breaks=seq(-1,1,length=26),
            key = T),
    silent = TRUE)
}
hm.2(cor(with(fit,U%*%t(C.))))
hm.2(cor(with(fit,U_Xosc%*%t(P_Xosc.))))
hm.2(cor(with(fit,Ff)))

## ----CPU timings, cache=TRUE---------------------------------------------
set.seed(2016^2)
fake_X <- scale(matrix(rnorm(1e2*1e4),1e2))
fake_Y <- scale(matrix(rnorm(1e2*1e2),1e2))
suppressMessages(
  scenario1 <- microbenchmark::microbenchmark(
    default=o2m(fake_X, fake_Y, 1, 1, 1),
    stripped=o2m(fake_X, fake_Y, 1, 1, 1, stripped=T),
    highD = o2m(fake_X, fake_Y, 1, 1, 1, stripped=T, p_thresh=1),
    times = 5, unit = 's',control=list(warmup=1))
)

fake_X <- scale(matrix(rnorm(1e2*2e3),1e2))
fake_Y <- scale(matrix(rnorm(1e2*2e3),1e2))
suppressMessages(
  scenario2 <- microbenchmark::microbenchmark(
    default=o2m(fake_X, fake_Y, 1, 1, 1),
    stripped=o2m(fake_X, fake_Y, 1, 1, 1, stripped=T),
    highD = o2m(fake_X, fake_Y, 1, 1, 1, stripped=T, p_thresh=1),
    times = 5, unit = 's',control=list(warmup=1))
)

fake_X <- scale(matrix(rnorm(1e2*5e4),1e2))
fake_Y <- scale(matrix(rnorm(1e2*5e4),1e2))
try(o2m(fake_X, fake_Y, 1, 1, 1, p_thresh=1e6))
try(o2m(fake_X, fake_Y, 1, 1, 1, stripped=T, p_thresh=1e6))
o2m(fake_X, fake_Y, 1, 1, 1, stripped=T)
rm(fake_X)
rm(fake_Y)


