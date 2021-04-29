library(OmicsPLS)
library(magrittr)

###############################
# Change settings here
N=50
p=1000 # please set to multiple of 100 to make sure the auto-generated group name is correct
q=100 # please set to multiple of 10 to make sure the auto-generated group name is correct
r=3
rx=1
ry=1
groupx = rep(paste0("g_", 1:100), p/100)
groupy = rep(paste0("g_", 1:10),q/10)

###############################
# Generate data
x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
params_true <- PO2PLS::generate_params(x, y, r, rx, ry, type='random')
dat <- PO2PLS::generate_data(N, params_true)
X <- dat$X %>% scale(scale = F)
Y <- dat$Y %>% scale(scale = F)
rm(x,y,dat)

colnames(X) <- paste0("x_", 1:p)
colnames(Y) <- paste0("y_", 1:q)
###############################
# Test

cv_nrcomp <- crossval_o2m_adjR2(X, Y, 1, 0:2, 0:2, nr_folds = 50)

print(cv_nrcomp)
crossval_o2m_adjR2(X, Y, 1, 1, 1, nr_folds = 50) %>% print

###############################
# use of O2PLS not affected
fit <- o2m(X, Y, r, rx, ry)
summary(fit)
scores(fit) %>% head
loadings(fit,loading_name = "Xjoint") %>% head
plot(fit, "Yjoint", i=1,j=3)

# SO2PLS
# Expected error, message to use O2PLS
fit <- o2m(X, Y, r, rx, ry, sparsity = T)

# SO2PLS with only one sparsity parameter specified
fit <- o2m(X, Y, r, rx, ry, sparsity = T, keepx = 30)
summary(fit)
fit$flags

# SO2PLS with non-integer sparsity parameter specified
fit <- o2m(X, Y, r, rx, ry, sparsity = T, keepx = 3.1, keepy = 10)

# Specify only groupx and keepx
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, keepx = 3)

# Specify only groupx and both sparsity level
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, keepx = 3, keepy = 3)

# Specify only groupx and both sparsity level
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, keepx = c(3,10,5), keepy = 3)

# Specify all, with non-integer sparisity
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, groupy=groupy, keepx = 3.1, keepy = 3)

# Expected error: keepx exceeds number of groups
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, groupy=groupy, keepx = 10000, keepy = 3)

#################################
# CV for sparsity
# SO2PLS
cv <- cv_sparsity(X, Y, r, rx, ry, lambda_kcv = 2, 
                  keepx_seq = seq(50,1000, by=150), keepy_seq = 3:5)
cv

# GO2PLS
cv <- cv_sparsity(X, Y, r, rx, ry, lambda_kcv = 5, groupx = groupx, groupy = groupy,
                  keepx_seq = seq(10,100, by=10), keepy_seq = 2:4)
cv
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, groupy=groupy, 
           keepx = cv$x_1sd, keepy = cv$y_1sd)

summary(fit)

#################################
# Generic functions
fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, groupy=groupy, keepx = 10, keepy = 3)

fit %>% str
summary(fit)
fit$flags$method
# loading at individual level
loadings(fit, "Yjoint") %>% head
# loading at group level
loadings(fit, "Yjoint_gr") %>% head
# scores
scores(fit) %>% head
# loading plot at individual level
plot(fit, "Xjoint", i=1,j=2, label = 'colnames')
# loading plot at group level
plot(fit, "Xjoint_gr", i=1,j=2, label = 'colnames')
