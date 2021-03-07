library(OmicsPLS)
source("/home/z/Mygit/Supervised_PO2PLS/SuPO2PLS.R")

N=100
p=10000
q=100
r=3
rx=1
ry=1

groupx = rep(1:100, 100)
groupy = rep(1:10,10)

x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
params_true <- generate_params(x, y, r, rx, ry, type='random')
dat <- generate_data(N, params_true)
X <- dat[1:N,1:p] %>% scale(scale = F)
Y <- dat[1:N,(p+1):(p+q)] %>% scale(scale = F)

fit <- o2m(X, Y, r, rx, ry)
summary(fit)
plot(fit, "Yjoint", i=1,j=3)


fit <- o2m(X, Y, r, rx, ry, sparsity = T, keepx = 3)

fit <- o2m(X, Y, r, rx, ry, sparsity = T, keepx = 3.1, keepy = 10)
fit %>% str
summary(fit)
plot(fit, "Yjoint_gr", i=1,j=3)

fit <- o2m(X, Y, r, rx, ry, sparsity = T, groupx=groupx, groupy=groupy, keepx = 3.1, keepy = 3)
fit %>% str
summary(fit)
loadings(fit, "Yjoint_gr")
plot(fit, "Yjoint_gr", i=1,j=3)
