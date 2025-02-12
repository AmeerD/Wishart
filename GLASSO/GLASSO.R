library(tidyverse)
library(patchwork)
library(mvtnorm)
library(purrr)
library(readxl)
library(glasso)
library(reshape2)
library(Matrix)
library(latex2exp)

set.seed(12345)

## Set input parameters
n <- 250
p <- 10
pblock <- floor(p/2)-1
Q <- diag(p)
Q <- bdiag(matrix(0.5, nrow=pblock, ncol=pblock) + 0.5*diag(pblock),
           matrix(0.25, nrow=pblock, ncol=pblock) + 0.75*diag(pblock),
           diag(p-2*pblock))
Sig <- solve(Q)

nsim <- 10
rho <- seq(0, 0.1, length.out=1001)

ncv <- 10

## Results array with dimensions for simulation, method, and rho
res <- array(0, dim=c(nsim,2,length(rho)))
rmin <- array(0, dim=c(nsim,2))
zeros <- array(0, dim=c(nsim,length(rho)))

dmat <- rmvnorm(n, sigma=as.matrix(Sig))
  
H <- eigen(diag(n) - matrix(1/n, nrow=n, ncol=n))$vectors[,-1]
  
for (i in 1:nsim) {
  ## Cross-validation
  ### Shuffle data
  dmat <- dmat[sample(n,n),]
  
  ### 10-fold Cross-validation
  cv.idx <- 1 + (1:n %% ncv)
  for (j in 1:ncv) {
    Str <- cov(dmat[which(cv.idx != j),])
    Ste <- cov(dmat[which(cv.idx == j),])
    
    glres <- glassopath(Str, rho, penalize.diagonal=F, trace=0)
    
    res[i,1,] <- res[i,1,] - apply(glres$wi, 3, function(x) {determinant(x)$modulus - sum(diag(x %*% Ste))})
  }
  
  rmin[i,1] <- which(res[i,1,] == min(res[i,1,]))
  
  ## Data thinning
  mhat <- rep(1,n) %*% t(colMeans(dmat))
  w <- crossprod(dmat - mhat)
  e <- eigen(w)
  dprime <- pracma::randortho(n-1)[,1:p] %*% diag(sqrt(e$values)) %*% t(e$vectors) 
  
  ### 10-fold Cross-validation
  cv.idx <- 1 + (1:(n-1) %% ncv)
  for (j in 1:ncv) {
    Str <- crossprod(dprime[which(cv.idx != j),])/length(which(cv.idx != j))
    Ste <- crossprod(dprime[which(cv.idx == j),])/length(which(cv.idx == j))
    
    glres <- glassopath(Str, rho, penalize.diagonal=F, trace=0)
    
    res[i,2,] <- res[i,2,] - apply(glres$wi, 3, function(x) {determinant(x)$modulus - sum(diag(x %*% Ste))})
  }
  
  rmin[i,2] <- which(res[i,2,] == min(res[i,2,]))
  
  glfinal <- glassopath(w/(n-1), rho, penalize.diagonal=F, trace=0)
  zeros[i,] <- apply(glfinal$wi, 3, function(x){sum(abs(x[lower.tri(x)]) <= 1e-8)})
}

lplot <- melt(res) %>%
  mutate(Var3 = rho[Var3]) %>%
  rename(sim = Var1, method = Var2, rho = Var3) %>%
  mutate(method = ifelse(method == 1, "Sample splitting", "Data thinning"),
         sim = paste0(method, "-", sim)) %>%
  # filter(method == "Sample splitting") %>%
  ggplot(aes(x=rho, y=value, group=sim, colour=method)) +
  geom_line(alpha=0.5) +
  theme(legend.position="right") +
  labs(colour="") +
  xlab(TeX(r"($\lambda$)")) +
  ylab("Test set negative log-likelihood")

lplot

zeros.final <- array(0, dim=c(nsim,2))

for (i in 1:nsim) {
  zeros.final[i,1] <- zeros[i,rmin[i,1]]
  zeros.final[i,2] <- zeros[i,rmin[i,2]]
}

print(sum(abs(Q[lower.tri(Q)]) <= 1e-8))
print(mean(zeros.final))

ggsave("Ameer's Notes/Wishart Mini Paper/glasso.jpg", lplot, width=6, height=3)






