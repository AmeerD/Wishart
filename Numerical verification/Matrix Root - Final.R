library(tidyverse)
library(patchwork)
library(mvtnorm)

set.seed(123)

nsim <- 10000

n <- 3
p <- 5
Sigma <- toeplitz(1/(1:p))
# for (i in 1:p) {
#   for (j in 1:p) {
#     Sigma[i,j] = 1/(1+abs(i-j))
#   }
# }

res <- array(NA, dim=c(2,n,p,nsim))

for (ns in 1:nsim) {
  x <- rmvnorm(n, sigma=Sigma) #matrix(rnorm(n*p),n,p)
  w <- crossprod(x)
  
  e <- eigen(w)
  e$values <- e$values[1:n]
  e$vectors <- e$vectors[,1:n]
  
  ## DV' Root
  res[1,,,ns] <- diag(sqrt(e$values)) %*% t(e$vectors)

  ## QDV' Root
  res[2,,,ns] <- pracma::randortho(n) %*% diag(sqrt(e$values)) %*% t(e$vectors) 
}

res2 <- res[,,,] %>% 
  reshape2::melt() %>%
  rename(method = Var1, i=Var2, j=Var3, idx=Var4)

p1 <- res2 %>% 
  filter(method == 1) %>%
  mutate(i = paste0("i=", i), j = paste0("j=", j)) %>%
  ggplot(aes(x=value)) +
  geom_histogram(aes(y=..density..), bins=100) +
  facet_grid(i~j) + 
  stat_function(fun=dnorm, args=list(mean=0, sd=1), col='seagreen', linewidth=0.5) +
  ggtitle("Standard matrix square root") +
  theme(text = element_text(size = 7.5))  +
  xlim(c(-5,5)) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4))

p2 <- res2 %>% 
  filter(method == 2) %>%
  mutate(i = paste0("i=", i), j = paste0("j=", j)) %>%
  ggplot(aes(x=value)) +
  geom_histogram(aes(y=..density..), bins=100) +
  facet_grid(i~j) +
  stat_function(fun=dnorm, args=list(mean=0, sd=1), col='seagreen', linewidth=0.5) +
  ggtitle("Algorithm 1 matrix square root") +
  theme(text = element_text(size = 7.5)) +
  xlim(c(-5,5)) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4))

ggsave("Ameer's Notes/Wishart Mini Paper/verification.jpg", p1+p2, width=8, height=3)

###################################################

set.seed(123)

nsim <- 10000

n <- 3
p <- 5
mu <- 1:5
Sigma <- toeplitz(1/(1:p))
# for (i in 1:p) {
#   for (j in 1:p) {
#     Sigma[i,j] = 1/(1+abs(i-j))
#   }
# }

H <- eigen(diag(n) - matrix(1/n, nrow=n, ncol=n))$vectors[,1:2]

ares <- array(NA, dim=c(n,p,nsim))

for (ns in 1:nsim) {
  x <- rmvnorm(n, mean=mu, sigma=Sigma) #matrix(rnorm(n*p),n,p)
  mhat <- rep(1,n) %*% t(colMeans(x))
  w <- crossprod(x - mhat)
  
  e <- eigen(w)
  e$values <- e$values[1:(n-1)]
  e$vectors <- e$vectors[,1:(n-1)]
  
  ## QDV' Root
  ares[,,ns] <- mhat + H %*% pracma::randortho(n-1) %*% diag(sqrt(e$values)) %*% t(e$vectors) 
}

ares2 <- ares[,,] %>% 
  reshape2::melt() %>%
  rename(i=Var1, j=Var2, idx=Var3)

sfdat <- expand_grid(i=1:3, j=1:5, x=seq(-5,5,by=0.01)) %>%
  mutate(x = x + j) %>%
  mutate(dens = dnorm(x, mean=j, sd=1)) %>% 
  mutate(i = paste0("i=", i), j = paste0("j=", j))

q1 <- ares2 %>% 
  mutate(i = paste0("i=", i), j = paste0("j=", j)) %>%
  ggplot(aes(x=value)) +
  geom_histogram(aes(y=..density..), bins=100) +
  geom_line(data=sfdat, aes(x=x, y=dens), linewidth=0.5, colour="seagreen") +
  facet_grid(i~j) + 
  ggtitle("Algorithm 2 matrix square root") +
  theme(text = element_text(size = 7.5))  #+
  # xlim(c(-5,5)) +
  # scale_x_continuous(breaks = c(-4,-2,0,2,4))

plot(q1)

ares2 %>% group_by(i,j) %>% summarise(m=mean(value, na.rm=T), v=var(value, na.rm=T))

ggsave("Ameer's Notes/Wishart Mini Paper/appvar.jpg", q1, width=4, height=3)
