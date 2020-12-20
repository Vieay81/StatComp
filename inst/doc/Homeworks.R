## -----------------------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
# estimated by OLS
lm.D9 <- lm(weight ~ group)
summary(lm.D9)

## -----------------------------------------------------------------------------
# residual plots
plot(lm.D9)

## ----iris---------------------------------------------------------------------
xtable::xtable(head(iris))

## -----------------------------------------------------------------------------
set.seed(1)
a <- 2
b <- 2
n <- 1000
u <- runif(n)
x <- 2*(1-u)^(-1/2)
hist(x, prob = TRUE, main = expression(f(x)==8*x^{-3}))
y <- seq(2, 100, .01)
lines(y, 8*y^(-3))

## -----------------------------------------------------------------------------
set.seed(2)
n <- 10000
u <- matrix(runif(3*n,-1,1), nrow = n)
x <- numeric(n)
L <- apply(abs(u),1,which.max) == 3
coord <- 3-L
for (i in 1:n)
{
  x[i] <- u[i,coord[i]]
}
hist(x, prob = TRUE, main = expression(f[e](x)==frac(3,4)*(1-x^2)))
y <- seq(-1, 1, .01)
lines(y, (3/4)*(1-y^2))

## -----------------------------------------------------------------------------
set.seed(4)
r <- 4
beta <- 2
n <- 1000
u <- runif(n)
y <- beta*(1-u)^(-1/r) - beta
hist(y, prob = TRUE, main = expression(f(y)==frac(64,(2+y)^5)))
x <- seq(0, 10, .01)
lines(x, 64/((2+x)^5))

## -----------------------------------------------------------------------------
set.seed(1)
n <- 10000
u <- runif(n)
x <- (pi/3)*u
theta_ <- (pi/3)*mean(sin(x))
theta <- 1 - cos(pi/3)

## -----------------------------------------------------------------------------
set.seed(2)
n <- 10000
x1 <- runif(n)
theta_1 <- mean(exp(x1))
x2 <- runif(n/2)
theta_2 <- mean(exp(x2)+exp(1-x2))/2
m <- 1000
MC1 <- MC2 <- numeric(m)
for (i in 1:m)
{
  x1 <- runif(n)
  x2 <- runif(n/2)
  MC1[i] <- mean(exp(x1))
  MC2[i] <- mean(exp(x2)+exp(1-x2))/2
}
var1 <- var(MC1)
var2 <- var(MC2)
theta <- exp(1)-1

## -----------------------------------------------------------------------------
x <- seq(1,5,0.01)
g <- function(x) {x^2*exp(-x^2/2)/sqrt(2*pi)}
f1 <- function(x) {2*exp(-x^2/2)/sqrt(2*pi)}
f2 <- function(x) {exp(-x/2)/sqrt(2*pi)}
gs <- c(expression(g(x)==x^2*e^{(-x^2/2)}/sqrt(2*pi)),
        expression(f[1](x)==2*e^{(-x^2/2)}/sqrt(2*pi)),
        expression(f[2](x)==e^{(-x/2)}/sqrt(2*pi)))
plot(x, g(x), type = "l", ylab = "")
lines(x, f1(x), col = "red", lty = 2)
lines(x, f2(x), col = "blue", lty = 3)
legend("topright", legend = gs,
           lty = 1:3, inset = 0.02, col = 1:3)

## -----------------------------------------------------------------------------
n <- 10000
k <- 5
N <- 100
set.seed(2)
theta_ <- numeric(k)
est <- numeric(N)
for (i in 1:N)
{
  for (j in 1:k)
 {
  u <- runif(n/k)
  x <- -log((exp(-1)-1)*(u+j-1)/5+1)
  theta_[j] <- mean((1-exp(-1))/(1+x^2))
 }
 est[i] = mean(theta_)
}
theta <- mean(est)
sd <- sd(est)

## -----------------------------------------------------------------------------
set.seed(3)
n <- 20
mu <- 0
sigma <- 1
N <- 1000
LCL <- UCL <- numeric(N)
for (i in 1:N)
{
  x <- rlnorm(n, 0, 1)
  LCL[i] <- mean(log(x)) - qt(0.975,n-1)*sd(log(x))/sqrt(n)
  UCL[i] <- mean(log(x)) + qt(0.975,n-1)*sd(log(x))/sqrt(n)
}
ECP <- mean(LCL <= mu & UCL >= mu)

## -----------------------------------------------------------------------------
set.seed(4)
n <- 20
mu <- 2
N <- 1000
LCL <- UCL <- numeric(N)
for (i in 1:N)
{
  x <- rchisq(n, 2)
  LCL[i] <- mean(x) - qt(0.975,n-1)*sd(x)/sqrt(n)
  UCL[i] <- mean(x) + qt(0.975,n-1)*sd(x)/sqrt(n)
}
ECP <- mean(LCL <= mu & UCL >= mu)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 100
m <- 10000
a <- c(1, 5, 20) # alternatives
v <- c(1, 10, 100) # alternatives
alpha <- 0.05
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2)/((n+1)*(n+3))))
sk <- function(x)
{
  x_ <- mean(x)
  m3 <- mean((x-x_)^3)
  m2 <- mean((x-x_)^2)
  m3/m2^1.5
}
power_beta <- power_t <- numeric(length(a))
for (i in 1:length(a))
{
  # Beta(a,a)
  sktests_beta <- replicate(m, expr = {
    x <- rbeta(n, a[i], a[i])
    as.integer(abs(sk(x)) >= cv)
  })
  power_beta[i] <- mean(sktests_beta)
  # t(v)
  sktests_t <- replicate(m, expr = {
    x <- rt(n, v[i])
    as.integer(abs(sk(x)) >= cv)
  })
  power_t[i] <- mean(sktests_t)
}

## -----------------------------------------------------------------------------
set.seed(2)
sigma_1 <- 1
sigma_2 <- 1.5
n <- c(10, 50, 100)
alpha <- 0.055
m <- 10000
count5test <- function(x, y)
{
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  as.integer(max(c(outx, outy)) > 5)
}
power_5 <- power_F <- numeric(length(n))
for (i in 1:length(n))
{
  power_5[i] <- mean(replicate(m, expr = {
  x <- rnorm(n[i], 0, sigma_1)
  y <- rnorm(n[i], 0, sigma_2)
  count5test(x,y)
}))
  pvalues <- replicate(m, expr = {
  x <- rnorm(n[i], 0, sigma_1)
  y <- rnorm(n[i], 0, sigma_2)
  Ftest <- var.test(x, y)
  Ftest$p.value
})
  power_F[i] <- mean(pvalues <= alpha)
}

## -----------------------------------------------------------------------------
library(MASS)
set.seed(3)
# Example 6.8
n <- c(10, 20, 30, 50, 100, 500)
mu0 <- rep(0,2)
Sigma <- matrix(c(1,.5,.5,1), nrow = 2)
d <- length(mu0)
cv <- qchisq(c(1-.975,.975), d*(d+1)*(d+2)/6)
multisk <- function(x)
{
  X <- as.matrix(x)
  n <- dim(X)[1]
  Sigma <- (n-1)*var(X)/n
  S_1 <- solve(Sigma)
  X_ <- as.matrix(apply(X, 2, mean))
  b <- 0
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      b <- b + (t(X[i,]-X_)%*%S_1%*%(X[j,]-X_))^3
    }
  }
  b/n^2
}
p_reject <- numeric(length(n))
m <- 1000
for (i in 1:length(n))
{
  p_reject[i] <- mean(replicate(m, expr = {
    x <- mvrnorm(n = n[i], mu0, Sigma)
    sk <- multisk(x)
    as.integer(n[i]*sk/6 < cv[1] | n[i]*sk/6 > cv[2])
  }))
}

## -----------------------------------------------------------------------------
# Example 6.10
alpha <- .1
n <- 30
m <- 2500
mu1 <- mu2 <- rep(0,2)
Sigma1 <- matrix(c(1,.5,.5,1), nrow = 2)
Sigma2 <- matrix(c(100,.5,.5,100), nrow = 2) 
d <- length(mu0)
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
power_multi <- numeric(N)
cv <- qchisq(c(alpha/2,1-alpha/2), d*(d+1)*(d+2)/6)
for (i in 1:N)
{
  e <- epsilon[i]
  power_multi[i] <- mean(replicate(m, expr = {
    index <- sample(c(1,2), replace = T, size = n, prob = c(1-e,e))
    x <- matrix(nrow = n, ncol = d)
    for (j in 1:n)
    {
      if (index[j] == 1)
        x[j,] <- mvrnorm(1, mu1, Sigma1)
      else
        x[j,] <- mvrnorm(1, mu2, Sigma2)
    }
    sk <- multisk(x)
    as.integer(n*sk/6 < cv[1] | n*sk/6 > cv[2])
  }))
}
plot(epsilon, power_multi, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = alpha, lty = 3)
se <- sqrt(power_multi*(1-power_multi)/m)
lines(epsilon, power_multi+se, lty = 3)
lines(epsilon, power_multi-se, lty = 3)

## -----------------------------------------------------------------------------
data(law, package = "bootstrap")
n <- nrow(law)
cor_ <- cor(law$LSAT, law$GPA)
cor_jack <- numeric(n)
for (i in 1:n) 
{
  cor_jack[i] <- cor(law$LSAT[-i], law$GPA[-i])
}
bias_jack <- (n-1) * (mean(cor_jack)-cor_) 
se_jack <- sd(cor_jack) / sqrt(n) * (n-1)

## ----warning=FALSE------------------------------------------------------------
library(boot)
data(aircondit, package = "boot")
set.seed(15798)
B <- 1000
boot_mean <- function(x, i) mean(x[i])
m <- boot(data = aircondit$hours, statistic = boot_mean, R = B)
CI <- boot.ci(m, type = c("norm", "basic", "perc", "bca"))
print(CI)

## ----warning=FALSE------------------------------------------------------------
set.seed(54964)
data(scor, package = "bootstrap")
B <- 2000
n <- nrow(scor)
Sigma_ <- var(scor)
eigen_ <- eigen(Sigma_, symmetric = T)
theta_ <- eigen_$values[1] / sum(eigen_$values)
Jackknife <- function(data)
{
  theta_jack <- numeric(n)
  for (i in 1:n)
  {
    Sigma_jack <- var(scor[-i,])
    eigen_jack <- eigen(Sigma_jack, symmetric = T)
    theta_jack[i] <- eigen_jack$values[1] / sum(eigen_jack$values)
  }
  theta_jack
}
theta_jack <- Jackknife(data = scor)
bias_jack <- (n-1) * (mean(theta_jack)-theta_)
sd_jack <- sqrt((n-1) * mean((theta_jack-mean(theta_jack))^2))

## ----message=FALSE, warning=FALSE---------------------------------------------
library(lattice)
library(DAAG)
data(ironslag)
n <- length(ironslag$magnetic)
e1 <- e2 <- e3 <- e4 <- matrix(nrow = n*(n-1)/2, ncol = 2)
flag <- 1
for (i in 1:(n-1)) 
{
  for (j in (i+1):n)
  {
    y <- ironslag$magnetic[-c(i,j)]
    x <- ironslag$chemical[-c(i,j)]
    lm1 <- lm(y ~ x, data = ironslag)
    y_ <- predict(lm1, newdata = data.frame(x = ironslag$chemical[c(i,j)]))
    e1[flag,] <- ironslag$magnetic[c(i,j)] - y_

    lm2 <- lm(y ~ x + I(x^2), data = ironslag)
    y_ <- predict(lm2, newdata = data.frame(x = ironslag$chemical[c(i,j)]))
    e2[flag,] <- ironslag$magnetic[c(i,j)] - y_
    
    lm3 <- lm(log(y) ~ x, data = ironslag)
    y_ <- predict(lm3, newdata = data.frame(x = ironslag$chemical[c(i,j)]))
    e3[flag,] <- ironslag$magnetic[c(i,j)] - y_
    
    lm4 <- lm(log(y) ~ log(x), data = ironslag)
    y_ <- predict(lm4, newdata = data.frame(x = ironslag$chemical[c(i,j)]))
    e4[flag,] <- ironslag$magnetic[c(i,j)] - y_
    
    flag <- flag + 1
  }
}
cat("model1 =", mean(apply(e1^2, 1, sum)),
"model2 =", mean(apply(e2^2, 1, sum)),
"model3 =", mean(apply(e3^2, 1, sum)),
"model4 =", mean(apply(e4^2, 1, sum))
)

## -----------------------------------------------------------------------------
counttest <- function(x, y)
{
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  max(c(outx, outy))
}

n1 <- 20
n2 <- 30
mu1 <- 1
mu2 <- 0
sigma1 <- sigma2 <- 1
B <- 1000
m <- 100
reject <- 0
for (i in 1:m)
{
  x1 <- rnorm(n1, mu1, sigma1)
  x2 <- rnorm(n2, mu2, sigma2)
  z <- c(x1, x2)
  t0 <- counttest(x1, x2)
  t <- replicate(B, expr = {
    k <- sample(1:length(z), size = length(x1), replace = F)
    x <- z[k]
    y <- z[-k]
    counttest(x, y)
  })
  p <- mean(c(t0, t) >= t0)
  reject <- reject + (p < .05)
}
reject / m

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) 
{
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data = z, k = k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

# Unequal variances and equal expectations
n1 <- 50
n2 <- 50
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 4
m <- 1e3
k <- 3
R <- 999 
n <- n1 + n2
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k)
{
  boot.obj <- boot(data = z, statistic = Tn, R = R, sim = "permutation", sizes = sizes, k = k)
  ts <- c(boot.obj$t0, boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic = ts[1], p.value = p.value)
}
p.values <- matrix(NA, m ,3)
for(i in 1:m)
{
  x <- matrix(rnorm(n1, mu1, sigma1), ncol = 1)
  y <- matrix(rnorm(n2, mu2, sigma2), ncol = 1)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z, N, k)$p.value
  p.values[i,2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i,3] <- bd.test(x = x, y = y, R = 999, seed = i*1165)$p.va
}
alpha <- .05;
power1 <- colMeans(p.values < alpha)

## -----------------------------------------------------------------------------
# Unequal variances and unequal expectations
n1 <- 50
n2 <- 50
mu1 <- 0
mu2 <- 5
sigma1 <- 1
sigma2 <- 4
m <- 1e3
k <- 3
R <- 999 
n <- n1 + n2
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k)
{
  boot.obj <- boot(data = z, statistic = Tn, R = R, sim = "permutation", sizes = sizes, k = k)
  ts <- c(boot.obj$t0, boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic = ts[1], p.value = p.value)
}
p.values <- matrix(NA, m ,3)
for(i in 1:m)
{
  x <- matrix(rnorm(n1, mu1, sigma1), ncol = 1)
  y <- matrix(rnorm(n2, mu2, sigma2), ncol = 1)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z, N, k)$p.value
  p.values[i,2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i,3] <- bd.test(x = x, y = y, R = 999, seed = i*7864)$p.va
}
alpha <- .05;
power2 <- colMeans(p.values < alpha)

## -----------------------------------------------------------------------------
# Non-normal distributions
n1 <- 50
n2 <- 50
mu1 <- 0
mu2 <- 5
sigma1 <- 1
sigma2 <- 4
m <- 1e3
k <- 3
R <- 999 
n <- n1 + n2
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k)
{
  boot.obj <- boot(data = z, statistic = Tn, R = R, sim = "permutation", sizes = sizes, k = k)
  ts <- c(boot.obj$t0, boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic = ts[1], p.value = p.value)
}
p.values <- matrix(NA, m ,3)
for(i in 1:m)
{
  x <- matrix(rt(n1, df = 1), ncol = 1)
  y <- matrix(c(rnorm(n2/2, mu1, sigma1), rnorm(n2/2, mu2, sigma2)), ncol = 1)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z, N, k)$p.value
  p.values[i,2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i,3] <- bd.test(x = x, y = y, R = 999, seed = i*94654)$p.va
}
alpha <- .05;
power3 <- colMeans(p.values < alpha)

## -----------------------------------------------------------------------------
# Unbalanced samples
n1 <- 50
n2 <- 5
mu1 <- 0
mu2 <- 5
sigma1 <- 1
sigma2 <- 4
m <- 1e3
k <- 3
R <- 999 
n <- n1 + n2
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k)
{
  boot.obj <- boot(data = z, statistic = Tn, R = R, sim = "permutation", sizes = sizes, k = k)
  ts <- c(boot.obj$t0, boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic = ts[1], p.value = p.value)
}
p.values <- matrix(NA, m ,3)
for(i in 1:m)
{
  x <- matrix(rnorm(n1, mu1, sigma1), ncol = 1)
  y <- matrix(rnorm(n2, mu2, sigma2), ncol = 1)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z, N, k)$p.value
  p.values[i,2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i,3] <- bd.test(x = x, y = y, R = 999, seed = i*78324)$p.va
}
alpha <- .05;
power4 <- colMeans(p.values < alpha)

## -----------------------------------------------------------------------------
rw_Metropolis <- function(sigma, x0, N)
{
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N)
  {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-abs(y))/exp(-abs(x[i-1])))
      x[i] <- y
    else
    {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  list(x = x, k = k)
}

set.seed(324)
N <- 2000
sigma <- c(.05, .5, 2, 4)
x0 <- 0
rw1 <- rw_Metropolis(sigma[1], x0, N)
rw2 <- rw_Metropolis(sigma[2], x0, N)
rw3 <- rw_Metropolis(sigma[3], x0, N)
rw4 <- rw_Metropolis(sigma[4], x0, N)
plot(1:N, rw1$x, type = "l", xlab = expression(sigma==0.05))
plot(1:N, rw2$x, type = "l", xlab = expression(sigma==0.5))
plot(1:N, rw3$x, type = "l", xlab = expression(sigma==02))
plot(1:N, rw4$x, type = "l", xlab = expression(sigma==4))
accept <- 1 - c(rw1$k, rw2$k, rw3$k, rw4$k) / N
print(accept)

## -----------------------------------------------------------------------------
GR <- function(psi)
{
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  m <- rowMeans(psi)
  B <- n * var(m)
  w <- apply(psi, 1, "var")
  W <- mean(w)
  v_ <- W*(n-1)/n + B/n
  r_ <- v_ / W
  r_
}

set.seed(324)
N <- 2000
burn <- 200
x0 <- c(-10, -5, 5, 10)
X <- matrix(nrow = 4, ncol = N)
for (i in 1:4)
{
  X[i,] <- rw_Metropolis(sigma[2], x0[i], N)$x
}
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
R_ <- GR(psi)
for (i in 1:4)
  plot(psi[i, (burn+1):N], type = "l", xlab = paste("x0 =", x0[i]), ylab = bquote(psi))
r_ <- numeric(N)
for (i in (burn+1):N)
  r_[i] <- GR(psi[,1:i])
plot(r_[(burn+1):N], type = "l", xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)

N <- 5000
burn <- 500
x0 <- c(-10, -5, 5, 10)
X <- matrix(nrow = 4, ncol = N)
for (i in 1:4)
{
  X[i,] <- rw_Metropolis(sigma[2], x0[i], N)$x
}
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
for (i in 1:4)
  plot(psi[i, (burn+1):N], type = "l", xlab = paste("x0 =", x0[i]), ylab = bquote(psi))
r_ <- numeric(N)
for (i in (burn+1):N)
  r_[i] <- GR(psi[,1:i])
plot(r_[(burn+1):N], type = "l", xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)
index <- min(which(r_[r_ != 0] < 1.2))
rw <- rw_Metropolis(sigma[2], x0[2], index+burn)
plot(1:(index+burn), rw$x, type = "l", xlab = "", main = "Runs until converges")

## -----------------------------------------------------------------------------
k <- c(4:25, 100, 500, 1000)
roots <- numeric(length(k))
for (i in 1:length(k))
{
  f <- function(a) pt(sqrt(a^2*(k[i]-1)/(k[i]-a^2)), df = k[i]-1) - pt(sqrt(a^2*k[i]/(k[i]+1-a^2)), df = k[i])
  roots[i] <- uniroot(f, c(1e-10,sqrt(k[i])-1e-10))$root
}
print(roots)
plot(1:length(k), type = "l", roots, xlab = "", ylab = expression(A(k)), main = "Intersection points")

## ----warning=FALSE------------------------------------------------------------
Eloglike <- function(p = p_t, q = q_t) 
    -(2*nOO*log(1-p-q) + nA.*log(p*(1-p-q)) + nB.*log(q*(1-p-q)) + nAB*log(p*q) + (p_t)/(p_t+2*r_t)*nA.*log(p/(1-p-q)) +   (q_t)/(q_t+2*r_t)*nB.*log(q/(1-p-q)))
likelihood <- function(p_t1, q_t1, r_t1, p_t, q_t, r_t)
  log(((p_t1/p_t)^2+2*(p_t1/p_t)*(r_t1/r_t))^nA. * ((q_t1/q_t)^2+2*(q_t1/q_t)*(r_t1/r_t))^nB. * ((r_t1/r_t)^2)^nOO * (2*(p_t1/p_t)*(q_t1/q_t))^nAB) 

library(stats4)
nA. <- 444
nB. <- 132
nOO <- 361
nAB <- 63
N <- 10000
p_t <- q_t <- r_t <- 1/3
eps <- .Machine$double.eps^.5
likelihood_o <- numeric(N)
for (i in 1:N)
{
  fit <- mle(Eloglike)
  p_t1 <- fit@coef[["p"]]
  q_t1 <- fit@coef[["q"]]
  r_t1 <- 1 - p_t1 - q_t1
  likelihood_o[i] <- likelihood(p_t1, q_t1, r_t1, p_t, q_t, r_t)
  res <- sum(abs(c(p_t1, q_t1, r_t1)-c(p_t, q_t, r_t)) / c(p_t, q_t, r_t))
  if (res < eps)
    break
  p_t <- p_t1
  q_t <- q_t1
  r_t <- r_t1
}
print(list("p" = p_t1, "q" = q_t1, "iterations" = i, "res" = res))
plot(2:i, likelihood_o[2:i], type = "l", xlab = "iteration", ylab = "log-likelihood difference", main = "The difference of log-likelihoods")

## -----------------------------------------------------------------------------
data(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
# for loops
for (i in formulas)
{
  model <- lm(i, data = mtcars)
  print(summary(model))
}
# lapply
lm4 <- lapply(formulas, function(f) lm(f, data = mtcars))
print(lm4)

## -----------------------------------------------------------------------------
set.seed(3)
random <- replicate(100, list(rpois(10, 10), rpois(7, 10)))
locs <- 1:100
p <- sapply(locs, function(i) t.test(random[[2*i-1]], random[[2*i]])$p.value)
print(p)

## -----------------------------------------------------------------------------
data(mtcars)
std <- function(x, x_bar) sum((x-x_bar)^2) / (length(x)-1)
mapva <- function(x)
{
  m <- vapply(x, mean, numeric(1))
  vapply(Map(std, x, m), unlist, numeric(1))
}
mapva(mtcars)

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
dir_cpp <- "../src/"
sourceCpp(paste0(dir_cpp, "Metropolis.cpp"))
rw_Metropolis <- function(x0, sigma, N)
{
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  for (i in 2:N)
  {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-abs(y))/exp(-abs(x[i-1])))
      x[i] <- y
    else
      x[i] <- x[i-1]
  }
  x
}

set.seed(123)
sigma <- 2
x0 <- 0
N <- 2000
xr <- rw_Metropolis(x0, sigma, N)
xcpp <- Metropolis(x0, sigma, N)
qqplot(xr[1001:2000], xcpp[1001:2000], main = "QQ plot")

## -----------------------------------------------------------------------------
ts <- microbenchmark(rw_Metropolis(x0, sigma, N), Metropolis(x0, sigma, N))
summary(ts)

