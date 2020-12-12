# Fubini norm
normF <- function(x) sum(x^2)
# 2-norm
norm2 <- function(x) sqrt(sum(x^2))
# 1-norm
norm1 <- function(x) sum(abs(x))
# update u
update_u <- function(yk, X, v)
{
  X_u <- kronecker(v, X)
  n <- nrow(X)
  lar <- Lars(yk, X_u)
  u_ <- as.matrix(lar$Beta[[which.min(lar$Cp),]], ncol = 1)
  d <- norm2(X%*%u_) / sqrt(n)
  u <- u_ / d
  list(d = d, u = u)
}
# update v
update_v <- function(yk, X, u)
{
  X_v <- kronecker(diag(1,q), X%*%u)
  lar <- Lars(yk, X_v)
  v_ <- as.matrix(lar$Beta[[which.min(lar$Cp),]], ncol = 1)
  d <- norm2(v_)
  v <- v_ / d
  list(d = d, v = v)
}

# Lars-CURE
#' This is some descriptio of this function.
#' @title CURE using Lars algorithm 
#' 
#' @description CURE using Lars algorithm 
#' 
#' @param yk the vector of current response variable
#' @param X design matrix
#' @param d0 initial value of d
#' @param u0 initial value of u
#' @param v0 initial value of v
#' @param lambda tuning parameter
#' @param mu parameter
#' @param eps iteration residual bound
#' @param max_iter maximal number of iteration
#'
#' @return a list including estimated coefficients (d,u,v) and iteration
#' 
#' @export

CURE <- function(yk, X, d0, u0, v0, lambda, mu, eps, max_iter)
{
  iter <- 0
  res <- 1
  d <- d0
  u <- u0
  v <- v0
  while ((res >= eps) & (iter < max_iter))
  {
    Ct <- d * u %*% t(v)
    # update u
    d <- update_u(yk, X, v)$d
    u <- update_u(yk, X, v)$u
    # update v
    d <- update_v(yk, X, v)$d
    v <- update_v(yk, X, v)$v
    C <- d * u %*% t(v)
    res <- normF(C-Ct) / normF(Ct)
    iter = iter + 1
  }
  list(d = d, u = u, v = v, iter = iter)
}