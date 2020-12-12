# Lars update
#' This is some descriptio of this function.
#' @title update step for Lars algorithm for linear regression
#' 
#' @description update step for Lars algorithm for linear regression
#' 
#' @param y a vector of response variable
#' @param X design matrix
#' @param bt current coefficient vector
#'
#' @return coefficient vector after update
#' 
#' @examples
#' \dontrun{
#' y < rnorm(100)
#' X < matrix(runif(500), nrow = 100)
#' bt < rnorm(5)
#' Lars_update(y, X, bt)
#' } 
#' @export

Lars_update <- function(y, X, bt)
{
  p <- ncol(X)
  # current correlations
  c_ <- t(X) %*% (y-X%*%bt)
  C_ <- max(abs(c_))
  # active set
  A <- which(abs(abs(c_)-C_) < 1e-5, arr.ind = TRUE)[,1]
  s <- sign(c_)
  XA <- as.matrix(X[,A]) %*% diag(s[A], nrow = length(A))
  gA <- t(XA) %*% XA
  m <- length(A)
  A1 <- matrix(1, nrow = m, ncol = 1)
  AA <- as.numeric(sqrt(t(A1) %*% solve(gA) %*% A1))
  wA <- AA * solve(gA) %*% A1
  # equiangular vector
  uA <- XA %*% wA
  # inner product vector
  a <- t(X) %*% uA
  if (length(A) < p)
  {
    index <- setdiff(c(1:p), A)
    gamma_set <- c((C_-c_[index])/(rep(AA,length(index))-a[index]), (C_+c_[index])/(rep(AA,length(index))+a[index]))
    gamma_set <- gamma_set[gamma_set > 0]
    # step size
    gamma_ <- min(gamma_set)
  }
  else
    gamma_ <- 1
  w <- matrix(0, nrow = p, ncol = 1)
  w[A,] <- diag(s[A], nrow = length(A)) %*% wA   
  # update
  b <- bt + gamma_*w
  b
}

# Lars
#' This is some descriptio of this function.
#' @title Lars algorithm for linear regression
#' 
#' @description Lars algorithm for linear regression
#' 
#' @param y a vector of response variable
#' @param X design matrix
#'
#' @return a list including estimated coefficients and Cps
#' 
#' @examples
#' \dontrun{
#' y <- rnorm(100)
#' X <- matrix(runif(500), nrow = 100)
#' Lars(y, X)
#' } 
#' 
#' @export

Lars <- function(y, X)
{
  n <- nrow(X)
  p <- ncol(X)
  bt <- matrix(0, nrow = p, ncol = 1)
  res <- list(Beta = list(), Cp = list())
  for (t in 1:p)
  {
    bt <- Lars_update(y, X, bt)
    res$Beta <- c(res$Beta, list(bt))
  }
  sigma <- sum((y - X%*%res$Beta[[p]])^2) / (n - p - 1)
  for (i in 1:p)
  {
    RSS <- sum((y - X%*%res$Beta[[i]])^2)
    cp <- RSS/sigma - (n-2*i) 
    res$Cp <- c(res$Cp, list(cp))
  }
  res
}
