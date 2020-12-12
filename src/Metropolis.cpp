# include <Rcpp.h>
# include <math.h>
using namespace Rcpp;

//' @title A Metropolis sampler using Rcpp
//' 
//' @description A Metropolis sampler using Rcpp
//' 
//' @param x0 the initial point of the series
//' @param sigma the standard deviation
//' @param N the sample size
//' 
//' @return a random sample of size \code{N}
//' 
//' @examples
//' \dontrun{
//' x0 = 0;
//' sigma = 2;
//' N = 1000;
//' x = Metropolis(x0, sigma, N)
//' }
//' 
//' @export
// [[Rcpp::export]]

NumericVector Metropolis(double x0, double sigma, int N) 
{
  NumericVector x(N);
  
  x[0] = x0;
  for (int i = 1; i < N; i++)
  {
    double y; 
    double u; 
    y = rnorm(1, x[i-1], sigma)[0];
    u = runif(1)[0];
    if (u <= (exp(-abs(y)) / exp(-abs(x[i-1]))))
      x[i] = y;
    else 
    x[i] = x[i-1];
  }
  return (x);
}
