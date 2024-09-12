
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Validity

<!-- badges: start -->
<!-- badges: end -->

The goal of the “Validity” package is to implement the method of the
paper “A Test for the Validity of Regression Models”.

Citation: *Frahm, Gabriel, A Test for the Validity of Regression Models
(2024). Available on SSRN:* <https://ssrn.com/abstract=4610329> or
<http://dx.doi.org/10.2139/ssrn.4610329>

## Disclaimer

This package is a pre-alpha version. The execution time of the
“validity” function in particular is still in need of improvement,
especially if the aim is to run a comprehensive simulation. At present,
it takes about 2.5 seconds to perform the validity test with n=100. With
n\>=1000 or more than one regressor you also need a lot of RAM. In these
cases, it is recommended to adjust the default bootstrap of N=1000 in
the “validity” function (e.g. to N=100). Especially since R uses only
one CPU core in standard mode.

Comments, help, and suggestions for improvement are greatly appreciated
at [schuetze@hsu-hh.de](schuetze@hsu-hh.de).

## Installation

You can install the development version of Validity from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("floschuetze/Validity")
```

## Primary Function

The “validity” function represents the fundamental component of the
package. It performs a statistical test with the objective of
determining the validity of a regression model that requires a dependent
variable, which is represented as a vector of length “n”. Additionally,
it requires regressors, which could be a vector or a matrix of length
“m\*n”. By default, an OLS regression is run using a constant and all
regressors to explain y.

For additional input options, please refer to the help file for the
validity function. This is particularly the case when a regression model
other than ordinary least squares (OLS) is to be employed.

The following example illustrates the application of the validity
function in the context of a standard OLS regression, whereby the null
hypothesis (H0) is rejected.

``` r
library(Validity)
#It is evident that the model "y=constant+beta*x" is not valid.
x <- rnorm(100, mean=0, sd=1)
y <- x^3-x^2 # true model
validity(y,x) 
#> The validity test was successfully completed. 
#> H0: The regression model is valid. 
#> H1: The regression model is invalid. 
#> p-Value: 0
#> $p_value
#> [1] 0
```

## Additional optional functions within the package

### Default Regression function

It is possible to specify your own regression function and pass it to
the “validity” function. It is important that your custom function has
as output the estimated values based on the custom regression for y as a
vector.

The following regression function is the default simple OLS regression
used in the function “validity”.

``` r
reg<-function(y, X) {
  d<-data.frame(X,y1=y)
  model <- lm(y1 ~ ., data = d)
  return(as.vector(model$fitted.values))
}
```

### Example 1

The following function can be used to replicate the first example from
the “A Test for the Validity of Regression Models” paper.

``` r
#c and n must be specified; the function returns the p-value; tau, N and v is optional
Example1 <- function(c,n,tau,N,v) {
  if (missing(tau)) tau <- 1
  if (missing(N)) N <- 1000
  if (missing(v)) v <- 0
  x <- rnorm(n, mean=0, sd=1)
  y <- rnorm(n, mean=0, sd=tau)
  y[which(x>0)]<- y[which(x>0)]+c
  y[which(x<=0)]<- y[which(x<=0)]-c
  p<-validity(y,x,v=v,N=N,text=0)$p_value
  return(p)
}
```

### Example 2

The following function can be used to replicate the second example from
the “A Test for the Validity of Regression Models” paper.

``` r
#c and n must be specified; the function returns the p-value; tau, N and v is optional
Example2 <- function(c,n,tau,N,v) {
  if (missing(tau)) tau <- 1
  if (missing(N)) N <- 1000
  if (missing(v)) v <- 0
  x <- rnorm(n, mean=0, sd=1)
  y<-(-1)+x+c*((x^2)-1)+rnorm(n, mean=0, sd=tau)
  p<-validity(y,x,v=v,N=N,text=0)$p_value
  return(p)
}
```

### Example 3

The following function can be used to replicate the third example from
the “A Test for the Validity of Regression Models” paper.

``` r
#c and n must be specified; the function returns the p-value; tau, N and v is optional
Example3 <- function(c,n,tau,N) {
  if (missing(tau)) tau <- 1
  if (missing(N)) N <- 1000
  if (missing(v)) v <- 0
  require(MASS)
  sigma<-rbind(c(1,0.5), c(0.5,1))
  mu<-c(0, 0) 
  LK<-as.matrix(mvrnorm(n=n, mu=mu, Sigma=sigma))
  y<-0.25*LK[,1]+0.75*LK[,2]+c*LK[,1]*LK[,2]+rnorm(n, mean=0, sd=tau)
  p<-validity(y,LK,v=v,N=N,text=0)$p_value
  return(p)
}
```
