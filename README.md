
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Validity

<!-- badges: start -->
<!-- badges: end -->

The goal of the “Validity” package is to implement the method of the
paper “A Test for the Validity of Regression Models”.

Citation: *Frahm, Gabriel, A Test for the Validity of Regression Models
(October 23, 2023). Available on SSRN:*
<https://ssrn.com/abstract=4610329> or
<http://dx.doi.org/10.2139/ssrn.4610329>

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
it requires the regressors, which could be a vector or a matrix of
length “m\*n”. By default, an OLS regression is run using a constant and
all regressors to explain y.

For additional input options, please refer to the help file for the
validity function. This is particularly the case when a regression model
other than ordinary least squares (OLS) is to be employed.

This is a simple example that shows you how to use the validity function
where a standard OLS regression is valid.

``` r
library(Validity)
#It is evident that the model "y=constant+beta*x" is valid.
x <- rnorm(100, mean=0, sd=1)
y <- 10*x
validity(y,x)
#> The validity test was successfully completed. 
#> H0: The model is considerd to be valid. 
#> H1: The model is not considered to be valid. 
#> t-Value: 0 
#> p-Value: 1
#> $t_value
#> [1] 0
#> 
#> $p_value
#> [1] 1
```

This is a simple example that shows you how to use the validity function
where a standard OLS regression is not valid.

``` r
#It is evident that the model "y=constant+beta*x" is not valid.
x <- rnorm(100, mean=0, sd=1)
y <- x^3-x^2
validity(y,x)
#> The validity test was successfully completed. 
#> H0: The model is considerd to be valid. 
#> H1: The model is not considered to be valid. 
#> t-Value: 140.2125 
#> p-Value: 0
#> $t_value
#> [1] 140.2125
#> 
#> $p_value
#> [1] 0
```

## Additional optional functions within the package

### Default Regression function

It is possible to specify your own regression function and pass it to
the “validity” function. It is important that your custom function has
as output the estimated values based on the custom regression for y. The
default regression function is a simple OLS regression.

``` r
reg<-function(y, X) {
  d<-data.frame(X,y1=y)
  model <- lm(y1 ~ ., data = d)
  return(model$fitted.values)
}
```

### Default Regression function

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
