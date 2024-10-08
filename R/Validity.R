validity <- function(y, X, reg, v = 0, s = 0,l=0,r=0,N=1000,text=0) {
  #' Test for the Validity of Regression Models 
  #'
  #' validity tests if a regression is valid. An own regression can be specified for this purpose, otherwise OLS is the default option. The function returns the corresponding p-value and t-value.
  #' @param y Dependent variable, vector of (nx1)
  #' @param X Regressors, can be either a vector (nx1) or a matrix (nxm) if there is more than one regressor
  #' @param reg (optional, default is OLS) It is possible to define your own regression function. The output of the function must be a vector of the fitted values. The input must be the dependent variable and the independent variables. See the function "reg" in the readme file as an example. Link: https://github.com/floschuetze/Validity
  #' @param v (optional, default is v = 0) Either 0 or 1. v = 0: Validity-Test from Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #' 
  #' v = 1: Test from Stute,W., 1997, Nonparametric model checks for regression. The Annals of Statistics 25, pp. 613 to 641.
  #' @param s (optional, default is s = 0) Either 0 or 1. s = 0: Cramer-von-Mises; s = 1: Kolmogorov-Smirnov; This does only matter for v = 1
  #' @param l (optional, default is l = 0) Controls how much weight is given to the left side of [0,1] of the beta distribution
  #' @param r (optional, default is r = 0) Controls how much weight is given to the right side of [0,1] of the beta distribution
  #' @param N (optional, default is N = 1000) The number of bootstrap replicates. Code may run faster if N<1000
  #' @param text (optional, default is text = 1) Either 0 or 1. If 0, output text is suppressed
  
  #' @details The "validity" function represents the fundamental component of the package. It performs a statistical test with the objective of determining the validity of a regression model that requires a dependent variable, which is represented as a vector of length "n". Additionally, it requires the regressors, which could be a vector or a matrix of length "m*n". By default, an OLS regression is run using a constant and all regressors to explain y. 
  #' @return The function returns the hypothesis to be tested and the associated p-value.
  #' @examples It is evident that the model "y=constant+beta*x" is valid.
  #' x <- rnorm(100, mean=0, sd=1)
  #' y <- 10*x
  #' validity(y,x)
  #' 
  #' It is evident that the model "y=constant+beta*x" is not valid.
  #' x <- rnorm(100, mean=0, sd=1)
  #' y <- x^3-x^2
  #' validity(y,x)
  #' @references 
  #'Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'
  #'Stute,W., 1997, Nonparametric model checks for regression. The Annals of Statistics 25, pp. 613 to 641.
  #' @author Florian Schuetze 26.08.2024
  #'@export
  
  #  Test for the validity of a regression model.
  # ---
  # Input:
  # y (n x 1) - dependent variable
  # X (n x m) - regressors
  # reg (@)   - regression function delifering fitted values
  # v (1 x 1) - version of the test
  # v == 0: Frahm-Schütze
  # v == 1: Stute
  # s (1 x 1) - test statistic (only relevant if v == 1)
  # s == 0: Cramer-von-Mises
  # s == 1: Kolmogorov-Smirnov
  # The choice of s does not matter for v = 0.
  # l (1 x 1) - weight on the left of [0,1]
  # r (1 x 1) - weight on the right of [0,1]
  # The choice of l and r does not matter for v = 1.
  # ---
  # Output:
  # p (1 x 1) - p-value
  # t (1 x 1) - test statistic
  #Florian Schütze 08/2024 based on Gabriel Frahm 02.06.2024
  
  if (missing(v)) v <- 0
  if (missing(s)) s <- 0
  if (missing(l)) l <- 0
  if (missing(r)) r <- 0
  if (missing(N)) N <- 1000
  if (missing(text)) text <- 1
  
  if (missing(reg)) reg <-function(y, X) {
    d<-data.frame(X,y1=y)
    model <- lm(y1 ~ ., data = d)
    return(model$fitted.values)
  }
  
  X<-as.matrix(X)
  n <- nrow(X)
  m <- ncol(X)
  f <- reg(y, X)
  Data <- (y - f) - mean(y - f)
  E <- matrix(sample(Data, n * N, replace = TRUE), nrow = n, ncol = N)
  
  rm(Data)
  Y <- f + E
  F <- apply(Y, 2, function(y) reg(y, X))
  Z1 <- matrix(rep(X, m), nrow = n, ncol = m * n)
  interleaved <- as.vector(t(X))[order(rep(1:nrow(X), each=m), rep(1:m, times=n))]
  Z2 <- matrix(rep(interleaved, n), nrow = n, byrow = TRUE)
  A<-Z1-Z2
  
  rm(Z1,Z2,interleaved)
  hhh<-cbind(rep(1:n, times = n), rep(seq(1,ncol(A),m), each = n))
  O<-matrix(0,nrow=n*n,ncol=1)
  hhh1<-hhh
  for (i in 1:m){
    hhh1[,2] <- hhh[, 2] + i-1
    O <- A[hhh1]
    if(i>1){
      C<-cbind(C,O)
    }else{
      C<-O
    }}
  B<-as.matrix(C*-1)
  
  rm(O,C,hhh,hhh1)
  k <- matrix(rep((1:n), each = n), nrow = n, byrow = TRUE)
  beta_pdf_values <- dbeta(((1:n) / n), shape1 = r + 1, shape2 = l + 1)
  w <- matrix(rep(beta_pdf_values, each = n), nrow = n, byrow = TRUE)
  K <- array(rep(1:n, n), dim = c(n, n, N))
  W <- array(rep(w, each = N), dim = c(dim(w), N))
  
  rm(beta_pdf_values)
  if (v == 0){
    C <- matrix(0, nrow = n, ncol = n)
    C <- sapply(1:n, function(i) seq(i, n*n, n))
    C <- matrix(as.vector(C), nrow = length(C), ncol = 1)
    if (m>1){
      kkkk<-matrix(sqrt(rowSums(B[C, ]^2)), nrow = n,ncol=n, byrow = TRUE)
    }else{
      kkkk<-matrix(sqrt(B[C,1 ]^2), nrow = n,ncol=n, byrow = TRUE)
    }
    I <- apply(kkkk, 2, order)
    
    rm(kkkk)
    d <- y-f
    e <- matrix(d[I], nrow = n, ncol = n)
    g <- apply(e, 2, cumsum)
    t <- mean(apply((w * g * g / k), 2, mean))
    R <- array((Y - F), dim = c(n, 1, N))
    E <- array(R[matrix(I, n^2, 1),,], dim = c(n, n, N))
    G <- apply(E, c(2, 3), cumsum)
    T <- apply(apply(((W * G * G) / K), c(1, 3), mean),2,mean)
  } else if (v == 1){
    I <- array(rowSums(B <= 0) == m, dim = c(n, n))
    e <- matrix(rep(y - f, n), nrow = length(y), ncol = n)
    e[I] <- 0
    Ex <- array((Y - F), dim = c(nrow(Y - F), 1, ncol(Y - F)))
    E <- aperm(array(rep(Ex, each = n), dim = c(dim(Ex)[1], n, dim(Ex)[3])), c(2, 1, 3))
    E[array(I, dim = c(nrow(I), ncol(I), N))]=0
    if (s == 0) {
      t <- sum(colSums(e)^2) / n^2
      T <- colSums(apply(E, c(3), function(x) colSums(x))^2, dims = 1)/n^2
    } else if (s == 1) {
      t <- max(abs(colSums(e)))
      T<-apply(abs(apply(E, c(3), function(x) colSums(x))),2,max)
    } else {
      stop("Test statistic not properly specified!")
    }
  } else {
    
    print("Error:Test version not properly specified!")
  }
  
  if (t <= .Machine$double.eps){
    t = 0
    p = 1
  }else{
    p <- 1 - sum(T <= t) / N}
  if (text ==1){
    message1 <- "The validity test was successfully completed."
    message2 <- "H0: The regression model is valid."
    message3 <- "H1: The regression model is invalid."
    cat(message1, "\n")
    cat(message2, "\n")
    cat(message3, "\n")
    cat("p-Value:", p, "\n")
  }
  return(list(p_value = p))
}

reg<-function(y, X) {
  #' The default OLS regression
  #'
  #' This function generates fitted values based on a set of regressors X and the independent variable y.
  #' @param y Dependent variable, vector of (nx1)
  #' @param X Regressors, can be either a vector (nx1) or a matrix (nxm) if there is more than one regressor
  #' @return The function returns the fitted values based on the input and the use of standard OLS (including a constant)
  #' @examples 
  #' x <- rnorm(100, mean=0, sd=1)
  #' y <- 10*x
  #' reg(y,x)
  #' @author Florian Schuetze 26.08.2024
  #' @export
  
  d<-data.frame(X,y1=y)
  model <- lm(y1 ~ ., data = d)
  return(as.vector(model$fitted.values))
}

Example1 <- function(c,n,tau,N,v) {
  #' The first example of Frahm (2024)
  #'
  #' This function replicates the first example of Frahm (2024)
  #' @param c A constant
  #' @param n The number of observations
  #' @param tau (optional, default is tau = 1) tau is the standard deviation of the error e in the regression Y = g(X) + e
  #' @param N (optional, default is N = 1000) N is the number of bootstrap replications in the validity function
  #' @param v (optional, default is v = 0) Either 0 or 1. v = 0: Validity-Test from Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #' 
  #' v = 1: Test from Stute,W., 1997, Nonparametric model checks for regression. The Annals of Statistics 25, pp. 613 to 641.
  #' @return The function returns the p-value of the validity test.
  #' @examples 
  #' Example1(0.25,100)
  #' @author Florian Schuetze 26.08.2024
  #' @references 
  #'Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'@export
  
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

Example2 <- function(c,n,tau,N,v) {
  #' The second example of Frahm (2024)
  #'
  #' This function replicates the second example of Frahm (2024)
  #' @param c A constant
  #' @param n The number of observations
  #' @param tau (optional, default is tau = 1) tau is the standard deviation of the error e in the regression Y = g(X) + e
  #' @param N (optional, default is N = 1000) N is the number of bootstrap replications in the validity function
  #' @param v (optional, default is v = 0) Either 0 or 1. v = 0: Validity-Test from Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #' 
  #' v = 1: Test from Stute,W., 1997, Nonparametric model checks for regression. The Annals of Statistics 25, pp. 613 to 641.
  #' @return The function returns the p-value of the validity test.
  #' @examples 
  #' Example2(0.25,100)
  #' @author Florian Schuetze 26.08.2024
  #' @references 
  #'Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'@export
  
  if (missing(tau)) tau <- 1
  if (missing(N)) N <- 1000
  if (missing(v)) v <- 0
  
  x <- rnorm(n, mean=0, sd=1)
  y<-(-1)+x+c*((x^2)-1)+rnorm(n, mean=0, sd=tau)
  p<-validity(y,x,v=v,N=N,text=0)$p_value
  return(p)
}

Example3 <- function(c,n,tau,N,v) {
  #' The third example of Frahm (2024)
  #'
  #' This function replicates the third example of Frahm (2024)
  #' @param c A constant
  #' @param n The number of observations
  #' @param tau (optional, default is tau = 1) tau is the standard deviation of the error e in the regression Y = g(X) + e
  #' @param N (optional, default is N = 1000) N is the number of bootstrap replications in the validity function
  #' @param v (optional, default is v = 0) Either 0 or 1. v = 0: Validity-Test from Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #' 
  #' v = 1: Test from Stute,W., 1997, Nonparametric model checks for regression. The Annals of Statistics 25, pp. 613 to 641.
  #' @return The function returns the p-value of the validity test.
  #' @examples 
  #' Example3(0.25,100)
  #' @author Florian Schuetze 26.08.2024
  #' @references 
  #'Frahm, G., 2024, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'@export
  
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