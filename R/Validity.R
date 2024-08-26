
validity <- function(y, X, reg, v = 0, s = 0,l=0,r=0,N=1000) {
  #' Test for the Validity of Regression Models 
  #'
  #' It tests if a regression is valid. A specific regression can be specified for this, otherwise OLS is the default option.
  #' @param y Dependent variable, vector of (nx1)
  #' @param X Regressors, can be either a vector (nx1) or a matrix (nxm) if there is more than one regressor
  #' @param reg (optional, default is OLS) It is possible to define your own regression function. The output of the function must be the fitted values. The input must be the dependent variable and the independent variables. See the function "reg" in the readme file as an example.
  #' @param v (optional, default is v = 0) v = 0: Validity-Test from Frahm, G., 2023, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #' 
  #' v = 1: Test from Stute,W., 1997, Nonparametric model checks for regression. The Annals of Statistics 25, pp. 613 to 641.
  #' @param s (optional, default is s = 0) s = 0: Cramer-von-Mises; s = 1: Kolmogorov-Smirnov; This does only matter for v = 1
  #' @param l (optional, default is l = 0) Controls how much weight is given to the left side of [0,1] of the beta distribution
  #' @param r (optional, default is r = 0) Controls how much weight is given to the right side of [0,1] of the beta distribution
  #' @param N (optional, default is N = 1000) The number of bootstrap replicates. Code may run faster if N<1000
  
  #' @details The "validity" function represents the fundamental component of the package. It performs a statistical test with the objective of determining the validity of a regression model that requires a dependent variable, which is represented as a vector of length "n". Additionally, it requires the regressors, which could be a vector or a matrix of length "m*n". By default, an OLS regression is run using a constant and all regressors to explain y. 
  #' @return The function returns the hypothesis to be tested. The t-value and p-value are also returned
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
  #'Frahm, G., 2023, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
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
  
  Y <- f + E

  F <- matrix(0, nrow = n, ncol = N)
  for (i in 1:N) {
    F[, i] <- reg(Y[, i], X)
  }
  F <- apply(Y, 2, function(y) reg(y, X))

  Z1 <- matrix(NA, nrow = n, ncol = m * n) 
  for (i in 1:n) {
    Z1[, (m * i - 1):(m * i)] <- X
  }
  
  
  
  repeated_rows <- matrix(rep(X, n), nrow = n, byrow = TRUE)
  
  interleaved <- as.vector(t(X))
  interleaved <- interleaved[order(rep(1:nrow(X), each=m), rep(1:m, times=n))]
  Z2 <- matrix(rep(interleaved, n), nrow = n, byrow = TRUE)
  A<-Z1-Z2
  
  
  
  
  O<-matrix(0,nrow=n,ncol=m)
  C<-matrix(0,nrow=0,ncol=m)
  
  for (ff in seq(1,ncol(A),m)){
    for (j in 1:m){
      O[,j]<-A[,(ff-1+j)]}
    
    C<-rbind(C,O)
  }
  B<-C*-1
  col_vector <- 1:n
  
  k <- matrix(rep(col_vector, each = n), nrow = n, byrow = TRUE)
  
  beta_pdf_values <- dbeta(((1:n) / n), shape1 = r + 1, shape2 = l + 1)
  
  w <- matrix(rep(beta_pdf_values, each = n), nrow = n, byrow = TRUE)
  K <- array(rep(1:n, n), dim = c(n, n, N))
  W <- array(rep(w, each = N), dim = c(dim(w), N))
  
  if (v == 0){
    C<-matrix(0,nrow=0,ncol=1)
    
    for (i in 1:n){
      o<-seq(i,n*n,n)
      o<-as.matrix(o)
      C<-rbind(C,o)
    }
    if (m>1){
      kkkk<-matrix(sqrt(rowSums(B[C, ]^2)), nrow = n,ncol=n, byrow = TRUE)
    }else{
      kkkk<-matrix(sqrt(B[C,1 ]^2), nrow = n,ncol=n, byrow = TRUE)
      
    }
    I<-matrix(0,nrow=n,ncol=n)
    for (i in 1:ncol(kkkk)){
      I[,i]<-order(kkkk[,i])
    }
    
    d <- y-f
    e<-matrix(0,nrow=n,ncol=n)
    for (i in 1:n){
      e[,i] = d[I[,i]]
    }
    g<-matrix(0,nrow=n,ncol=n)
    for (i in 1:n){
      g[,i] = cumsum(e[,i])
    }
    
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
  message1 <- "The validity test was successfully completed."
  message2 <- "H0: The model is considerd to be valid."
  message3 <- "H1: The model is not considered to be valid."
  
  cat(message1, "\n")
  cat(message2, "\n")
  cat(message3, "\n")
  cat("t-Value:", t, "\n")
  cat("p-Value:", p, "\n")
  
  return(list(t_value = t, p_value = p))
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
  return(model$fitted.values)
}

Example1 <- function(c,n,tau) {
  #' The first example of Frahm (2023)
  #'
  #' This function replicates the first example of Frahm (2023)
  #' @param c A constant
  #' @param n The number of observations
  #' @param tau (optional, default is tau = 1) tau is the variance of the error ϵ in the regression Y = g(X) + ϵ
  #' @return The function returns the p-value of the validity test.
  #' @examples 
  #' Example1(0.25,100)
  #' @author Florian Schuetze 26.08.2024
  #' @references 
  #'Frahm, G., 2023, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'@export
  
  if (missing(tau)) tau <- 1
  x <- rnorm(n, mean=0, sd=1)
  y <- rnorm(n, mean=0, sd=tau)
  y[which(x>0)]<- y[which(x>0)]+c
  y[which(x<=0)]<- y[which(x<=0)]-c
  p<-validity(x,y)
  return(p)
}

Example2 <- function(c,n,tau) {
  #' The second example of Frahm (2023)
  #'
  #' This function replicates the second example of Frahm (2023)
  #' @param c A constant
  #' @param n The number of observations
  #' @param tau (optional, default is tau = 1) tau is the variance of the error ϵ in the regression Y = g(X) + ϵ

  #' @return The function returns the p-value of the validity test.
  #' @examples 
  #' Example2(0.25,100)
  #' @author Florian Schuetze 26.08.2024
  #' @references 
  #'Frahm, G., 2023, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'@export
  
if (missing(tau)) tau <- 1
  x <- rnorm(n, mean=0, sd=1)
  y<-(-1)+x+c*((x^2)-1)+rnorm(n, mean=0, sd=tau)
  p<-validity(x,y)
  return(p)
}

Example3 <- function(c,n,tau) {
  #' The third example of Frahm (2023)
  #'
  #' This function replicates the third example of Frahm (2023)
  #' @param c A constant
  #' @param n The number of observations
  #' @param tau (optional, default is tau = 1) tau is the variance of the error ϵ in the regression Y = g(X) + ϵ
  #' @return The function returns the p-value of the validity test.
  #' @examples 
  #' Example3(0.25,100)
  #' @author Florian Schuetze 26.08.2024
  #' @references 
  #'Frahm, G., 2023, A Test for the Validity of Regression Models. Available on SSRN: https://ssrn.com/abstract=4610329
  #'@export
  
  if (missing(tau)) tau <- 1
  require(MASS)
  sigma<-rbind(c(1,0.5), c(0.5,1))
  mu<-c(0, 0) 
  LK<-as.matrix(mvrnorm(n=n, mu=mu, Sigma=sigma))
  y<-0.25*LK[,1]+0.75*LK[,2]+c*LK[,1]*LK[,2]+rnorm(n, mean=0, sd=tau)
  p<-validity(LK,y)
  return(p)
}