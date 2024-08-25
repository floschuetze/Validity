validity <- function(y, X, reg, v = 0, s = 0,l=0,r=0) {
  
  #  Test for the validity of a regression model.
  # ---
  # Input:
  # y (n x 1) - dependent variable
  # X (n x m) - regressors
  # reg (@)   - regression function delifering fitted values
  # v (1 x 1) - version of the test
  # v == 0: Frahm-Sch?tze
  # v == 1: Stute
  # s (1 x 1) - test statistic (only relevant if v == 1)
  # s == 0: Cram?r-von-Mises
  # s == 1: Kolmogorov-Smirnov
  # The choice of s does not matter for v = 0.
  # l (1 x 1) - weight on the left of [0,1]
  # r (1 x 1) - weight on the right of [0,1]
  # The choice of l and r does not matter for v = 1.
  # ---
  # Output:
  # p (1 x 1) - p-value
  # t (1 x 1) - test statistic
  #Florian Sch?tze 08/2024 based on Gabriel Frahm 02.06.2024
  
  if (missing(v)) v <- 0
  if (missing(s)) s <- 0
  if (missing(l)) l <- 0
  if (missing(r)) r <- 0
  if (missing(reg)) reg <-function(y, X) {
    d<-data.frame(X,y1=y)
    model <- lm(y1 ~ ., data = d)
    return(model$fitted.values)
  }
  
  N <- 1000
  X<-as.matrix(X)
  n <- nrow(X)
  m <- ncol(X)
  
  f <- reg(y, X)
  Data <- (y - f) - mean(y - f)
  
  E <- matrix(sample(Data, n * N, replace = TRUE), nrow = n, ncol = N)
  
  Y <- f + E
  #Y <- readMat("Y.mat")$Y
  
  F <- matrix(0, nrow = n, ncol = N)
  for (i in 1:N) {
    F[, i] <- reg(Y[, i], X)
  }
  F <- apply(Y, 2, function(y) reg(y, X))
  # install.packages("R.matlab")
  # library(R.matlab)
  #F <- readMat("F.mat")$F
  Z1 <- matrix(NA, nrow = n, ncol = m * n) #vorversion f?r z (repmat(X,1,n);)
  for (i in 1:n) {
    Z1[, (m * i - 1):(m * i)] <- X
  }
  
  
  
  # Repeat the row 7 times
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
  
  # Replicate the column vector across columns
  k <- matrix(rep(col_vector, each = n), nrow = n, byrow = TRUE)
  
  # Compute the Beta PDF values
  beta_pdf_values <- dbeta(((1:n) / n), shape1 = r + 1, shape2 = l + 1)
  
  # Replicate the vector horizontally
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
    
    # Compute the mean along the second dimension (columns)
    t <- mean(apply((w * g * g / k), 2, mean))
    
    # Reshape to a 3D array (100x1x1000)
    R <- array((Y - F), dim = c(n, 1, N))
    #####hier#####
    
    #D ist gesampelt, achtung!
    E <- array(R[matrix(I, n^2, 1),,], dim = c(n, n, N))
    
    G <- apply(E, c(2, 3), cumsum)
    
    
    # Compute the mean across the second dimension (rows in MATLAB)
    T <- apply(apply(((W * G * G) / K), c(1, 3), mean),2,mean)
  } else if (v == 1){
    I <- array(rowSums(B <= 0) == m, dim = c(n, n))
    
    # e is a matrix where rows are repeated 'n' times
    e <- matrix(rep(y - f, n), nrow = length(y), ncol = n)
    # Set elements of e to 0 where I is TRUE
    e[I] <- 0
    
    # E is a 3D array created by replicating Y - F across dimensions
    
    
    # Convert the 2D matrix into a 3D array with dimensions (100, 1, 1000)
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
  # Compute the mean along the first dimension (rows)
  
  # Mean along the first dimension (rows)
  # library(pracma)
  # h <- repmat(matrix(1:n, n, 1), 1, n)
  # t <- mean(apply(g * g / h, 2, mean))
  # H <- array(rep(1:n, times = n * N), dim = c(n, n, N))
  # T <- apply(G * G / H, c(2, 3), mean) 
  # T <- apply(T, 2, mean)  
  
  
  if (t <= .Machine$double.eps){
    t = 0
    p = 1
  }else{
    p <- 1 - sum(T <= t) / N}
  #return(list(P-Value= p , t-Value = t))
  message1 <- "The validity test was successfully completed."
  message2 <- "H0: The model is considerd to be valid."
  message3 <- "H1: The model is not considered to be valid."
  
  
  # Ergebnisse anzeigen
  cat(message1, "\n")
  cat(message2, "\n")
  cat(message3, "\n")
  cat("t-Value:", t, "\n")
  cat("p-Value:", p, "\n")
  
  # Optionale Rückgabe der Werte
  return(list(t_value = t, p_value = p))
}

reg<-function(y, X) {
  d<-data.frame(X,y1=y)
  model <- lm(y1 ~ ., data = d)
  return(model$fitted.values)
}

onlyconstant <- function(c,n) {
  x <- rnorm(n, mean=0, sd=1)
  y <- rnorm(n, mean=0, sd=1)
  y[which(x>0)]<- y[which(x>0)]+c
  y[which(x<=0)]<- y[which(x<=0)]-c
  p<-validity(x,y,reg,0,0,0,0)
  return(p)
}

quadratic <- function(c,n) {
  
  x <- rnorm(n, mean=0, sd=1)
  y<-(-1)+x+c*((x^2)-1)+rnorm(n, mean=0, sd=1)
  p<-validity(x,y,reg,0,0,0,0)
  return(p)
}

cobbdouglas <- function(c,n) {
  
  require(MASS)
  sigma<-rbind(c(1,0.5), c(0.5,1))
  mu<-c(0, 0) 
  LK<-as.matrix(mvrnorm(n=n, mu=mu, Sigma=sigma))
  y<-0.25*LK[,1]+0.75*LK[,2]+c*LK[,1]*LK[,2]+rnorm(n, mean=0, sd=1)
  p<-validity(LK,y,reg,0,0,0,0)
  return(p)
}