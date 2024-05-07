validity <- function(X, y, r) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  N <- 1e3
  n <- nrow(X)
  m <- ncol(X)

  if (r == 0) {
    r <- max(floor(log(n)), 2)
  }

  D <- cbind(1, X)
  p <- solve(t(D) %*% D) %*% t(D) %*% y
  alpha <- p[1]
  beta <- p[-1]
  f <- alpha + X %*% beta
  Z <- cbind(X, y, f, y - f)
  Z <- Z[order(Z[, m + 2]), ]

  T2 <- numeric(r)
  for (j in 1:r) {
    l <- ceiling((j - 1) * n / r) + 1
    u <- ceiling(j * n / r)
    e <- Z[l:u, m + 3]
    k <- length(e)
    T2[j] <- sum(e)^2 / k
  }
  T <- sum(T2)

  theta <- valtheta(X, y, alpha, beta, r)
  X_sum <- rowSums(theta^2)

  if (T > 0) {
    P <- 1 - sum(X_sum <= T) / N
  } else {
    P <- 1
  }

  return(list(P = P, T = T, alpha = alpha, beta = beta))
}

valtheta <- function(X, y, alpha, beta, r) {
  N <- 1e3
  n <- nrow(X)
  m <- ncol(X)
  f <- alpha + X %*% beta
  E <- matrix(sample(y - f, n * N, replace = TRUE), nrow = n, ncol = N)
  Y <- alpha + X %*% matrix(rep(beta, N), ncol = N) + E
  D <- cbind(1, X)
  P <- solve(t(D) %*% D) %*% t(D) %*% Y

  theta <- matrix(NA, nrow = N, ncol = r)
  for (i in 1:N) {
    a <- P[1, i]
    b <- P[-1, i]
    F <- a + X %*% b
    Z <- cbind(X, Y[, i], F, Y[, i] - F)
    Z <- Z[order(Z[, m + 2]), ]
    for (j in 1:r) {
      l <- ceiling((j - 1) * n / r) + 1
      u <- ceiling(j * n / r)
      e <- Z[l:u, m + 3]
      k <- length(e)
      theta[i, j] <- sqrt(k) * mean(e)
    }
  }

  return(theta)
}


constant <- function(c, n) {
  x <- rnorm(n, mean=0, sd=1)
  y <- rnorm(n, mean=0, sd=1)
  y[which(x>0)]<- y[which(x>0)]+c
  y[which(x<=0)]<- y[which(x<=0)]-c
  z<-cbind(x,y)
  return(z)
}

quadratic <- function(c, n) {
  x <- rnorm(n, mean=0, sd=1)
  y<-rnorm(n, mean=0, sd=1)-1+x+c*((x^2)-1)
  z<-cbind(x,y)
  return(z)
}

cobbdouglas <- function(c, n) {
  require(MASS)
  sigma<-rbind(c(1,0.5), c(0.5,1))
  mu<-c(0, 0)
  LK<-as.matrix(mvrnorm(n=n, mu=mu, Sigma=sigma))
  y<- rnorm(n, mean=0, sd=1)+0.25*LK[,1]+0.75*LK[,2]+c*LK[,1]*LK[,2]
  d <- data.frame(Y=y,X1=LK[,1],X2=LK[,2])
  return(d)
}

stute1997 <- function(a, sigma) {
  x <- runif(400, min = 0, max = 1)
  y<-  5*x+a*x^2+rnorm(400, mean=0, sd=sqrt(sigma))
  d <- data.frame(Y=y,X1=x)
  return(d)
}
