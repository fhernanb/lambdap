#' Normality and White Noise Testing
#'
#' @param x is a vector.
#'
#' @return The p-value for the test.
#'
#' @details This function is a modification of the normwhn.test::normality.test1 function
#'
#' @export
#'
#' @example examples/examples_dh_test.R
#'
dh.test <- function (x){
  x <- as.data.frame(x)
  n <- as.double(nrow(x))
  nvars <- as.double(ncol(x))
  m1 <- m2 <- m3 <- m4 <- numeric(nvars)
  sk <- k <- numeric(nvars)
  for (p in 1:nvars){
    m1[p] <- mean(x[, p])
    m2[p] <- (n^(-1)) * sum((x[, p] - m1[p])^2)
    m3[p] <- (n^(-1)) * sum((x[, p] - m1[p])^3)
    m4[p] <- (n^(-1)) * sum((x[, p] - m1[p])^4)
    sk[p] <- m3[p]/(m2[p]^(3/2))
    k[p] <- m4[p]/m2[p]^2
  }
  mean.vector <- t(as.matrix(cov.wt(x)$center))
  covariance.matrix <- cov.wt(x)$cov
  diagcov <- diag(covariance.matrix)
  idiagcov <- 1/sqrt(diagcov)
  v.matrix <- matrix(data = 0, nrow = nvars, ncol = nvars)
  for (i in 1:nvars) v.matrix[i, i] <- idiagcov[i]
  correlation.matrix <- cor(x)
  lambda.matrix <- diag(eigen(correlation.matrix)$values)
  h.matrix <- eigen(correlation.matrix)$vectors
  xhat.matrix <- x - t(matrix(rep(mean.vector, n), nvars))
  xhatprime.matrix <- t(xhat.matrix)
  rprime.matrix <- h.matrix %*% solve(lambda.matrix)^(1/2) %*%
    t(h.matrix) %*% v.matrix %*% xhatprime.matrix
  rprime <- rprime.matrix
  trprime <- t(rprime)
  v1 <- v2 <- v3 <- v4 <- rtb1 <- b2 <- numeric(nvars)
  for (j in 1:nvars){
    v1[j] <- mean(trprime[, j])
    v2[j] <- (n^(-1)) * sum((trprime[, j] - v1[j])^2)
    v3[j] <- (n^(-1)) * sum((trprime[, j] - v1[j])^3)
    v4[j] <- (n^(-1)) * sum((trprime[, j] - v1[j])^4)
    rtb1[j] <- v3[j]/(v2[j]^(3/2))
    b2[j] <- v4[j]/v2[j]^2
  }
  beta <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3))/((n -
                                                            2) * (n + 5) * (n + 7) * (n + 9))
  w2 <- (-1) + sqrt(2 * (beta - 1))
  delta <- 1/sqrt(log(sqrt(w2)))
  f <- (w2 - 1)/2
  g <- (n + 1) * (n + 3)/(6 * (n - 2))
  h <- sqrt(f * g)
  y <- rtb1 * h
  z1 <- delta * log(y + sqrt(y^2 + 1))
  del <- ((n - 3) * (n + 1) * (n^2 + 15 * n - 4))
  aye <- ((n - 2) * (n + 5) * (n + 7) * (n^2 + 27 * n - 70))/(6 *
                                                                del)
  cee <- ((n - 7) * (n + 5) * (n + 7) * (n^2 + 2 * n - 5))/(6 *
                                                              del)
  alp <- aye + ((rtb1^2) * cee)
  kap <- ((n + 5) * (n + 7) * (n^3 + 37 * n^2 + 11 * n - 313))/(12 *
                                                                  del)
  chi <- (b2 - 1 - rtb1^2) * (2 * kap)
  chi <- abs(chi)
  z2 <- (((chi/(2 * alp))^(1/3)) - 1 + (1/((9 * alp)))) * sqrt(9 *
                                                                 alp)
  pvalsk <- c(numeric(nvars))
  for (p in 1:nvars) pvalsk[p] <- pnorm(z1[p], lower.tail = FALSE)
  for (p in 1:nvars) pvalsk[p] <- 2 * pvalsk[p]
  for (p in 1:nvars) if (pvalsk[p] > 1)
    pvalsk[p] <- 2 - pvalsk[p]
  pskneg <- c(numeric(nvars))
  for (p in 1:nvars) pskneg[p] <- pnorm(z1[p])
  pskpos <- 1 - pskneg
  pvalk <- c(numeric(nvars))
  for (p in 1:nvars) pvalk[p] <- pnorm(z2[p], lower.tail = FALSE)
  for (p in 1:nvars) pvalk[p] <- 2 * pvalk[p]
  for (p in 1:nvars) if (pvalk[p] > 1)
    pvalk[p] <- 2 - pvalk[p]
  pkneg <- c(numeric(nvars))
  for (p in 1:nvars) pkneg[p] <- pnorm(z2[p])
  pkpos <- 1 - pkneg
  z1 <- matrix(z1, nrow = 1)
  z2 <- matrix(z2, nrow = 1)
  Ep <- z1 %*% t(z1) + z2 %*% t(z2)
  dof <- 2 * nvars
  sig.Ep <- 1 - pchisq(Ep, dof)
  as.numeric(sig.Ep)
}
