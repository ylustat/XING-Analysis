# Adapted from "bindata" package.
library(stats)
library(mvtnorm)
rmvbin <- function (n, margprob, commonprob = diag(margprob), bincorr = diag(length(margprob)), 
                    sigma = diag(length(margprob)), colnames = NULL, simulvals = NULL) {
  if (missing(sigma)) {
    if (!missing(commonprob)) {
      if (missing(margprob)) 
        margprob <- diag(commonprob)
      sigma <- commonprob2sigma(commonprob, simulvals)
    } else if (!missing(bincorr)) {
      commonprob <- bincorr2commonprob(margprob, bincorr)
      sigma <- commonprob2sigma(commonprob, simulvals)
    }
  }else if (any(eigen(sigma)$values < 0)) 
    stop("Sigma is not positive definite.")
  retval <- rmvnorm(n, qnorm(margprob), as.matrix(sigma))
  retval <- ra2ba(retval)
  dimnames(retval) <- list(NULL, colnames)
  retval
}
bincorr2commonprob <- function (margprob, bincorr) 
{
  retval <- 0 * bincorr
  for (k in 1:ncol(retval)) {
    for (l in 1:ncol(retval)) {
      retval[k, l] <- bincorr[k, l] * sqrt(margprob[k] * 
                                             (1 - margprob[k]) * margprob[l] * (1 - margprob[l])) + 
        margprob[k] * margprob[l]
    }
  }
  retval
}

ra2ba <- function (x) 
{
  retval <- as.numeric(x > 0)
  dim(retval) <- dim(x)
  retval
}

commonprob2sigma <- function (commonprob, simulvals = NULL) {
  if (is.null(simulvals)) {
    simulvals <- SimulVals
  }
  margprob <- diag(commonprob)
  if (!(check <- check.commonprob(commonprob))) {
    cat(attr(check, "message"), sep = "\n")
    stop("Matrix commonprob not admissible.")
  }
  sigma <- diag(nrow(commonprob))
  for (m in 1:(ncol(commonprob) - 1)) {
    for (n in (m + 1):nrow(commonprob)) {
      x <- cbind(margprob[m], margprob[n], as.numeric(dimnames(simulvals)[[3]]))
      y <- interpolate(x, simulvals)
      f <- approxfun(y, x[, 3])
      sigma[m, n] <- sigma[n, m] <- f(commonprob[m, n])
    }
  }
  if (any(is.na(sigma))) 
    stop("Extrapolation occurred ... margprob and commonprob not compatible?")
  if (any(eigen(sigma)$values < 0)) {
    cat("Warning: Resulting covariance matrix is not positive definite.\n")
    cat("         Smallest eigenvalue equals", min(eigen(sigma)$values), 
        ".\n")
    cat("         Please check whether the results are still useful.\n")
  }
  sigma
}


interpolate <- function (x, a, adims = lapply(dimnames(a), as.numeric), method = "linear") 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  if (!is.array(a)) 
    stop("a is not an array")
  ad <- length(dim(a))
  method <- pmatch(method, c("linear", "constant"))
  if (is.na(method)) 
    stop("invalid interpolation method")
  if (any(unlist(lapply(adims, diff)) < 0)) 
    stop("dimensions of a not ordered")
  retval <- rep(0, nrow(x))
  bincombi <- bincombinations(ad)
  convexcoeff <- function(x, y) {
    ok <- y > 0
    x[ok] <- y[ok] - x[ok]
    x
  }
  for (n in 1:nrow(x)) {
    leftidx <- rep(0, ad)
    xabstand <- rep(0, ad)
    aabstand <- rep(0, ad)
    for (k in 1:ad) {
      if (x[n, k] < min(adims[[k]]) || x[n, k] > max(adims[[k]])) 
        stop("No extrapolation allowed")
      else {
        leftidx[k] <- max(seq(adims[[k]])[adims[[k]] <= 
                                            x[n, k]])
        if (leftidx[k] == length(adims[[k]])) 
          leftidx[k] <- leftidx[k] - 1
        xabstand[k] <- x[n, k] - adims[[k]][leftidx[k]]
        aabstand[k] <- adims[[k]][leftidx[k] + 1] - 
          adims[[k]][leftidx[k]]
      }
    }
    coefs <- list()
    if (method == 1) {
      for (k in 1:(2^ad)) {
        retval[n] <- retval[n] + element(a, leftidx + 
                                           bincombi[k, ]) * prod((aabstand - convexcoeff(xabstand, 
                                                                                         aabstand * bincombi[k, ]))/aabstand)
      }
    }
    else if (method == 2) {
      retval[n] <- element(a, leftidx)
    }
  }
  names(retval) <- rownames(x)
  retval
}
element <- function (x, i) 
{
  if (!is.array(x)) 
    stop("x is not an array")
  ni <- length(i)
  dx <- dim(x)
  if (length(i) != length(dx)) 
    stop("Wrong number of subscripts")
  if (ni == 1) {
    return(x[i])
  }
  else {
    m1 <- c(i[1], i[2:ni] - 1)
    m2 <- c(1, cumprod(dx)[1:(ni - 1)])
    return(x[sum(m1 * m2)])
  }
}
check.commonprob <- function (commonprob) 
{
  retval <- TRUE
  message <- character(0)
  nm <- 0
  if ((any(commonprob < 0)) || (any(commonprob > 1))) {
    retval <- FALSE
    message[nm <- nm + 1] <- "Not all probabilities are between 0 and 1."
  }
  n <- dim(commonprob)[1]
  if (n != dim(commonprob)[2]) {
    retval <- FALSE
    message[nm <- nm + 1] <- "Matrix of common probabilities is not quadratic."
  }
  for (i in 1:(n - 1)) {
    for (j in 2:n) {
      ul <- min(commonprob[i, i], commonprob[j, j])
      ll <- max(commonprob[i, i] + commonprob[j, j] - 
                  1, 0)
      if ((commonprob[i, j] > ul) || (commonprob[i, j] < 
                                      ll)) {
        retval <- FALSE
        message[nm <- nm + 1] <- paste("Error in Element (", 
                                       i, ",", j, "): Admissible values are in [", 
                                       ll, ",", ul, "].")
      }
    }
  }
  if (n > 2) 
    for (i in 1:(n - 2)) for (j in (i + 1):(n - 1)) for (k in (j + 
                                                               1):n) {
      l <- commonprob[i, i] + commonprob[j, j] + commonprob[k, 
                                                            k] - 1
      if (commonprob[i, j] + commonprob[i, k] + commonprob[j, 
                                                           k] < l) {
        retval <- FALSE
        message[nm <- nm + 1] <- paste("The sum of the common probabilities of", 
                                       i, ",", j, ",", k, "must be at least", l, 
                                       ".")
      }
    }
  attr(retval, "message") <- message
  retval
}

bincombinations <- function (p) 
{
  retval <- matrix(0, nrow = 2^p, ncol = p)
  for (n in 1:p) {
    retval[, n] <- rep(c(rep(0, (2^p/2^n)), rep(1, (2^p/2^n))), 
                       length.out = 2^p)
  }
  retval
}

