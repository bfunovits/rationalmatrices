# Power series parameters ######################################################

#' Power series Parameters
#' 
#' This function returns the coefficients of the power series expansion (around \eqn{z=0}) of  
#' a rational matrix 
#' \deqn{c(z) = c_0 + c_1 z + c_2 z^2 + \cdots}{c(z) = c[0] + c[1] z + c[2] z^2 + \dots} 
#' Note that for a left matrix fraction description \eqn{c(z) = a(z)^{-1} b(z)} the 
#' matrix \eqn{a(0)}  must be an invertible matrix. 
#' If the rational functions corresponds to a VARMA or statespace model, 
#' then this sequence of coefficients is also 
#' called \emph{impulse response function}. 
#' 
#' A \code{pseries} object is simply a 3-dimensional (numeric or complex) array 
#' of dimension c(m,n,lag.max+1) with a class attribute c('pseries', 'ratm'). 
#' The dimensions \eqn{(m,n)} may also be zero. 
#' 
#'
#' @param obj (rational) matrix object, i.e. a \code{\link{polm}}, \code{\link{lmfd}},  \code{\link{rmfd}}, \code{\link{stsp}} object 
#'            or an object which may be coerced to a polynomial matrix with \code{polm(obj)}. 
#'            The default \code{S3} method first coerces the input argument \code{obj} to a \code{polm} object. 
#'            If this fails an error is thrown. 
#' @param lag.max (integer) maximum lag.
#' @param ... not used.
#'
#' @return object of type \code{pseries}
#' 
#' @rdname pseries
#' @export
pseries = function(obj, lag.max, ...){
  UseMethod("pseries", obj)
}

#' @rdname pseries
#' @export
pseries.default = function(obj, lag.max = 5, ...) {
  # try to coerce to polm
  obj = try(polm(obj), silent = TRUE)
  if (inherits(obj, 'try-error')) stop('could not coerce argument "obj" to "polm" object')
  
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }
  
  obj = unclass(obj)
  m = dim(obj)[1]
  n = dim(obj)[2]
  p = dim(obj)[3] - 1
  
  ir = array(0, dim = unname(c(m,n,lag.max + 1)))
  
  if (p > -1) ir[,,1:min(lag.max+1, p+1)] = obj[,, 1:min(lag.max+1,p+1)]
  
  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}

#' @rdname pseries
#' @export
pseries.polm = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }
  
  obj = unclass(obj)
  m = dim(obj)[1]
  n = dim(obj)[2]
  p = dim(obj)[3] - 1
  
  ir = array(0, dim = unname(c(m,n,lag.max + 1)))
  
  if (p > -1) ir[,,1:min(lag.max+1, p+1)] = obj[,, 1:min(lag.max+1,p+1)]
  
  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}

#' @rdname pseries
#' @export
pseries.lmfd = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }
  
  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  ir = array(0, dim = unname(c(m,n,lag.max +1)))
  
  if ((m*n*(q+1)) == 0) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }
  
  ab = unclass(obj)
  a = ab[,1:(m*(p+1)), drop = FALSE]
  b = ab[,(m*(p+1)+1):(m*(p+1)+n*(q+1)), drop = FALSE]
  
  # check a(z)
  if (p < 0) stop('left factor "a(z)" is not a valid polynomial matrix')
  
  a0 = a[,1:m, drop = FALSE]
  junk = try(solve(a0)) 
  if (inherits(junk, 'try-error')) stop('left factor "a(0)" is not invertible')
  
  # b[,,i] -> a0^{-1} b[,,i]
  b = solve(a0, b)
  dim(b) = c(m, n, q+1)
  
  ir[ , , 1:min(lag.max+1, q+1)] = b[ , , 1:min(lag.max+1, q+1)]
  
  if ((p == 0) || (lag.max == 0)) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }
  
  # a[,,i] -> a0^{-1} a[,,i]
  a = solve(a0, a[,(m+1):ncol(a), drop = FALSE])
  dim(a) = c(m, m, p)
  
  for (lag in (1:lag.max)) {
    for (i in (1:min(p,lag))) {
      ir[,,lag+1] = matrix(ir[,,lag+1], nrow = m, ncol = n) - 
        matrix(a[,,i], nrow = m, ncol = m) %*% matrix(ir[,,lag+1-i], nrow = m, ncol = n)
    }
  }  
  
  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}


#' @rdname pseries
#' @examples 
#' 
#' (obj = rmfd(c = polm(array(c(c(diag(2)), c(0.5, 0.25, 0.125, 0.5)), dim = c(2,2,2))), d = NULL))
#' (k = pseries(obj, lag.max = 4))
#' pseries2rmfd(k)
#' 
#' (obj = lmfd(a = polm(array(c(c(diag(2)), c(0.5, 0.25, 0.125, 0.5)), dim = c(2,2,2)))))
#' (k = pseries(obj, lag.max = 4))
#' pseries2lmfd(k)
#'
#' @export
pseries.rmfd = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }
  
  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  ir = array(0, dim = unname(c(m,n,lag.max +1)))
  
  if ((m*n*(q+1)) == 0) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }
  
  # check c(z)
  if (p < 0) stop('right factor "c(z)" is not a valid polynomial matrix')
  
  # Separate c(z) and d(z)
  cd = unclass(obj) 
  c = cd[1:(n*(p+1)), , drop = FALSE]
  d = cd[(n*(p+1)+1):(n*(p+1)+m*(q+1)), , drop = FALSE]
  
  c0 = matrix(c[1:n, , drop = FALSE], nrow = n, ncol = n)
  c0inv = tryCatch(solve(c0),
                   error = function(cnd) stop(' "c(0)" is not invertible'))
  
  c = c %*% c0inv
  d = d %*% c0inv
  
  # Convert d(z) to an array 
  d = t(d)
  dim(d) = c(n, m, q+1)
  d = aperm(d, perm = c(2,1,3))
  
  # Initialize impulse response with d(z)
  ir[ , , 1:min(lag.max+1, q+1)] = d[ , , 1:min(lag.max+1, q+1)]
  
  if ((p == 0) || (lag.max == 0)) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }
  
  # Convert c(z) an array and discard zero-lag coefficient (equal to I_q)
  c = c[(n+1):nrow(c),, drop = FALSE]
  c = t(c)
  dim(c) = c(n, n, p)
  c = aperm(c, perm = c(2,1,3))
  
  # Calculate impulse response
  for (lag in (1:lag.max)) {
    for (i in (1:min(p,lag))) {
      ir[,,lag+1] = matrix(ir[,,lag+1], nrow = m, ncol = n) - matrix(ir[,,lag+1-i], nrow = m, ncol = n) %*% matrix(c[,,i], nrow = n, ncol = n)
    }
  }  
  
  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}



#' @rdname pseries
#' @export
pseries.stsp = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }
  
  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  ir = array(0, dim = unname(c(m,n,lag.max +1)))
  
  if ((m*n) == 0) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }
  
  ABCD = unclass(obj)
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]
  
  ir = array(0, dim = unname(c(m,n,lag.max + 1)))
  ir[,,1] = D
  
  if ((s == 0) || (lag.max == 0)) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }
  
  for (lag in (1:lag.max)) {
    ir[,,lag+1] = C %*% B
    B = A %*% B
  }  
  
  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}

#' @rdname pseries
#' @export
pseries.lpolm = function(obj, ...) {
  min_deg = attr(obj, which = "min_deg")
  if (min_deg >= 0){
    obj = lpolm(obj) # returns a polm object
    return(pseries(obj))
  } else {
    stop("A lpolm object cannot be coerced to a pseries object.")
  }
}
  
