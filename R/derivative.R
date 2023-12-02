# derivative.R #######################################
# compute derivatives of rational matrices



#' Derivative of a rational Matrix
#'
#' Computes the derivative of a rational matrix \eqn{k(z)} (with repect to the complex variable \eqn{z}). 
#' Note that computing the derivative for an impulse response object decreases the number of lags by one!
#' 
#' @param obj an object of class \code{\link{polm}},  \code{\link{stsp}} or  \code{\link{pseries}}.
#' @param ... not used
#'
#' @return an object of the same class as the argument \code{obj}
#' @export
#'
#' @examples
#' # create random (3 by 2) polynomial matrix with degree 2
#' K = test_polm(dim = c(3,2), degree = 2)
#' derivative(K)
#' derivative(derivative(K))
#' derivative(derivative(derivative(K)))
#' 
#' # note: computing the derivative of the impulse response 
#' # decreases "lag.max" by one!
#' all.equal(pseries(derivative(K)), derivative(pseries(K, lag.max = 6)))
#' 
#' # create statespace realization of a random (3 by 2) rational matrix
#' # with statespace dimension s = 4
#' K = test_stsp(dim = c(2,2), s = 4, bpoles = 1)
#' all.equal(pseries(derivative(K)), derivative(pseries(K, lag.max = 6)))
#' 
#' \dontrun{
#' # 'lmfd' objects and 'zvalues' objects are not supported
#' derivative(test_lmfd(dim = c(3,3), degrees = c(1,1)))
#' derivative(zvalues(K))
#' }
derivative = function(obj, ...){
  UseMethod("derivative", obj)
}

#' @rdname derivative
#' @export
derivative.lpolm = function(obj, ...) {
  stop('computation of the derivative is not implemented for "lpolm" objects')
}


#' @rdname derivative
#' @export
derivative.polm = function(obj, ...) {
  obj = unclass(obj)
  dim = dim(obj)
  m = dim[1]
  n = dim[2]
  p = dim[3] - 1
  
  if (p <= 0) {
    obj = array(0, dim = c(m,n,0))
    class(obj) = c('polm','ratm')
    return(obj)
  }
  
  obj = obj[,,-1,drop = FALSE]
  k = array(1:p, dim = c(p,m,n))
  k = aperm(k, c(2,3,1))
  obj = obj*k
  class(obj) = c('polm','ratm')
  return(obj)
}


#' @rdname derivative
#' @export
derivative.lmfd = function(obj, ...) {
  stop('computation of the derivative is not implemented for "lmfd" objects')
}

#' @rdname derivative
#' @export
derivative.rmfd = function(obj, ...) {
  stop('computation of the derivative is not implemented for "rmfd" objects')
}

#' @rdname derivative
#' @export
derivative.stsp = function(obj, ...) {
  dim = dim(obj)
  m = dim[1]
  n = dim[2]
  s = dim[3]
  if (min(dim) <= 0) {
    obj = stsp(matrix(0, nrow = 0, ncol = 0), matrix(0, nrow = 0, ncol = n),
               matrix(0, nrow = m, ncol = 0), matrix(0, nrow = m, ncol = n))
    return(obj)    
  }
  A = obj$A
  B = obj$B
  C = obj$C
  
  obj = stsp(rbind( cbind(A, diag(s)), cbind(matrix(0, nrow = s, ncol = s), A) ),
             rbind( B, A %*% B), cbind(C %*%A, C), C %*% B)
  return(obj)
}


#' @rdname derivative
#' @export
derivative.pseries = function(obj, ...) {
  obj = unclass(obj)
  dim = dim(obj)
  m = dim[1]
  n = dim[2]
  p = dim[3] - 1
  
  if (p <= 0) {
    obj = array(0, dim = c(m,n,1))
    class(obj) = c('pseries','ratm')
    return(obj)
  }
  
  obj = obj[,,-1,drop = FALSE]
  k = array(1:p, dim = c(p,m,n))
  k = aperm(k, c(2,3,1))
  obj = obj*k
  class(obj) = c('pseries','ratm')
  return(obj)
}


#' @rdname derivative
#' @export
derivative.zvalues = function(obj, ...) {
  stop('computation of the derivative is not implemented for "zvalues" objects')
}
