# transpose for rational matrices ##############################################  

#' Rational Matrix Transpose
#'
#' Compute the transpose of a rational matrix \eqn{x(z)}. 
#' 
#' @param x rational matrix object, i.e. a \code{\link{polm}}, \code{\link{lpolm}}, 
#'          \code{\link{lmfd}}, \code{\link{rmfd}}, \code{\link{stsp}}, 
#'          \code{\link{pseries}} or \code{\link{zvalues}} object. 
#'
#' @return A rational matrix object,
#'    which represents the transposed rational matrix \eqn{x'(z)}. 
#'    The output is of the same class as the input \code{x} unless 
#'    \code{x} is an \code{\link{lmfd}} or an \code{\link{rmfd}} object: 
#'    The transposition of 
#'    an \code{lmfd} object is an \code{rmfd} object, and vice versa.

#' 
#' @name transpose
#' @rdname transpose
#' @export
#' 
#' @examples 
#' x = test_polm(dim = c(2,3), degree = 3, random = TRUE)
#' all.equal(pseries(t(x)), t(pseries(x)))
#' 
#' x = test_stsp(dim = c(3,2), s = 1)
#' all.equal(zvalues(t(x)), t(zvalues(x)))
#' 
#' # the transpose of an LMFD object is RMFD
#' (x = test_lmfd(dim = c(3,2), degrees = c(1,1)))
#' t(x)
#' all.equal(x, t(t(x)))
#' 
t.polm = function(x) {
  x = unclass(x)
  x = aperm(x, c(2,1,3))
  return(polm(x))
}

#' @rdname transpose
#' @export
t.lpolm = function(x) {
  x = unclass(x)
  min_deg = attr(x, which = "min_deg")
  x = aperm(x, c(2,1,3))
  return(lpolm(x, min_deg = min_deg))
}


#' @rdname transpose
#' @export
t.lmfd = function(x) {
  rmfd(c = t(x$a), d = t(x$b))
}  

#' @rdname transpose
#' @export
t.rmfd = function(x) {
  lmfd(a = t(x$c), b = t(x$d))
}  

  
#' @rdname transpose 
#' @export
t.stsp = function(x) {
  
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  x = unclass(x)
  A = t(x[iseq(1,s), iseq(1,s), drop = FALSE])
  C = t(x[iseq(1,s), iseq(s+1,s+n), drop = FALSE])
  B = t(x[iseq(s+1,s+m), iseq(1,s), drop = FALSE])
  D = t(x[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE])
  
  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}

#' @rdname transpose
#' @export
t.pseries = function(x) {
  
  x = aperm(unclass(x), c(2,1,3))
  x = structure(x, class = c('pseries','ratm'))
  return(x)
}

#' @rdname transpose
#' @export
t.zvalues = function(x) {
  
  z = attr(x, 'z')
  x = aperm(unclass(x), c(2,1,3))
  x = structure(x, z = z, class = c('zvalues','ratm'))
  return(x)
}

# Hermitian Transpose ####

#' Hermitean Transpose 
#'
#' The \emph{Hermitean transpose} of a rational matrix \eqn{x(z)} is defined as 
#' \deqn{x^{*}(z)=\overline{x(\bar{z}^{-1})}'.}{x^*(z)=Conj(x(Conj(z)^{-1}))'.} 
#' This means e.g. for a polynomial with real coefficients 
#' \eqn{x(z)=a_0 + a_1 z +\cdots + a_p z^p}{x(z)=a[0] + a[1] z +\dots + a[p] z^p} 
#' that the coefficient matrices are transposed and that \eqn{z} is replaced 
#' by \eqn{z^{-1}}: 
#' \eqn{x^*(z)=a_0' + a_1' z^{-1} +\cdots + a_p' z^{-p}}{
#'      x^*(z)=a[0]' + a[1]' z^{-1} +\dots + a[p]' z^{-p}}, i.e. 
#' the result is (in general) a Laurent polynomial. 
#' 
#' The Hermitean transpose is only implemented for 
#' polynomial matrices (\code{\link{polm}} objects), 
#' Laurent polynomials \code{\link{lpolm}} objects), 
#' rational matrices in statespace form (\code{\link{stsp}} objects) and for 
#' the frequency response of rational matrices (\code{\link{zvalues}} objects).  
#' 
#' The case of rational matrices in LMFD or RMFD form has not (yet) been implemented. 
#' It is not possible to directly construct the Hermitean transpose from 
#' given power series coefficients. Therefore \code{Ht()} does 
#' not support \code{\link{pseries}} objects.
#' 
#' Finally note that the Hermitean transpose of a rational matrix 
#' \eqn{x(z)=D+zB(I-Az)^{-1}B} in statespace form 
#' (in general) has a pole at zero, if the state transition matrix \eqn{A} 
#' is singular. Since rational matrices with a pole at \eqn{z=0} have \bold{no} 
#' statespace realization, the procedure \code{Ht()} throws an error in this case.
#' 
#' @param x rational matrix object of class \code{\link{polm}}, \code{\link{lpolm}}, 
#'          \code{\link{stsp}} or \code{\link{zvalues}}. 
#'
#' @return A rational matrix object which represents the Hermitean transpose \eqn{x^*(z)}. 
#'         The Hermitean transpose of a polynomial matrix (in general) is a Laurent polynomial, 
#'         therefore the output is of class \code{lpolm} if \code{x} is a \code{polm} object. 
#'         In all other cases the output is of the same class as the input \code{x}. 
#' 
#' @rdname Ht
#' @aliases Hermitean-transpose
#' @export
#' 
#' @examples 
#' # rational matrix in statespace form 
#' x = test_stsp(dim = c(4,4), s = 2)
#' Ht(x)
#' all.equal(zvalues(Ht(x)), Ht(zvalues(x)))
#' 
#' # Note (x(z)  x^*(z)) is Hermitean for |z| = 1
#' xx = zvalues(x) %r% Ht(zvalues(x))
#' # print(xx, digits = 2)
#' apply(unclass(xx), MARGIN = 3, FUN = isSymmetric)
#' 
#' # polynomial matrix 
#' x = test_polm(dim = c(2,3), degree = 1)
#' Ht(x)
#' 
#' \dontrun{
#' Ht(test_lmfd(dim = c(2,2), degrees = c(3,3)))
#' Ht(test_rmfd(dim = c(2,2), degrees = c(3,3)))
#' Ht(pseries(test_lmfd(dim = c(2,2), degrees = c(3,3))))
#' }
#' 
Ht = function(x) {
  UseMethod("Ht", x)
}

#' @rdname Ht
#' @export
Ht.polm = function(x) {
  x = unclass(x)
  max_deg = dim(x)[3] - 1
  x = aperm(Conj(x), c(2, 1, 3))
  if (dim(x)[3] > 1) x = x[ , , (dim(x)[3]):1, drop = FALSE]
  return(structure(x, class = c("lpolm", 'ratm'), min_deg = -max_deg)) 
}


#' @rdname Ht
#' @export
Ht.lpolm = function(x) {
  x = unclass(x)
  min_deg = attr(x, which = "min_deg")
  max_deg = dim(x)[3] - 1 + min_deg
  x = aperm(Conj(x), c(2, 1, 3))
  if (dim(x)[3] > 1) x = x[ , , (dim(x)[3]):1, drop = FALSE]
  return(structure(x, class = c("lpolm", 'ratm'), min_deg = -max_deg)) 
}

# #' @rdname Ht
# #' @export
# Ht.lmfd = function(x) {
#   stop('Hermitean transpose of "lmfd" objects is not implemented')
# }
# 
# #' @rdname Ht
# #' @export
# Ht.rmfd = function(x) {
#   stop('Hermitean transpose of "rmfd" objects is not implemented')
# }
# 



#' @rdname Ht
#' @export
Ht.stsp = function(x) {
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  if ((m*n) == 0) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), 
             B = matrix(0, nrow = 0, ncol = m), 
             C = matrix(0, nrow = n, ncol = 0), 
             D = matrix(0, nrow = n, ncol = m))
    return(x)
  }
  
  ABCD = Conj(unclass(x))
  # if (is.complex(ABCD)) {
  #   stop('Hermitean transpose is only implemented for "real" state space models')
  # }
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]
  
  if (s == 0) {
    x = stsp(A, B = t(C), C = t(B), D = t(D))
    return(x)
  }
  
  A = try(solve(A), silent = TRUE)
  if (inherits(A, 'try-error')) {
    stop('The state transition matrix "A" is singular.')
  }
  
  junk = -t( A %*% B)         # C -> - B' A^{-T}
  D = t(D) + junk %*% t(C)    # D -> D' - B' A^{-T} C'
  B = t( C %*% A )            # B -> A^{-T} C'
  A = t(A)                    # A -> A^{-T}
  C = junk
  
  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}


#' @rdname Ht
#' @export
Ht.zvalues = function(x) {
  z = 1/Conj(attr(x,'z'))
  
  x = Conj(aperm(unclass(x), c(2,1,3)))
  x = structure(x, z = z, class = c('zvalues','ratm'))
  
  return(x)
}

