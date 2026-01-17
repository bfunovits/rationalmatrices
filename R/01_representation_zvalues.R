# ===================================================================
# Functional Area #1: Representation Management (CLAUDE.md)
#
# Purpose: Create and manage frequency response values
#   - zvalues(): Store function values evaluated at z points
#   - Frequency response on unit circle or arbitrary z points
#
# Related Files: 01_representation_classes.R (constructors), 01_representation_conversions.R (conversions to/from zvalues)
# ===================================================================

# frequency response function ######################################################

#' Frequency Response Function
#' 
#' This function evaluates a rational matrix at given (complex) arguments. 
#' If the rational matrix corresponds to a \emph{linear, dynamic filter}, then 
#' this rational function is also called \emph{transfer function} of the filter 
#' and the values of the function on the complex unit circle are called 
#' \emph{frequency response} of the filter.  
#'
#' A \code{zvalues} object is simply a (complex valued) 3-dimensional array of dimension \eqn{(m,n,n.f)} 
#' with attribute \code{z} and class attribute \code{c('zvalues','ratm')}. 
#' The dimensions \eqn{(m,n,n.f)} may also be zero, where \code{n.f} is the length of \code{z}.
#' 
#' \code{z} and \code{f} can be accessed through \code{zvalues_obj$z} and \code{zvalues_obj$z}, respectively.
#' 
#' If you only need the value of the matrix evaluated at one point, then use 
#' \code{\link{zvalue}(obj, z)}.
#' 
#' @param obj (rational) matrix object, i.e. a \code{\link{polm}}, \code{\link{lpolm}}, \code{\link{lmfd}}, \code{\link{rmfd}}, \code{\link{stsp}} object. 
#'            or an object which may be coerced to a polynomial matrix with \code{polm(obj)}.
#'            The default \code{S3} method acts as a constructor for zvalues objects: and array object and \code{z} are supplied, corresponding to the evaluated object at the values of \code{z}. (Numeric or complex) vectors or matrices are coerced arrays. If \code{f} is supplied instead of \code{z}, then \code{z} is considered to be \code{exp(complex(imaginary = (-2*pi)*f))}.
#' @param z (numeric or complex) vector, points at which to evaluate the rational matrix.
#'   Theses values are stored as attribute and can be accessed through \code{zvalues_obj$z}.
#' @param f (numeric) vector of frequencies. If \code{z = NULL} then \code{z} is set to 
#'          \code{z = exp(complex(imaginary = (-2*pi)*f))}. If \code{z} is non \code{NULL}, then 
#'          \code{f} is set to \code{f = -Arg(z)/(2*pi)}.
#'          Theses values can be accessed through \code{zvalues_obj$f}.
#' @param n.f (integer) number of frequencies. If \code{z = f = NULL} then a grid of frequencies \code{f = (0:(n.f-1))/n.f} 
#'        is used (and \code{z} is generated as explained above).
#' @param sort.frequencies boolean, sort frequencies
#' @param ... optional additional parameters
#'
#' @return Object of type \code{zvalues}
#' 
#' @rdname zvalues
#' @export
zvalues = function(obj, z, f, n.f, sort.frequencies, ...){
  UseMethod("zvalues", obj)
}

# Internal function used by zvalues.____ methods
# Intended to check the inputs to zvalues() and bring them in a normalized form
get_z = function(z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE) {
  
  # If there are no values in z, 
  # then either take the frequencies directly,
  # or create them from the "number of frequencies" n.f
  if (is.null(z)) {
    if (is.null(f)) {
      n.f = as.vector(as.integer(n.f))[1]
      stopifnot("Number of frequencies must be a non negative integer" = (n.f >= 0))
      if (n.f == 0) {
        f = double(0)
      } else {
        f = (0:(n.f-1))/n.f
      }
    } 
    
    if (!is.numeric(f)) {
      stop('"f" must be a vector of (real) frequencies')
    } else {
      f = as.vector(f)
      n.f = length(f)
      z = exp(complex(imaginary = (-2*pi)*f))
    }
  } 
  
  # Check non-NULL inputs z
  if (! (is.numeric(z) || is.complex(z))) stop('"z" must be a numeric or complex valued vector')
  z = as.vector(z)
  
  # Bring potential other arguments in line with z values
  n.f = length(z)
  f = -Arg(z)/(2*pi)
  
  # Sort if requested
  if ((sort.frequencies) && (n.f > 0)) {
    i = order(f)
    z = z[i]
    f = f[i]
  }
  
  return(list(z = z, f = f, n.f = n.f))
}

#' @rdname zvalues
#' @export
#' @examples 
#' (vv = zvalues(1:3, z = 4:6))
#' (mm = zvalues(matrix(1:4, 2, 2), z = 1))
#' (aa = zvalues(array(1:8, dim = c(2,2,2)), f = c(0.2, 0.8)))
zvalues.default = function(obj, z = NULL, f = NULL, n.f = NULL, sort.frequencies = FALSE,  ...) {
  
  stopifnot("For the default method (constructing a zvalues object), *obj* needs to be an array, matrix, or vector." = any(class(obj) %in% c("matrix", "array")) || is.vector(obj),
            "Input *obj* must be numeric or complex" = (is.numeric(obj) || is.complex(obj)),
            "For the default method (constructing a zvalues object), *z*, *f*, or *n.f* needs to be supplied. Use of *z* is recommended." = any(!is.null(z), !is.null(f), !is.null(f)),
            "For the default method (constructing a zvalues object), *sort.frequencies* needs to be FALSE" = !sort.frequencies)
  
  if (is.vector(obj)) {
    dim(obj) = c(1,1,length(obj))
  }
  if (is.matrix(obj)) {
    dim(obj) = c(dim(obj),1)
  }
  
  stopifnot("could not coerce input parameter *obj* to a valid 3-D array!" = (is.array(obj)) && (length(dim(obj)) == 3))
  
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f
  
  d = dim(obj)
  m = d[1]
  n = d[2]

  stopifnot("Length of *obj* (i.e. length of vector, third dimension of array) needs to be equal to *n.f* or length of *z*, *f*." = (d[3] == n.f))
  
  # Empty zvalues object
  if (prod(d)*n.f == 0) {
    fr = array(0, dim = unname(c(m, n, n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  fr = structure(obj, z = z, class = c('zvalues','ratm'))
  return(fr)
}


#' @rdname zvalues
#' @export
zvalues.polm = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  
  # Obtain checked and normalized values
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z # values on the unit circle (if z,f NULL)
  f = zz$f # between -0.5 and 0.5
  n.f = zz$n.f
  
  # Integer-valued parameters
  obj = unclass(obj)
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  
  # Check if dimension, degree, or number of frequencies is zero
  if ((m*n*(p+1)*n.f) == 0) {
    fr = array(0, dim = unname(c(m, n, n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  # Copy highest coefficient n.f times in an array as initial value
  fr = array(obj[,, p+1], dim = unname(c(m, n, n.f)))
  
  if (p == 0) {
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
 
  # Array of dim ( out, in, n_points): every element (i,j) of a polynomial will be evaluated at n_f values
  zz = aperm(array(z, dim = c(n.f, m, n)), c(2,3,1))
  
  # Horner-scheme for evaluation
  for (j in (p:1)) {
    fr = fr*zz + array(obj[,,j], dim = unname(c(m,n,n.f)))
  }
  
  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

#' @rdname zvalues
#' @export
zvalues.lpolm = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  
  # Obtain checked and normalized values
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z # values on the unit circle (if z,f NULL)
  f = zz$f # between -0.5 and 0.5
  n.f = zz$n.f
  
  # Integer-valued parameters
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] # Dimensions are retrieved from polm object, NOT from unclassed polm obj = array
  min_deg = d[4]
  
  obj = obj %>% unclass()
  
  # Check if a dimension or number of frequencies is zero (degree not part of this in contrast to ___.polm() function)
  if ((m*n*n.f) == 0) {
    fr = array(0, dim = unname(c(m, n, n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  # Copy highest coefficient n.f times in an array as initial value
  fr = array(obj[,, p+1-min_deg], dim = unname(c(m, n, n.f)))
  
  if (p == 0 && min_deg == 0) {
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  # Array of dim ( out, in, n_points): every element (i,j) of a polynomial will be evaluated at n_f values
  zz = aperm(array(z, dim = c(n.f, m, n)), c(2,3,1))
  
  # Horner-scheme for evaluation
  for (j in ((p-min_deg):1)) {
    fr = fr*zz + array(obj[,,j], dim = unname(c(m,n,n.f)))
  }
  
  # Every polynomial (i,j) - evaluated at n_f points - needs to be premultiplied by z_0^{min_deg} (where z_0 denotes one particular point of evaluation)
  zz_factor_min_deg = aperm(array(z, dim = c(n.f, m, n)), perm = c(2,3,1))^min_deg
  fr = fr * zz_factor_min_deg
  
  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

#' @rdname zvalues
#' @export
#' @examples 
#' (ll = test_lmfd(dim = c(2,2), bpoles = 1, bzeroes = 1))
#' zvalues(ll)
zvalues.lmfd = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f

  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  if ((m<1) || (p<0)) stop('obj is not a valid "rldm" object')
  
  if ((m*n*(q+1)*n.f) == 0) {
    fr = array(0, dim = unname(c(m,n,n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  fr_a = unclass(zvalues(obj$a, z = z))
  fr = unclass(zvalues(obj$b, z = z))
  
  for (i in (1:n.f)) {
    # x = try(solve(matrix(fr_a[,,i], nrow = m, ncol = m), 
    #               matrix(fr[,,i], nrow = m, ncol = n)), silent = TRUE)
    # # print(x)
    # if (!inherits(x, 'try-error')) {
    #   fr[,,i] = x
    # } else {
    #   fr[,,i] = NA
    # }
    fr[,,i] = tryCatch(solve(matrix(fr_a[,,i], nrow = m, ncol = m), 
                             matrix(fr[,,i], nrow = m, ncol = n)), 
                       error = function(e) NA_real_)
  }
  
  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

#' @rdname zvalues
#' @export
zvalues.rmfd = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f
  
  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  if ((m<1) || (p<0)) stop('obj is not a valid "rldm" object')
  
  if ((m*n*(q+1)*n.f) == 0) {
    fr = array(0, dim = unname(c(m,n,n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  # compute first transposed values: t(c(z))^(-1) t(d(z))
  fr_c = aperm(unclass(zvalues(obj$c, z = z)), c(2,1,3))
  fr = aperm(unclass(zvalues(obj$d, z = z)), c(2,1,3))
  
  for (i in (1:n.f)) {
    # fr[,,i] = tryCatch(matrix(fr[,,i], nrow = m, ncol = n) %*% 
    #                      solve(matrix(fr_c[,,i], nrow = n, ncol = n)),
    #                    error = function(e) NA)
    fr[,,i] = tryCatch(solve(matrix(fr_c[,,i], nrow = n, ncol = n),  
                             matrix(fr[,,i], nrow = n, ncol = m)),
                       error = function(e) NA_real_)
  }
  # undo transposition
  fr = aperm(fr, c(2,1,3))
  
  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

#' @rdname zvalues
#' @export
zvalues.stsp = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  if ((m*n*n.f) == 0) {
    fr = array(0, dim = unname(c(m,n,n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  ABCD = unclass(obj)
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]
  
  fr = array(D, dim = unname(c(m,n,n.f)))
  
  if (s == 0) {
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }
  
  Is = diag(s)
  for (i in (1:n.f)) {
    # zz = z[i]^{-1}
    # x = try(C %*% solve(diag(zz, nrow = s, ncol = s) - A, B), silent = TRUE)
    # if (!inherits(x, 'try-error')) {
    #   fr[,,i] = fr[,,i] + x
    # } else {
    #   fr[,,i] = NA
    # }
    fr[,,i] = fr[,,i] + z[i]*tryCatch(C %*% solve(Is - z[i]*A, B), 
                                 error = function(e) NA_real_)
  }
  
  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

# Evaluate at one point only ####

#' Frequency Response Function
#' 
#' This function evaluates a rational matrix at one given (complex) arguments. 
#' It is the non-plural version of \code{\link{zvalues}}.
#' 
#' @param obj (rational) matrix object, i.e. a \code{\link{polm}}, \code{\link{lpolm}}, \code{\link{lmfd}}, \code{\link{rmfd}}, \code{\link{stsp}} object 
#'            or an object which may be coerced to a polynomial matrix with \code{polm(obj)}.
#'            The default \code{S3} method first coerces the input argument \code{obj} to a \code{polm} object. 
#'            If this fails an error is thrown. 
#' @param z (numeric or complex) vector of length one at which to evaluate the rational matrix.
#' @param ... optional additional parameters
#'
#' @return Matrix of dimensions of input object
#' 
#' @export
#' @rdname zvalue
zvalue = function(obj, z, ...){
  UseMethod("zvalue", obj)
}

#' @rdname zvalue
#' @export
zvalue.lpolm = function(obj, z = NULL, ...) {
  
  stopifnot(length(z) == 1)
  stopifnot(inherits(obj, "ratm"))
  
  # Integer-valued parameters
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] # Dimensions are retrieved from polm object, NOT from unclassed polm obj = array
  min_deg = d[4]
  
  obj = obj %>% unclass()
  
  # Check if a dimension or number of frequencies is zero (degree not part of this in contrast to ___.polm() function)
  if ((m*n) == 0) {
    return(matrix(0, m, n))
  }
  
  # Copy highest coefficient in a matrix as initial value
  out = matrix(obj[,, p+1-min_deg], m, n)
  
  if (p == 0 && min_deg == 0) {
    return(out)
  }
  
  # Horner-scheme for evaluation
  for (j in ((p-min_deg):1)) {
    out = out*z + matrix(obj[,,j], m, n)
  }
  
  # Every polynomial (i,j) needs to be premultiplied by z_0^{min_deg} (where z_0 denotes one particular point of evaluation)
  factor_min_deg = z^min_deg
  return(out * factor_min_deg)
}

#' @rdname zvalue
#' @export
zvalue.polm = function(obj, z = NULL, ...) {
  
  stopifnot(length(z) == 1)
  stopifnot(inherits(obj, "ratm"))
  
  # Integer-valued parameters
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] # Dimensions are retrieved from polm object, NOT from unclassed polm obj = array

  obj = obj %>% unclass()
  
  # Check if a dimension or number of frequencies is zero (degree not part of this in contrast to ___.polm() function)
  if ((m*n*(p+1)) == 0) {
    return(matrix(0, m, n))
  }
  
  # Copy highest coefficient in a matrix as initial value
  out = matrix(obj[,, p+1], m, n)
  
  if (p == 0) {
    return(out)
  }
  
  # Horner-scheme for evaluation
  for (j in (p:1)) {
    out = out*z + matrix(obj[,,j], m, n)
  }
  
  return(out)
}

#' @rdname zvalue
#' @export
zvalue.lmfd = function(obj, z = NULL, ...) {

  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if (p == -1) {
    stop("a(z) matrix is identically zero.")
  }
  
  if ((m*n*(q+1)) == 0) {
    return(matrix(0, m, n))
  }

  out_a = unclass(zvalue(obj$a, z = z))
  out_b = unclass(zvalue(obj$b, z = z))
  
  out_a_inv = try(solve(out_a))
  if(inherits(out_a_inv, "try-error")) {
    stop("a(z) not invertible for input z.")
  }
  
  return(out_a_inv %*% out_b)
}

#' @rdname zvalue
#' @export
zvalue.rmfd = function(obj, z = NULL, ...) {
  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  if (p == -1) {
    stop("c(z) matrix is identically zero.")
  }
  
  if ((m*n*(q+1)) == 0) {
    return(matrix(0, m, n))
  }
  
  out_c = unclass(zvalue(obj$c, z = z))
  out_d = unclass(zvalue(obj$d, z = z))
  
  out_c_inv = try(solve(out_c))
  if(inherits(out_c_inv, "try-error")) {
    stop("c(z) not invertible for input z.")
  }
  
  return(out_d %*% out_c_inv)
}

#' @rdname zvalue
#' @export
zvalue.stsp = function(obj, z = NULL,  ...) {

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  if ((m*n) == 0) {
    return(matrix(0, m, n))
  }
  
  ABCD = unclass(obj)
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]
  
  out = matrix(D, m, n)
  
  if (s == 0 || z == 0) {
    return(out)
  }

  z = z^{-1}
  x = try(C %*% solve(diag(z, nrow = s, ncol = s) - A, B), silent = TRUE)
  if (!inherits(x, 'try-error')) {
    out = out + x
  } else {
    stop("State space object has a pole at input z.")
  }
  
  return(out)
}
