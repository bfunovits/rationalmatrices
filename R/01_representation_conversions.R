# ===================================================================
# Functional Area #1: Representation Management (CLAUDE.md)
#
# Purpose: Convert between all rational matrix representations
#   - as.polm(), as.lpolm(): Convert to polynomial forms
#   - as.lmfd(), as.rmfd(): Convert to matrix fraction descriptions
#   - as.stsp(): Convert to state-space form
#   - as.pseries(), as.zvalues(): Convert to power series and frequency values
#
# Related Files: 01_representation_classes.R (constructors), 03_arithmetic_operations.R (upgrade logic)
# ===================================================================

# Conversion between polm() and lpolm() ####

#' Coerece to polynomial object
#' 
#' @param obj \code{\link{lpolm}} object
#' @param ... other arguments 
#'
#' @return \code{\link{polm}} object
#' @export
#'
#' @examples
#' lp = lpolm(1:3, min_deg = 1) 
#' as.polm(lp)
as.polm = function(obj, ...){
  UseMethod("as.polm", obj)
}

#' @rdname as.polm
#' @export
as.polm.lpolm = function(obj, ...){
  obj = prune(obj)
  min_deg = attr(obj, "min_deg")
  stopifnot("The *min_deg* attribute needs to be non-negative for coercion to polm-obj! Use get_bwd() for discarding negative powers." = min_deg >= 0)
  polm_offset = array(0, dim = c(dim(obj)[1:2], min_deg))
  return(polm(dbind(d = 3, polm_offset, unclass(obj))))
}

#' Coerce to Laurent polynom object
#' 
#' The attribute \code{min_deg} is set to zero for the given function input.
#'
#' @param obj \code{\link{polm}} object
#' @param ... other arguments
#'
#' @return \code{\link{lpolm}} object
#' @export
#'
#' @examples
#' p = test_polm(degree = 2)
#' as.lpolm(p)
as.lpolm = function(obj, ...){
  UseMethod("as.lpolm", obj)
}

#' @rdname as.lpolm
#' @export
as.lpolm.polm = function(obj, ...){
  attr(obj, "min_deg") = 0
  class(obj)[1] = "lpolm"
  return(obj)
}



# as.stsp.____ methods ##################################################

#' Coerce to Statespace Realization
#'
#' The function \code{\link{as.stsp.pseries}} calls
#' \code{\link{pseries2stsp}} with default parameters.
#' Of course the \code{\link{pseries}} object must contain
#' sufficiently many lags.
#'
#'
#' @param obj object
#' @param method character string
#' @param ... optional additional parameters
#'
#' @return object of class \code{\link{stsp}}.
#'
#' @rdname as.stsp
#' @export
as.stsp = function(obj, ...){
  UseMethod("as.stsp", obj)
}

#' @rdname as.stsp
#' @export
as.stsp.lpolm = function(obj, ...){
  stop("A lpolm object cannot be coerced to a state space model.")
}

#' @rdname as.stsp
#' @export
as.stsp.polm = function(obj, ...){
  obj = unclass(obj)
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if (p < 0) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n))
    return(x)
  }
  if (p == 0) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(obj, nrow = m, ncol = n))
    return(x)
  }
  if (m >= n) {
    x = stsp(A = rbind(matrix(0, nrow = n, ncol = p*n), diag(x = 1, nrow = (p-1)*n, ncol = p*n)),
             B = diag(x = 1, nrow = p*n, ncol = n),
             C = matrix(obj[,,-1], nrow = m, ncol = p*n), D = matrix(obj[,,1], nrow = m, ncol = n))
    return(x)
  }
  B = obj[,,-1,drop = FALSE]
  B = aperm(B, c(1,3,2))
  dim(B) = c(p*m, n)
  x = stsp(A = cbind(matrix(0, nrow = p*m, ncol = m), diag(x = 1, nrow = p*m, ncol = (p-1)*m)),
           B = B, C = diag(x = 1, nrow = m, ncol = p*m), D = matrix(obj[,,1], nrow = m, ncol = n))
  return(x)
}

#' @rdname as.stsp
#' @export
as.stsp.lmfd = function(obj, ...){
  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  # note for a valid lmfd object, m > 0, p >= 0 must hold!
  if ((m*(p+1)) == 0) stop('input object is not a valid "lmfd" object')

  if ( (n*(q+1)) == 0 ) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n))
    return(x)
  }

  ab = unclass(obj)
  a = ab[,1:(m*(p+1)), drop = FALSE]
  b = ab[,(m*(p+1)+1):(m*(p+1)+n*(q+1)), drop = FALSE]

  # check a(z)
  a0 = matrix(a[, 1:m, drop = FALSE], nrow = m, ncol = m)
  junk = try(solve(a0))
  if (inherits(junk, 'try-error')) stop('left factor "a(0)" is not invertible')

  if ((p == 0) && (q == 0)) {
    # static system
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = solve(a0, matrix(b, nrow = m, ncol = n)))
    return(x)
  }

  # a[,,i] -> a0^{-1} a[,,i], convert to matrix and append zeroes if p < q
  if (p > 0) {
    a = solve(a0, a[, (m+1):(m*(p+1)), drop = FALSE])
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }
  if (p < q) a = cbind(a, matrix(0, nrow = m, ncol = (q-p)*m))
  p = max(p,q)

  # compute impulse response
  # this is not very efficient,
  # e.g. the scaling of the coefficients a0^(-1)a[,,i] and a0^(-1)b[,,i] is done twice
  k = unclass(pseries(obj, lag.max = p))

  A = rbind(-a, diag(x = 1, nrow = m*(p-1), ncol = m*p))
  D = matrix(k[,,1], nrow = m, ncol = n)
  C = cbind(matrix(0, nrow = m, ncol = (p-1)*m), diag(m))
  B = aperm(k[,,(p+1):2,drop = FALSE], c(1,3,2))
  dim(B) = c(p*m, n)
  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}


#' @rdname as.stsp
#' @export
#' 
#' @examples
#' (rr = test_rmfd(dim = c(3,2), degrees = c(2,1)))
#' as.stsp(rr)
as.stsp.rmfd = function(obj, ...){
  d = attr(obj, 'order')
  m = d[1] # output dimension
  n = d[2] # input dimension
  p = d[3] # degree of c(z)
  q = d[4] # degree of d(z)

  # otherwise c(z) is not invertible
  stopifnot("Input object is not a valid rmfd object" = (m*(p+1)) > 0)

  # if d(z) is identically zero or has no column, return early
  if ( (n*(q+1)) == 0 ) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n))
    return(x)
  }

  # c(0) must be the identity matrix. Note that this trafo only changes the covariance matrix of the inputs
  c0 = matrix(c(unclass(obj$c))[1:(n^2)], nrow = n, ncol = n)
  c0inv = tryCatch(solve(c0),
                   error = function(cnd) stop(' "c(0)" is not invertible'))
  cd = unclass(obj) %*% c0inv
  
  # Separate c(z) and d(z)
  c = cd[1:(n*(p+1)),, drop = FALSE]
  d = cd[(n*(p+1)+1):(n*(p+1)+m*(q+1)),, drop = FALSE]

  # Case 1: Static System
  if ((p == 0) && (q == 0)) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(d, nrow = m, ncol = n))
    return(x)
  }

  # Case 2: (p>0 or q>0): Construct state space system  
  if (p < q+1){
    c = rbind(c,
              matrix(0, nrow = (q+1-p)*n, ncol = n))
  }
  p = max(p,q+1)

  # Transpose and reshuffle c(z) such that the coefficients are in a wide matrix, then create stsp matrices (A,B)
  c = t(c)[,-(1:n), drop = FALSE] %>% array(dim = c(n,n,p)) %>% aperm(perm = c(2,1,3)) %>% matrix(nrow = n)

  A = rbind(-c,
            diag(x = 1, nrow = n*(p-1), ncol = n*p))
  B = rbind(c0inv,
            matrix(0, nrow = n*(p-1), ncol = n))

  # Same for d(z), and add zeros to d(z) if p>q+1, then create stsp matrices (C,D) 
  d = t(d) %>% array(dim = c(n,m,q+1)) %>% aperm(perm = c(2,1,3)) %>% matrix(nrow = m)
  C = cbind(d, matrix(0, nrow = m, ncol = n*(p-q-1)))
  D = C %*% B
  C = C %*% A

  # Create output
  x = stsp(A = A, B = B,
           C = C, D = D)
  return(x)
}




#' Ho-Kalman Realization Algorithm
#'
#' This helper function implements the Ho-Kalman algorithm.
#'
#' The procedure(s) may be used for model reduction (with some care).
#'
#' There are a number of restrictions on the number of lags \eqn{l}
#' of the impulse response, the number of block rows (\eqn{f}),
#' block columns (\eqn{p}) of the Hankel matrix and the Kronecker
#' indices \eqn{\nu_i}{\nu[i]}. We require that:
#' \eqn{p>0}, \eqn{f>1},
#' \eqn{l \geq f+p-1}{l\ge f+p-1} and \eqn{\nu_i <f}{\nu[i] < f}.
#' If these restrictions are not satisfied an error is thrown.
#'
#' @param obj \code{\link{pseries}} object or 3-D array with dimension \eqn{(m,n,l+1)}.
#' @param method Character string, which determines the method and the "parametrization" type
#'             of the state space model. See below for more details.
#' @param Hsize integer vector \code{c(f,p)}, number of block rows and block columns
#'              of the Hankel matrix which is used to construct the statespace realization.
#'              If NULL a default choice is made.
#' @param s desired state dimension. Only used for \code{method = "balanced"}.
#'   Note however, if \eqn{s} is larger than the rank of the Hankel matrix,
#'   then the procedure will break down.
#'   If \code{s} is missing, then the state dimension is determined from
#'   the singular values of the Hankel matrix.
#'   To be precise the state dimension is chosen as the number of singular values which are
#'   greater than or equal to \code{tol} times the maximum singular value.
#' @param nu Kronecker indices. Only used for  \code{method = "echelon"}.
#'   If missing, then \code{nu} is computed with a QR decomposition of the
#'   transpose of the Hankel matrix of the impulse response coefficients.
#' @param tol tolerance parameter used for the QR decomposition or the SVD decomposition
#'            of the Hankel matrix \eqn{H} of the impulse response coefficients.
#' @param Wrow,Wcol weighting matrices (default is no weighting, i.e. identity matrices).
#'                  These weighting matrices are only used for \code{method="balanced"}, where the
#'                  SVD of the weighted Hankel matrix \code{Wrow \%*\% H \%*\% t(Wcol)} is computed.
#'
#' @aliases Ho-Kalman
#'
#' @return List with slots
#' \item{Xs}{\code{\link{stsp}} object, the rational matrix in statespace form}
#' \item{Hsv}{Singular values of the Hankel matrix for \code{method='balanced'} and \code{NULL} else.}
#' \item{nu}{Kronecker indices for \code{method='echelon'} and \code{NULL} else.}
#' @export
#'
#' @examples
#' # generate random rational matrix X(z) in statespace form
#' # make sure that the A matrix is stable
#' m = 3
#' n = 2
#' s = 7
#' A = matrix(rnorm(s*s), nrow = s, ncol = s)
#' A = A / (1.1 * max(abs(eigen(A, only.values = TRUE)$values)))
#' Xs = stsp(A, B = matrix(rnorm(s*n), nrow = s, ncol = n),
#'           C = matrix(rnorm(s*m), nrow = m, ncol = s),
#'           D = diag(1, nrow = m, ncol = n))
#' Xi = pseries(Xs, lag.max = 20)
#'
#' out = pseries2stsp(Xi, method = 'balanced')
#' print(out)
#' # check impulse response
#' all.equal(pseries(out$Xs, lag.max = 20), Xi)
#'
#' Xs1 = as.stsp(Xi)
#' all.equal(Xs1, out$Xs)
#'
#' out = pseries2stsp(Xi, method = 'echelon')
#' print(out)
#' # check impulse response
#' all.equal(pseries(out$Xs, lag.max = 20), Xi)
#'
#' Xs1 = as.stsp(Xi, method = 'echelon')
#' all.equal(Xs1, out$Xs)
pseries2stsp = function(obj, method = c('balanced', 'echelon'),
                        Hsize = NULL, s = NULL, nu = NULL,
                        tol = sqrt(.Machine$double.eps), Wrow = NULL, Wcol = NULL) {
  # no input checks

  method = match.arg(method)

  # construct Hankel matrix with the helper function pseries2hankel
  H = try(pseries2hankel(obj, Hsize = Hsize))
  if (inherits(H, 'try-error')) {
    stop('computation of Hankel matrix failed')
  }

  k0 = attr(H, 'k0')
  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  p = d[4]

  # take care of the case of an empty impulse response (m*n=0)
  if ((m*n) == 0) {
    s = 0
    Xs = stsp(A = matrix(0, nrow = s, ncol = s), B = matrix(0, nrow = s, ncol = n),
              C = matrix(0, nrow = m, ncol = s), D = matrix(0, nrow = m, ncol = n))
    return(list(Xs = Xs, Hsv = numeric(0), nu = integer(0)))
  }

  # compute statespace model in "balanced" form ###############
  if (method == 'balanced') {
    if (!is.null(Wrow)) {
      H = Wrow %*% H
    }
    if (!is.null(Wcol)) {
      H = H %*% t(Wcol)
    }

    svd.H = svd(H)
    Hsv = svd.H$d    # singular values of (weighted) Hankel matrix

    if (is.null(s)) {
      # determine state dimension from singular values
      s = ifelse(svd.H$d[1]>.Machine$double.eps, sum(svd.H$d >= (tol*svd.H$d[1])),0)
    }

    if (s>0) {
      if (s > m*(f-1)) {
        stop('number of block rows of "H" is too small for the (desired) state dimension "s"!')
      }
      sv2 = sqrt(Hsv[1:s])
      # (approximately) factorize H as H = U V
      U = svd.H$u[,1:s,drop = FALSE] * matrix(sv2, nrow = nrow(H), ncol = s, byrow = TRUE)
      V = t( svd.H$v[,1:s,drop = FALSE] * matrix(sv2, nrow = ncol(H), ncol = s, byrow = TRUE) )
      if (!is.null(Wrow)) {
        U = solve(Wrow, U)
      }
      if (!is.null(Wcol)) {
        V = t( solve(Wcol, t(V)) )
      }
      C = U[1:m,,drop = FALSE]
      B = V[,1:n,drop = FALSE]
      A = stats::lsfit(U[1:(m*(f-1)), ,drop = FALSE],
                       U[(m+1):(m*f),,drop = FALSE], intercept = FALSE)$coef
      A = unname(A)
    } else {
      A = matrix(0, nrow=s,  ncol=s)
      B = matrix(0, nrow=s,  ncol=n)
      C = matrix(0, nrow=m, ncol=s)
    }
    D = k0

    Xs = stsp(A = A, B = B, C = C, D = D)

    return( list(Xs = Xs, Hsv = Hsv, nu = NULL) )
  } # method == balanced

  # compute statespace model in "echelon" form ###############

  # if (n < m) {
  #   stop('method "echelon" for the case (n < m) is not yet implemented!')
  # }

  # compute Kronecker indices
  if (is.null(nu)) {
    nu = try(hankel2nu(H, tol = tol))
    if (inherits(nu, 'try-error')) stop('computation of Kronecker indices failed')
  }
  # check nu
  nu = as.vector(as.integer(nu))
  if ( (length(nu) != m) || (min(nu) < 0) || (max(nu) > (f-1)) ) {
    stop('Kronecker indices are not compatible with the impulse response')
  }
  basis = nu2basis(nu)
  s = length(basis) # state space dimension

  # consider the transposed Hankel matrix!
  H = t(H)

  if (s > 0) {
    AC = matrix(0, nrow = s+m, ncol = s)
    dependent = c(basis + m, 1:m)
    # cat('basis', basis,'\n')
    for (j in (1:length(dependent))) {
      i = dependent[j]
      if (min(abs(i - basis)) == 0) {
        # row i is in the basis!
        k = which(i == basis)
        AC[j, k] = 1
      } else {
        # explain row i by the previous basis rows
        k = which(basis < i)
        AC[j, k] = stats::lsfit(H[ , k, drop=FALSE], H[ , i], intercept=FALSE)$coef
      }
      # cat('j=', j, 'i=', i, 'k=', k, 'basis[k]=', basis[k], '\n')
    }
    B = t(H[(1:n), basis, drop=FALSE])
    A = AC[1:s,,drop = FALSE]
    C = AC[(s+1):(s+m), , drop = FALSE]
  } else {
    A = matrix(0, nrow=s,  ncol=s)
    B = matrix(0, nrow=s,  ncol=n)
    C = matrix(0, nrow=m, ncol=s)
  }
  D = k0

  Xs = stsp(A = A, B = B, C = C, D = D)

  return( list(Xs = Xs, Hsv = NULL, nu = nu) )
}


#' @rdname as.stsp
#' @export
as.stsp.pseries = function(obj, method = c('balanced','echelon'), ...){
  out = pseries2stsp(obj, method = method)
  return(out$Xs)
}

# as.lmfd.____ methods ##################################################


#' Construct a LMFD Representation from Impulse Response
#'
#' This (helper) function constructs a left matrix fraction description
#' of a rational matrix from its impulse response. 
#' Of course the impulse response must contain sufficiently many lags. 
#' The constructed LMFD is in \emph{echelon canonical} form.
#'
#' There are a number of restrictions on the dimension \eqn{(m,n)}
#' and number of lags \eqn{l} of the impulse response,
#' the number of block rows (\eqn{f}), block columns (\eqn{p})
#' of the Hankel matrix and the Kronecker indices \eqn{\nu_i}{\nu[i]}:
#'
#' We require that:
#' \eqn{m>0}, \eqn{p>0}, \eqn{f>1},
#' \eqn{l \geq f+p-1}{l\ge f+p-1} and \eqn{\nu_i <f}{\nu[i] < f}.
#' If these restrictions are not satisfied an error is thrown.
#'
#' @param obj \code{\link{pseries}} object or 3-D array with dimension \eqn{(m,n,l+1)}.
#' @param Hsize integer vector \code{c(f,p)}, number of block rows and block columns
#'              of the Hankel matrix which is used to construct the LMFD.
#'              If NULL a default choice is made.
#' @param nu integer vector with the Kronecker indices. If \code{NULL} then the
#'           Kronecker indices are computed via a QR decomposition of the
#'           transpose of the Hankel matrix, see \code{\link{qr}}.
#' @param tol tolerance parameter, used by \code{\link{qr}}.
#'
#' @return List with components
#' \item{Xl}{\code{\link{lmfd}} object which contains the LMFD representation
#'          of the rational matrix (in echelon form).}
#' \item{nu}{integer vector with the Kronecker indices.}
#' @export
#'
#' @examples
#' # generate a random LMFD object
#' m = 3
#' n = 2
#' p = 1
#' q = 1
#' a = test_polm(dim = c(m,m), degree = p, random = TRUE, bzeroes = 1)
#' b = test_polm(dim = c(m,n), degree = q, random = TRUE, bzeroes = 1)
#' Xl = lmfd(a,b)
#'
#' # compute impulse response of this matrix
#' Xi = pseries(Xl, lag.max = 2*(max(p,q)+1))
#'
#' # reconstruct a matrix from this impulse response
#' out = pseries2lmfd(Xi)
#' out
#'
#' # check that the lmfd object is a realization of the given impulse response
#' Xi1 = pseries(out$Xl, lag.max = 2*(max(p,q)+1))
#' all.equal(Xi, Xi1)
pseries2lmfd = function(obj, Hsize = NULL, nu = NULL, tol = sqrt(.Machine$double.eps)) {


  # construct Hankel matrix with the helper function pseries2hankel
  H = try(pseries2hankel(obj, Hsize = Hsize))
  if (inherits(H, 'try-error')) {
    stop('computation of Hankel matrix failed')
  }

  k0 = attr(H, 'k0')
  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  # p = d[4]

  if (m == 0) stop('the number of rows (m) must be positive!')

  # compute Kronecker indices
  if (is.null(nu)) {
    nu = try(hankel2nu(H, tol = tol))
    if (inherits(nu, 'try-error')) stop('computation of Kronecker indices failed')
  }
  # check nu
  nu = as.vector(as.integer(nu))
  if ( (length(nu) != m) || (min(nu) < 0) || (max(nu) > (f-1)) ) {
    stop('Kronecker indices are not compatible with the impulse response')
  }

  k = unclass(obj)
  H = t(H) # use the transposed of H in the following!

  p = max(nu)

  # a(z),b(z) have degree zero => c(z) is constant
  if (p == 0) {
    Xl = lmfd(b = k0)
    return(list(Xl = Xl, nu = nu))
  }

  # Note that the block diagonal is zero because 
  # the element R[,,1] in the argument R of btoeplitz() 
  # overwrites the block diagonal element from argument C 
  # (corresponding to the zero lag coefficient of the transfer function)
  Tk = btoeplitz(R = array(0, dim = c(m, n, p)), C = k[, , (1:(p+1)), drop=FALSE])
  # print(Tk)

  basis = nu2basis(nu)
  # index of the "dependent" rows of H
  dependent = m*nu + (1:m)
  # matrices with coefficients in reversed order
  a = matrix(0, nrow = m, ncol = m*(p+1))  # [a[p], a[p-1], ..., a[0]]
  # note b relates to b(z) - b(0) = b(z) - a(0)*k(0)
  b = matrix(0, nrow = m, ncol = n*(p+1))  # [b[p], b[p-1], ..., b[0]]
  for (i in (1:m)) {
    j = basis[basis < dependent[i]]
    # cat('i=', i,', dependent=', dependent[i], ', j=', j,'\n',sep = ' ')
    if (length(j) > 0) {
      ai = stats::lsfit(H[ , j,drop = FALSE], H[ , dependent[i]], intercept = FALSE)$coef
      a[i, j + m*(p-nu[i])] = -ai
    }
    a[i, dependent[i] + m*(p-nu[i])] = 1
    # print(a)
    # note b(0) = 0
    # cat(i,':',nu[i],':',p,':',iseq(1 + n*(p-nu[i]), n*p),'\n')
    b[i, iseq(1 + n*(p-nu[i]), n*p)] =
      a[i,] %*% Tk[, iseq(1 + n*(p-nu[i]), n*p), drop = FALSE]
    # print(b)
  }
  # print(t(H))
  # print(a)
  # print(a %*% t(H)[1:ncol(a),,drop = FALSE])
  # print(b)
  dim(a) = c(m, m, p+1)
  dim(b) = c(m, n, p+1)

  # ak0 <=> a(z)*k(0)
  ak0 = aperm(a, c(1,3, 2))
  dim(ak0) = c(m*(p+1), m)
  ak0 = ak0 %*% k0
  dim(ak0) = c(m, p+1, n)
  ak0 = aperm(ak0, c(1, 3, 2))

  # b -> b + a(z)*k(0)
  b = b + ak0

  # reshuffle
  a = a[,,(p+1):1]
  b = b[,,(p+1):1]

  Xl = lmfd(polm(a), polm(b))
  return(list(Xl = Xl, nu = nu))
}


#' Coerce to Left Matrix Fraction Description
#'
#' The function \code{\link{as.lmfd.pseries}} calls
#' \code{\link{pseries2lmfd}} with default parameters.
#' Of course the \code{\link{pseries}} object must contain
#' sufficiently many lags.
#'
#' @param obj object
#' @param method character string
#' @param ... optional additional arguments
#'
#' @return object of class \code{\link{lmfd}}
#'
#' @rdname as.lmfd
#' @export
as.lmfd = function(obj, method, ...){
  UseMethod("as.lmfd", obj)
}

#' @rdname as.lmfd
#' @export
as.lmfd.pseries = function(obj, method, ...){
  return(pseries2lmfd(obj)$Xl)
}

# as.rmfd.____ methods ####

#' Construct an RMFD Representation from Impulse Response
#'
#' This (helper) function constructs a right matrix fraction description of a rational matrix from its impulse response.
#' Of course, the impulse response must contain sufficently many lags.
#' The constructed RMFD is in canonical form.
#'
#' There are a number of restrictions on the dimension \eqn{(m,n)} and number of lags \eqn{l} of the impulse response,
#' the number of block rows (\eqn{f}), block columns (\eqn{p}) of the Hankel matrix and the right-Kronecker indices \eqn{\mu_i}{\mu[i]}:
#'
#' We require that:
#' \itemize{
#' \item \eqn{m>0}, \eqn{p>0}, \eqn{f>1},
#' \item \eqn{l \geq f+p-1}{l\ge f+p-1} and
#' \item \eqn{\nu_i <f}{\nu[i] < f}.
#' }
#' If these restrictions are not satisfied an error is thrown.
#'
#' @param obj \code{\link{pseries}} object or 3-D array with dimension \eqn{(m,n,l+1)}.
#' @param Hsize integer vector \code{c(f,p)}, number of block rows and block columns of the Hankel matrix which is used to construct the RMFD.
#'     If NULL a default choice is made.
#' @param mu integer vector with the right-Kronecker indices.
#'     If \code{NULL} then the right-Kronecker indices are computed via a QR decomposition of the transpose of the Hankel matrix, see \code{\link{qr}}.
#' @param tol tolerance parameter, used by \code{\link{qr}}.
#'
#' @return List with components
#' \item{Xr}{\code{\link{rmfd}} object which contains the RMFD resresentation
#'          of the rational matrix (in echelon form).}
#' \item{mu}{integer vector with the right-Kronecker indices.}
#' @importFrom purrr map
#' @export
#'
#' @examples
#' # generate a random RMFD object
#' m = 3
#' n = 2
#' p = 1
#' q = 1
#' cc = test_polm(dim = c(n,n), degree = p, random = TRUE)
#' dd = test_polm(dim = c(m,n), degree = q, random = TRUE)
#' Xr = rmfd(cc,dd)
#'
#' # compute impulse response of this matrix
#' Xi = pseries(Xr, lag.max = 2*(max(p,q)+1))
#'
#' # reconstruct a matrix from this impulse response
#' out = pseries2rmfd(Xi)
#' out
#'
#' # check that the lmfd object is a realization of the given impulse response
#' Xi1 = pseries(out$Xr, lag.max = 2*(max(p,q)+1))
#' all.equal(Xi, Xi1)
pseries2rmfd = function(obj, Hsize = NULL, mu = NULL, tol = sqrt(.Machine$double.eps)) {

  # Construct Hankel matrix with the helper function pseries2hankel ####
  H = try(pseries2hankel(obj, Hsize = Hsize))
  if (inherits(H, 'try-error')) {
    stop('computation of Hankel matrix failed')
  }

  # __ Dimensions of Hankel ####
  k0 = attr(H, 'k0')
  H_order = attr(H, 'order')
  dim_out = H_order[1]
  dim_in = H_order[2]
  block_rows = H_order[3]
  block_cols = H_order[4] # this was commented out for pseries2lmfd, here we need it though

  if (dim_out == 0) stop('The number of rows (m) must be positive!')

  # Compute right-Kronecker indices mu ####
  if (is.null(mu)) {
    mu = try(hankel2mu(H, tol = tol))
    if (inherits(mu, 'try-error')){
      stop('Computation of Kronecker indices failed')
    } 
  }
  # check mu
  mu = as.vector(as.integer(mu))
  if ( (length(mu) != dim_in) || (min(mu) < 0) || (max(mu) > (block_cols-1)) ) {
    stop('Kronecker indices are not compatible with the impulse response')
  }

  # Transfer function as array (i.e. power series coefficients) ####
  k = unclass(obj)
  max_mu = max(mu) # difficulty in notation: p is used as deg(c(z)) or deg(a(z)), but also as number of block columns of the Hankel matrix

  # Special case: Zero poly (c(z), d(z) have degree zero => c(z) is constant ) ####
  if (max_mu == 0) {
    Xr = rmfd(d = k0)
    return(list(Xr = Xr, mu = mu))
  }

  # Needed to determine the coefficients in d(z) ####
  # Note how the coefficients are stacked when using this Toeplitz matrix! ([\tilde{d}[0]', \tilde{d}[1]', ..., \tilde{d}[p]']')
  Tk = btoeplitz(R = array(c(rep(0,dim_out*dim_in),
                           k[,,2:(max_mu+1)]), dim = c(dim_out, dim_in, max_mu+1)),
                 C = array(rep(0, dim_out*dim_in*(max_mu+1)), dim =c(dim_out, dim_in, max_mu+1)))

  # Which columns of Hankel span its column space? ####
  basis = nu2basis(mu) # nu2basis is okay also for right-Kronecker indices. No need to create new function mu2basis

  # First dependent column of Hankel for each variable (index of the "dependent" columns of H) ####
  dependent = dim_in*mu + (1:dim_in)

  # Tilde coefficient matrices ####
  # \tilde{c}(z) = \tilde{c_0} + \tilde{c_1} z + ... for describing \tilde{k}(z) = k_1 z^(-1) + k_2 z^(-2) + ... (see Hannan/Deistler Chapter 2.4 and 2.5)
  
  # [\tilde{c}[0]', \tilde{c}[1]', ..., \tilde{c}[p]']' in echelon form, but the columns get shifted immediately with z^(p-mu_i) (in the for-loop)!
  c = matrix(0, nrow = dim_in*(max_mu+1), ncol = dim_in)  
  # [\tilde{d}[0]', \tilde{d}[1]', ..., \tilde{d}[p]']': also shifted in the for-loop below
  d = matrix(0, nrow = dim_out*(max_mu+1), ncol = dim_in)  

  # Obtain the coefficients for each (dependent) column of Hankel ####
  for (ix_var in (1:dim_in)) {

    # indices of basis columns of Hankel to the left of the first dependent column pertaining to variable ix_var
    ix_basis = basis[basis < dependent[ix_var]]

    # regression coefficients
    if (length(ix_basis) > 0) {
      ci = stats::lsfit(H[ , ix_basis,drop = FALSE], H[ , dependent[ix_var]], intercept = FALSE)$coef #lsfit(x,y)
      c[ix_basis, ix_var] = -ci
    }
    # Coefficient of dependent column corresponding to variable ix_var is set to one
    c[dependent[ix_var], ix_var] = 1
    
    # Obtain \tilde{d}(z) from the Toeplitz matrix (corresponding to non-negative coefficients in comparison)
    # and the corresponding column in the matrix representation of \tilde{c}(z)
    d[, ix_var] =  Tk %*% c[, ix_var, drop = FALSE]
    
    # Multiply z^(max_mu - mu[ix_var]) on column "ix_var" of matrix representation of \tilde{c}(z) and \tilde{d}(z) 
    # such that \tilde{c}_p has ones on the diagonal and is upper triangular
    c[, ix_var] = c(rep(0, dim_in*(max_mu - mu[ix_var])),
                    c[1:(dim_in*(mu[ix_var]+1)), ix_var])
    d[, ix_var] = c(rep(0, dim_out*(max_mu - mu[ix_var])),
                    d[1:(dim_out*(mu[ix_var]+1)), ix_var])
    
  }

  # Transform \tilde{c}(z) to array ####
  c = t(c)
  dim(c) = c(dim_in, dim_in, max_mu+1)
  c = aperm(c, perm = c(2,1,3))

  # k0c <=> k(0)*c(z) (to be added to \tilde{d}(z)) ####
  k0c = purrr::map(1:(max_mu+1), ~ k0 %*% c[,,.x]) %>% unlist() %>% array(dim = c(dim_out, dim_in, max_mu+1))

  # Transform \tilde{d}(z) to array ####
  d = t(d)
  dim(d) = c(dim_in, dim_out, max_mu+1) # (dim_in, dim_out) because we work here with tranposed!
  d = aperm(d, perm = c(2,1,3))

  # Add k_0 * \tilde(c)(z) to \tilde{d}(z) ####
  d = d + k0c

  # Obtain representation in backward shift (instead of forward shift) ####
  c = c[,,(max_mu+1):1, drop = FALSE]
  d = d[,,(max_mu+1):1, drop = FALSE]

  Xr = rmfd(polm(c), polm(d))
  return(list(Xr = Xr, mu = mu))
}


#' Coerce to Right Matrix Fraction Description
#'
#' The function \code{\link{as.rmfd.pseries}} calls \code{\link{pseries2rmfd}} with default parameters.
#' Of course the \code{\link{pseries}} object must contain sufficiently many lags.
#'
#' @param obj object
#' @param method character string
#' @param ... optional additional arguments
#'
#' @return object of class \code{\link{rmfd}}
#'
#' @rdname as.rmfd
#' @export
as.rmfd = function(obj, method, ...){
  UseMethod("as.rmfd", obj)
}

#' @rdname as.rmfd
#' @export
as.rmfd.pseries = function(obj, method, ...){
  return(pseries2rmfd(obj)$Xr)
}
