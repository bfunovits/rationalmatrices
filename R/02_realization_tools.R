# ===================================================================
# Functional Area #2: Realization Algorithms & #6: State-Space Tools (CLAUDE.md)
#
# Purpose: Hankel matrix operations and echelon canonical forms
#   - Hankel matrix construction and analysis
#   - Echelon form basis computation
#   - Kronecker indices and nu2basis conversion
#   - Realization from power series using Hankel approach
#   - Row/column reduction algorithms
#
# Related Files: 06_statespace_methods.R (state-space), 01_representation_conversions.R (conversions)
# NOTE: Consider renaming to realization_10_utilities_tools.R to better reflect purpose
# ===================================================================

# Kronecker indices and echelon form ---------------------------------------------------------
#' Tools related to Kronecker indices
#'
#' The Kronecker indices \eqn{(\nu_1,\ldots,\nu_m)}{(nu[1],\ldots,nu[m])} describe
#' a (nice) basis of the row space of the Hankel matrix of the impulse response
#' coefficients. These indices are e.g. used to construct a (unique) LMFD representation
#' in echelon canonical form for a given impulse response.
#'
#' The function \code{pseries2nu} first constructs a Hankel matrix of the impulse response coefficients 
#' with \eqn{f} and \eqn{p} block rows and columns respectively. Then a (nice) basis for the row space is 
#' computed via a QR decomposition (with pivoting) of the transposed Hankel matrix. See \code{\link{qr}}.
#' If the size of the Hankel matrix is not specified (the parameter \code{Hsize=c(f,p)} is missing),  
#' then a default choice is made such that \eqn{f+p-1 = l}, \eqn{p\geq 1}{p \ge 1} and 
#' \eqn{f \geq p+1}{f \ge p+1} holds.
#' \cr
#' The function \code{pseries2nu} throws an error if the conditions \eqn{p\geq 1}{p \ge 1},  
#' \eqn{f \geq 2}{f \ge 2} and \eqn{l \geq f+p-1}{l \ge f+p-1} 
#' are not satisfied. In particular, this implies that 
#' \eqn{l \geq 2}{l \ge 2} must hold. 
#
#' @param basis s-dimensional (integer) vector which contains the indices of the
#'   basis rows of the Hankel matrix.
#' @param m (integer) the number of rows of the underlying rational matrix. 
#' @param nu vector with the Kronecker indices, i.e. an m-dimensional integer vector.
#' @param obj \code{\link{pseries}} object or 3-D array with dimension \eqn{(m,n,l+1)}. 
#'        This object represents the impulse response of a rational matrix.
#' @param Hsize integer vector \code{c(f,p)}, number of block rows and block columns 
#'        of the Hankel matrix of the impulse response coefficients. 
#'              If NULL a default choice is made. 
#' @param tol tolerance parameter, used by \code{\link{qr}}.
#'
#' @return
#' \item{basis2nu}{returns the Kronecker indices (\code{nu}) for given indices of the basis rows.}
#' \item{nu2basis}{returns the indices of the basis rows for given Kronecker indices \code{nu}.}
#' \item{pseries2nu}{determines the Kronecker indices for given impulse response  coefficients.
#' The function uses a QR decomposition (with pivoting) of the transpose of 
#' the Hankel matrix to determine the rank and the basis for its row space. See \code{\link{qr}}.}
#'
#' @examples
#' basis = c(1,2,3,4,6,7)     # suppose rows 1,2,3,4,6,7 of the Hankel matrix form a basis
#' nu = basis2nu(basis, m=3)  # compute the corresponding Kronecker index
#' print(nu)
#' all.equal(basis,nu2basis(nu))  # nu2basis(nu) returns the indices of the basis rows (basis)
#'
#' # generate random rational matrix in statspace form
#' # make sure that the matrix is stable
#' m = 3
#' n = 2
#' s = 7
#' A = matrix(rnorm(s*s), nrow = s, ncol = s)
#' A = A / (1.1 * max(abs(eigen(A, only.values = TRUE)$values)))
#' x = stsp(A, B = matrix(rnorm(s*n), nrow = s, ncol = n),
#'             C = matrix(rnorm(s*m), nrow = m, ncol = s),
#'             D = diag(1, nrow = m, ncol = n))
#' k = pseries(x, lag.max = 20)
#'
#' # compute the Kronecker indices of this  rational matrix
#' pseries2nu(k)
#'
#' \dontrun{
#' # Suppose the rational matrix has dimension m=2. Then the rows 1,2,5
#' # do not form a "nice" basis for the row space of the Hankel matrix.
#' # Therefore "basis2nu" stops with an error message.
#' basis2nu(c(1,2,5), m=2)
#' }
#'
#' @references
#' \insertRef{Hannan.Deistler12}{rationalmatrices}
#'
#' @name Kronecker-Indices
NULL

#' @rdname Kronecker-Indices
#' @export
basis2nu = function(basis, m) {
  # no basis elements => rank of H is zero
  if (length(basis)==0) {
    return(integer(m))
  }

  # What is the highest Kronecker indices?
  p = ceiling(max(basis)/m)

  # Create a matrix with one row for each variable.
  # The element in the j-th column is one if the (n*(j-1)+i)-th row of the Hankel matrix is in the basis
  in_basis = matrix(0, nrow=m, ncol=p)
  in_basis[basis] = 1

  # Calculate the Kronecker indices by summing across columns
  nu  = apply(in_basis, MARGIN=1, FUN=sum)

  # Check if there are holes in the basis, i.e. the (n*j+i)-th row is in the basis while while the (n*(j-1)+i)-th is not
  nu1 = apply(in_basis, MARGIN=1, FUN=function (x) sum(cumprod(x)))
  if (any (nu != nu1)){
    stop('This is not a nice basis, i.e. there are holes.')
  }

  return(nu)
}


#' @rdname Kronecker-Indices
#' @export
nu2basis = function(nu) {
  m = length(nu) # number of variables
  p = max(nu)    # largest Kronecker index

  # all Kronecker indices are equal to zero => rank is zero
  if (p==0) return(integer(0))

  # Create as many ones (boolean) in a row as the Kronecker index indicates
  # (apply returns columns and cbinds them, therefore t())
  in_basis = t( apply(matrix(nu), MARGIN = 1, FUN = function(x) c(rep(TRUE,x),rep(FALSE,p-x)) ) )

  # Where are the ones?
  basis = which(in_basis==1)
  return(basis)
}

#' Construct Hankel Matrix from Impulse Response Coefficients
#' 
#' This (internal) helper function builds a Hankel matrix from a 
#' given impulse response of dimension \eqn{(m,n)} with 
#' \eqn{l} lags. If the parameter \code{Hsize=c(f,p)} is not given, 
#' then a default choice for the number of block rows (\eqn{f}) and 
#' block columns (\eqn{p}) is made such that 
#' \eqn{f+p-1 = l}, \eqn{p\geq 1}{p \ge 1} and 
#' \eqn{f \geq p+1}{f \ge p+1}.
#' \cr
#' The function throws an error if the conditions \eqn{p\geq 1}{p \ge 1},  
#' \eqn{f \geq 2}{f \ge 2} and \eqn{l \geq f+p-1}{l \ge f+p-1} 
#' are not satisfied. In particular, this implies that 
#' \eqn{l \geq 2}{l \ge 2} must hold. 
#'
#' @param obj \code{\link{pseries}} object or 3-D array with dimension \eqn{(m,n,l+1)}.
#' @param Hsize integer vector \code{c(f,p)}, number of block rows and block columns. 
#'              If NULL a default choice is made. 
#' @return Block Hankel matrix with dimension \eqn{(fm, pn)} and 
#'         attributes \code{order=c(m,n,f,p)} and \code{k0}. The \eqn{(m,n)}-dimensional 
#'         matrix \code{k0} is the lag zero coefficient of the impulse response.
#' 
#' @export
#' @keywords internal
pseries2hankel = function(obj, Hsize = NULL) {
  
  # check input parameter "obj"
  k = unclass(obj)
  if ( (!(is.numeric(k) || is.complex(k))) || (!is.array(k)) || (length(dim(k)) !=3) ) {
    stop('parameter "obj" must be an "pseries" object or a 3-D array')
  }
  d = dim(k)
  m = d[1]
  n = d[2]
  lag.max = d[3] - 1
  if (lag.max < 2) {
    stop('the impulse response contains less than 2 lags')
  }
  
  # check size of Hankel matrix
  if (is.null(Hsize)) {
    Hsize = ceiling((lag.max+2)/2)
  }
  Hsize = as.vector(as.integer(Hsize))
  if (length(Hsize) == 1) {
    Hsize[2] = lag.max + 1 - Hsize
  }
  f = Hsize[1] # number of block rows of the Hankel matrix <=> future
  p = Hsize[2] # number of block columns of the Hankel matrix <=> past
  # the default choice is:
  # (f+p) = lag.max + 1 and f >= p+1 
  
  if ((f < 2) || (p < 1) || (lag.max < (f+p-1)) ) {
    stop('the conditions (f>1), (p>0) and (lag.max >= f+p-1) are not satisfied')
  }
  
  k0 = matrix(k[,,1], nrow = m, ncol = n)
  
  # Hankel matrix
  H = bhankel(k[,,-1,drop=FALSE], d = c(f,p))
  attr(H,'order') = c(m,n,f,p)
  attr(H,'k0') = k0
  
  return(H)
}    

#' Compute left-Kronecker Indices for given Hankel Matrix
#' 
#' This (internal) helper functions determines the left-Kronecker indices 
#' given the Hankel matrix of the impulse response coefficients. 
#' There are no checks on parameters! In particular, note that 
#' the Hankel matrix must have an attribute \code{order=c(m,n,f,p)} 
#' which describes the block size \eqn{(m,n)} and the number 
#' of block rows (\eqn{f}) and block columns (\eqn{p}).
#' 
#' @param H Block Hankel matrix, as computed e.g. by \code{\link{pseries2hankel}}.
#' @param tol tolerance parameter, used by \code{\link{qr}}.
#'
#' @return Integer vector with the Kronecker indices.
#' 
#' @export
#' @keywords internal
hankel2nu = function(H, tol = sqrt(.Machine$double.eps)) {
  
  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  p = d[4]
  
  if ((m*n) == 0) return(integer(m))
  
  # compute 'nice' basis for row space of H, via QR decomposition
  # of the transposed matrix.
  qr.H = qr(t(H), tol = tol)
  
  # rank of H is zero!
  if (qr.H$rank ==0) {
    nu = integer(m)
  } else {
    basis = qr.H$pivot[1:qr.H$rank]
    nu = try(basis2nu(basis, m))
    if (inherits(nu,'try-error')) {
      stop(paste(paste(basis,  collapse=' '),' is not a nice basis for the (',
                 m*f,',',n*p,') Hankel matrix', sep = ''))
    } 
  }
  
  return(nu)
}    

#' Compute right-Kronecker Indices for given Hankel Matrix
#' 
#' This (internal) helper functions determines the right-Kronecker indices 
#' given the Hankel matrix of the impulse response coefficients. 
#' There are no checks on parameters! In particular, note that 
#' the Hankel matrix must have an attribute \code{order=c(m,n,f,p)} 
#' which describes the block size \eqn{(m,n)} and the number 
#' of block rows (\eqn{f}) and block columns (\eqn{p}).
#' 
#' @param H Block Hankel matrix, as computed e.g. by \code{\link{pseries2hankel}}.
#' @param tol tolerance parameter, used by \code{\link{qr}}.
#'
#' @return Integer vector with the right-Kronecker indices.
#' 
#' @export
#' @keywords internal
hankel2mu = function(H, tol = sqrt(.Machine$double.eps)) {
  
  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  p = d[4]
  
  if ((m*n) == 0) return(integer(n))
  
  # compute 'nice' basis for column space of H, via QR decomposition
  qr.H = qr(H, tol = tol)
  
  if (qr.H$rank == 0) {   # rank of H is zero!
    mu = integer(n)
  } else {
    basis = qr.H$pivot[1:qr.H$rank]
    mu = try(basis2nu(basis, n)) # basis2nu does what would be needed for basis2mu, so no additional function is necessary
    if (inherits(mu,'try-error')) {
      stop(paste(paste(basis,  collapse=' '),' is not a nice basis for the (',
                 m*f,',',n*p,') Hankel matrix', sep = ''))
    } 
  }
  
  return(mu)
}    


#' @rdname Kronecker-Indices
#' @export
pseries2nu = function(obj, Hsize = NULL, tol = sqrt(.Machine$double.eps)) {
  
  # call the helper functions pseries2hankel and hankel2nu
  nu = try(hankel2nu(pseries2hankel(obj, Hsize = Hsize), tol = tol))
  if (inherits(nu,'try-error')) {
    stop('computation of Kronecker indices failed')
  }
  return(nu)
}


# internal function
# return statespace realization template 
#' @keywords internal
nu2stsp_template = function(nu, D) {
  m = nrow(D)
  n = ncol(D)
  s = sum(nu)
  if ((length(nu) != m) || ( (s > 0) && (n==0) )) {
    # notw: n = 0 implies nu = (0,...,0)!
    stop('the parameters "nu" and "dim" are not compatible')
  }

  if (s == 0) {
    A = matrix(0, nrow = 0, ncol = 0)
    B = matrix(0, nrow = 0, ncol = n)
    C = matrix(0, nrow = m, ncol = 0)
    return(stsp(A, B, C, D))
  }
  
  basis = nu2basis(nu)
  AC = matrix(0, nrow = s + m, ncol = s)
  dependent = c(basis + m, 1:m)
  for (i in (1:length(dependent))) {
    d = abs(basis-dependent[i])
    if (min(d) == 0) {
      # dependent[i]-th row is in basis
      j = which(d == 0)
      AC[i,j] = 1
    } else {
      j = which(basis < dependent[i])
      AC[i,j] = NA_real_
    }
  }
  A = AC[1:s,,drop = FALSE]
  C = AC[(s+1):(s+m), , drop = FALSE]
  B = matrix(NA_real_, nrow = s, ncol = n)
  return(stsp(A, B, C, D))
}
