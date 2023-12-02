# Essentials (computing zeros!) ####

#' Left Prime and Left Coprime Polynomials
#'
#' Check whether a polynomial is \emph{left prime} or a pair of two polynomials is \emph{left coprime}.
#' This check uses a (singular) pencil associated with the polynomial(s). For more details 
#' see the vignette \href{../doc/technical_details.html}{Technical Details}.
#'
#' @param a,b If \code{a} is an \code{\link{lmfd}} or an \code{\link{rmfd}} object, 
#'     which represents a left matrix fraction description, say \eqn{p^{-1}(z) q(z)}, 
#'     or a right MFD, say \eqn{r(z) s^{-1}(z)}, 
#'     then the procedure tests whether the pair \eqn{(p(z),q(z))} or \eqn{(t(r(z)),t(q(z)))}is left coprime. 
#'            
#'     Otherwise the arguments \code{a} and \code{b} (if \code{b} is not NULL) must represent two compatible polynomial matrices, 
#'     i.e. \code{a}, \code{b} must be \code{\link{polm}} objects (or objects which may be coerced to \code{polm} objects). 
#'     If \code{b} is NULL, the procedures checks whether \eqn{a(z)} is left prime,
#'     otherwise the pair \eqn{(a(z),b(z))} is checked for left coprimeness.
#' @param tol a tolerance parameter, which is used to decide the rank of certain matrices.
#' @param only.answer if TRUE, just return a logical (\code{TRUE} or \code{FALSE}). 
#'     Otherwise a list with additional information is returned.
#' @param debug if TRUE, print some diagnostic information.
#'
#' @note
#' This procedure returns different objects, depending on the parameter \code{only.answer}.
#'
#' @return If \code{only.answer} is true then a logical (\code{TRUE} or \code{FALSE}) is returned.
#'   Otherwise, a list with the following slots is returned. A more detailed description 
#'   of these items is given in the vignette \href{../doc/technical_details.html}{Technical Details}.
#'   \item{answer}{A boolean as above}
#'   \item{A,B}{These matrices represent the pencil \eqn{(A-Bz)} (in staircase form) which is used 
#'              to check the left (co-)prime condition.}
#'   \item{m,n}{Two integer vectors which code the structure of the staircase form.}
#'   \item{zeroes}{If available, a vector of zeroes of the matrix \eqn{(a(z),b(z))}. If 
#'                 \eqn{(a,b)} have no common zeroes (the left coprime case) then 
#'                 \code{zeroes} is an empty numeric vector. The case that 
#'                 \eqn{(a(z),b(z))} is rank deficient for \emph{all} \eqn{z \in C}{z in C} 
#'                 is coded with \code{z=NA}.}
#'   
#' @export
#'
#' @examples
#' # Ex 1: Two coprime polynomials ##################################################
#'
#' # Generate two random (2 x 2) polynomial matrices with degree 2
#' set.seed(1803)
#' a = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)
#' b = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)
#'
#' # Output: "only.answer = TRUE"
#' is.coprime(a, b, debug = FALSE, only.answer = TRUE)
#'
#' # Output: "only.answer = FALSE"
#' out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
#' str(out)
#' out$answer
#' out$zeroes
#' 
#' # we could equivalently use the syntax: 
#' is.coprime(cbind(a,b))
#' is.coprime(lmfd(a,b))
#'
#' # Ex 2: Two non-coprime polynomials with a finite number of common zeros #############
#' # Dimensions of a, b, and the common factor r
#' dim = 3
#' deg_aa = 1
#' deg_bb = 1
#' deg_r = 1
#'
#' # Generate random polynomial matrices
#' a0 = a
#' b0 = b
#' # generate common factor 
#' r = test_polm(dim = c(2,2), degree = 1, random = TRUE, digits = 1)
#'
#' # Generate polynomials with a common factor
#' a = r %r% a0
#' b = r %r% b0
#'
#' out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
#' out$answer
#' out$zeroes
#'
#'
#' # Ex 3: Two non-coprime polynomials: Everywhere rank deficient ###################
#'
#' # generate a common factor of rank 1 
#' r = test_polm(dim = c(2,1), degree = 1, random = TRUE, digits = 1) %r% 
#'     test_polm(dim = c(1,2), degree = 1, random = TRUE, digits = 1)
#'
#' # Rank deficient matrices with common factor
#' a = r %r% a0
#' b = r %r% b0
#'
#' out = is.coprime(a,b, only.answer = FALSE)
#' out$answer
#' out$zeroes
#' 
#' 
#' # Ex 4: Right-MFD ####
#' 
#' c = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)
#' d = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)
#'
#' # Output: "only.answer = TRUE"
#' is.coprime(t(c), t(d), debug = FALSE, only.answer = TRUE)
#'
#' # Output: "only.answer = FALSE"
#' out = is.coprime(t(c), t(d), debug = FALSE, only.answer = FALSE)
#' str(out)
#' out$answer
#' out$zeroes
#' 
#' # we could equivalently use the syntax: 
#' is.coprime(rbind(c,d))
#' is.coprime(rmfd(c,d))
#' 
#' # reset seed
#' set.seed(NULL)
is.coprime = function(a, b = NULL, tol = sqrt(.Machine$double.eps), only.answer = TRUE, debug = FALSE) {
  
  # check inputs polm objects a(z), b(z)
  if (!inherits(a, 'ratm')) {
    a = try(polm(a), silent = TRUE)
    if (inherits(a, 'try-error')) {
      stop('could not coerce input "a" to a "polm" object')
    }
  } 
  if ( !( inherits(a,'polm') || inherits(a,'lmfd') || inherits(a, 'rmfd') ) ) {
    stop('input "a" must be a "polm", "lmfd", or "rmfd" object or an object which is coercible to a "polm" object')
  }
  
  # Concatenate/transpose MFDs
  if (inherits(a,'lmfd')) {
    c = cbind(a$a, a$b)
  } else if (inherits(a,'rmfd')){
    c = cbind(t(a$c), t(a$d))
  } else {
    if (is.null(b)) {
      c = a
    } else {
      c = try(cbind(a,b), silent = TRUE)
      if (inherits(c, 'try-error')) {
        stop('inputs "a" and "b" must be compatible "polm" objects')
      }
    }
  }
  
  # skip zero leading coefficients and "unclass"
  c = unclass(prune(c, tol = 0))
  d = dim(c)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  
  if (m*n == 0) {
    # c is an empty polynomial
    stop('the polynomial "c=[a,b]" is empty')
  }
  
  if (p == (-1)) {
    # c is a zero polynomial 
    if (only.answer) return(FALSE)
    return(list(answer = FALSE, A = NULL, B = NULL, zeroes = NA_real_, m = integer(0), n = integer(0)))
  }
  
  if (p == 0) {
    # c is a constant polynomial
    dim(c) = c(m,n)
    if (m > n) {
      answer = FALSE
    } else {
      # check singular values
      d = svd(c, nu = 0, nv = 0)$d
      answer = (d[m] >= tol)
    }
    if (only.answer) return(answer)
    if (answer) {
      zeroes = numeric(0)
    } else {
      zeroes = NA_real_
    }
    return(list(answer = answer, A = NULL, B = NULL, zeroes = zeroes, m = integer(0), n = integer(0)))
  }
  
  # construct Pencil A - Bz corresponding to the polynomial c(z)
  c0 = matrix(c[,,1], nrow = m, ncol = n)
  c = -c[,,-1,drop = FALSE]
  dim(c) = c(m, n*p)
  
  A = bdiag(c0, diag(n*(p-1)))
  B = rbind(c, diag(1, nrow = n*(p-1), ncol = n*p))
  
  if (m > n) {
    # c is a "tall" polynomial
    if (only.answer) return(FALSE)
    return(list(answer = FALSE, A = A, B = B, zeroes = NA_real_, m = m, n = n))
  }
  
  if (debug) {
    message('is.coprime() start:')
    cat('pencil (A - Bz):\n A=\n')
    print(round(A, 4))
    cat('B=\n')
    print(round(B, 4))
  }
  
  # dimensions of pencil (A - Bz)
  m = nrow(A)
  n = ncol(A)
  
  row = 1
  col = 1
  step = 1
  mm = integer(0)
  nn = integer(0)
  while ((row <= m) && (col<= n)) {
    # consider the lower, right block A22 = A[row:n, col:m], B22 = B[row:m, col:n]
    
    # Column Trafo: Separate columns of A22 and B22 ####
    # n1 = dimension of the right kernel of B22
    # make the first n1 columns of B22 equal to zero 
    svd_x = svd(B[row:m, col:n, drop = FALSE], nv = (n-col+1), nu = 0)
    n1 = (n-col+1) - sum(svd_x$d > tol)
    svd_x$v = svd_x$v[,(n - col + 1):1, drop = FALSE]
    A[, col:n] = A[, col:n, drop = FALSE] %*% svd_x$v
    B[, col:n] = B[, col:n, drop = FALSE] %*% svd_x$v
    # impose zeroes 
    if (n1 > 0) {
      B[row:m, col:(col+n1-1)] = 0
    }
    
    # Return early: If right-kernel of B22 is empty, i.e. n1 = 0  ####
    # (A22 - B22z) is a square or tall pencil (where B22 has full column rank) => pencil has zeroes
    if (n1 == 0) {
      if (only.answer) {
        return(FALSE)
      }
      mm = c(mm, m-row+1)
      nn = c(nn, n-col+1)
      
      # ( A22 - B22*z ) is a tall pencil => non trivial left kernel for all z!
      if ((m - row) > (n - col)) {
        return(list(answer = FALSE, A = A, B = B, zeroes = NA_real_, m = mm, n = nn))
      }
      
      # ( A22 - B22*z ) is a square, regular pencil
      zeroes = eigen(solve(B[row:m, col:n, drop = FALSE], A[row:m, col:n, drop = FALSE]), only.values = TRUE)$values
      return(list(answer = FALSE, A = A, B = B, zeroes = zeroes, m = mm, n = nn))
    }
    
    # Row Trafo: Separate rows of A and B ####
    # m1 = rank of A22
    # make the last ((m - row + 1) - m1) rows of A22 equal to zero
    # => the first m1 rows are linearly independent
    # Note: m1 may be zero
    svd_x = svd(A[row:m, col:(col+n1-1), drop = FALSE], nu = (m - row + 1), nv = 0)
    m1 = sum(svd_x$d > tol)
    A[row:m, col:n] = t(svd_x$u) %*% A[row:m, col:n, drop = FALSE]
    B[row:m, col:n] = t(svd_x$u) %*% B[row:m, col:n, drop = FALSE]
    # impose zeroes 
    if (m1 < (m - row +1)) {
      A[(row+m1):m, col:(col+n1-1)] = 0
    }
    
    if (debug) {
      message(paste('is.coprime() step ', step, ' (',
                    row, ':', row + m1 - 1, ' x ',
                    col, ':', col + n1 - 1, '):', sep=''))
      cat('pencil (A - Bz):\n A=\n')
      print(round(A, 4))
      cat('B=\n')
      print(round(B, 4))
    }
    
    row = row + m1
    col = col + n1
    mm = c(mm, m1)
    nn = c(nn, n1)
    
    step = step + 1
  }
  
  if (row > m) {
    # coprime!
    if (only.answer) {
      return(TRUE)
    }
    nn[length(nn)] = nn[length(nn)] + (n - sum(nn))
    return(list(answer = TRUE, A = A, B = B, zeroes = numeric(0), m = mm, n = nn))
  }
  # not coprime 
  if (only.answer) {
    return(FALSE)
  }
  mm[length(mm)] = mm[length(mm)] + (m - sum(mm))
  return(list(answer = FALSE, A = A, B = B, zeroes = NA_real_, m = mm, n = nn))
}

#' Companion Matrix of a Polynomial Matrix
#'
#' Computes a companion matrix for a square (\eqn{m,m)}-dimensional), matrix polynomial 
#' \deqn{a(z) = a_0 + a_1 z + \cdots + a_p z^p}{a(z) = a[0] + a[1] z + \dots + a[p] z^p}
#' The companion matrix is e.g. used to determine the zeroes of a polynomial matrix,
#' see \code{\link{zeroes}}. 
#' \cr
#' Note that the function throws an error, if the constant term \eqn{a_0}{a[0]} is singular.
#' There is no check whether some of the leading coefficients are zero. 
#' So the results is an \eqn{(mp,mp)}-dimensional matrix, even if \eqn{a_p}{a[p]} is zero. 
#'
#' @param a A square polynomial matrix, i.e. an object of class \code{\link{polm}}. 
#'
#' @return A (companion) matrix of dimensions \eqn{(mp,mp)}. 
#' @export
#'
#' @examples
#' companion_matrix(polm(c(1,0,0,0.5,0))) # scalar polynomial
#' companion_matrix(polm(diag(3)))        # zero degree polynomial 
#' companion_matrix(polm(dbind(d = 3, diag(2), -test_array(dim = c(2,2,1)))))
#' companion_matrix(polm(dbind(d = 3, diag(2), -test_array(dim = c(2,2,2)))))
#' 
#' \dontrun{
#' # the following examples throw an error
#' companion_matrix(polm(c(0,0,0,0.5))) # constant term is zero
#' companion_matrix(polm(test_array(dim = c(2,1,3)))) # non-square polynomial
#' companion_matrix(polm(test_array(dim = c(2,2,0)))) # zero polynomial
#' }
companion_matrix = function(a) {
  if (!inherits(a, 'polm')) {
    a = try(polm(a), silent = TRUE)
    if (inherits(a, 'try-error')) stop('argument "a" is not coercible to polm object!')
  }
  a = unclass(a)
  d = dim(a)
  if ((d[1] != d[2]) || (d[3] <= 0)) stop('argument "a" must represent a square, non-singular polynomial matrix')
  
  m = d[1]
  p = d[3] - 1
  
  if (m > 0) {
    # check a(0)
    a0 = try(solve(matrix(a[,,1], nrow = m, ncol = m)), silent = TRUE)
    if (inherits(a0, 'try-error')) stop('constant term a[0] is non invertible')
  }
  
  if ((m*p) == 0) return(matrix(0, nrow = 0, ncol = 0))
  
  # coerce to (m,m(p+1)) matrix
  dim(a) = c(m, m*(p+1))
  
  # normalize constant term a[0] -> I, a[i] -> - a[0]^{-1} a[i]
  a = (-a0) %*% a[, (m+1):(m*(p+1)), drop = FALSE]
  
  if (p == 1) {
    return(a)
  }
  return( rbind(a, diag(x = 1, nrow = m*(p-1), ncol = m*p)) )
}


# Degree related ####

#' Polynomial Degree
#'
#' Compute the polynomial degrees (of the elements) of a polynomial matrix. 
#' Note that for a (scalar) polynomial with zero coefficients the degree 
#' is set to \eqn{-1}.
#' 
#' The main advantage of setting the degree of a zero polynomial to \eqn{-1} rather than \code{-Inf}
#' is that indexing and assigning is programmatically easier: 
#' E.g., \code{p[,,(deg+2):(max_deg+1)] = 0} (in obvious notation) also works when \code{deg = -1}.
#' When multiplying with zero polynomials, it does not make a difference whether one needs to check for \code{deg = -Inf} or \code{deg = -1}.
#' 
#' @param x A polynomial matrix, i.e. an object of class \code{\link{polm}}.
#' @param which (character string) decides whether a matrix with the respectives degrees 
#'              of the entries of the matrix, or a vector with the respective maximal 
#'              degrees in each row or column, or simply the maximum degree of all 
#'              elements of the polynomial matrix is computed. 
#'
#' @return The outcome depends on the parameter \code{which}: 
#' \item{elements}{A matrix with the degrees of the respective elements 
#'                of the polynomial matrix.}
#' \item{rows}{A vector with the maximum degrees within each row.}
#' \item{columns}{A vector with the maximum degrees within each column.}
#' \item{matrix}{maximum of the degrees of the elements of the matrix.}
#'
#' @export
#'
#' @examples
#' x = polm(array(c(0,1,1,0,
#'                  0,0,1,0,
#'                  0,0,0,1,
#'                  0,0,0,0), dim = c(2,2,4)))
#' x
#' degree(x)
#' degree(x, 'rows')
#' degree(x, 'columns')
#' degree(x, 'matrix')
degree = function(x, which = c('elements', 'rows', 'columns', 'matrix')) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be a "polm" object!')
  }
  which = match.arg(which)
  x = unclass(x)
  # degree of a univariate polynomial (= vector): 
  # length of vector - 1 - number of zero leading coefficients
  deg_scalar = function(x) {
    length(x) - sum(cumprod(rev(x) == 0)) - 1
  }
  deg = apply(x, MARGIN = c(1,2), FUN = deg_scalar)
  if (which == 'matrix') return(max(deg))
  if (which == 'columns') return(apply(deg, MARGIN = 2, FUN = max))
  if (which == 'rows') return(apply(deg, MARGIN = 1, FUN = max))
  return(deg)
}

#' Column End Matrix of a Polynomial Matrix
#' 
#' The \emph{column end matrix} of an \eqn{(m,n)}-dimensional polynomial matrix 
#' \eqn{a(z)=a_0 + a_1 z + \cdots + a_p z^p}{a(z)=a[0] + a[1] z + \dots + a[p] z^p} is defined as follows. 
#' Suppose that the maximum degree of the elements in the \eqn{i}-th column is \eqn{p_i}{p[i]}. Then 
#' the column end matrix is the \eqn{(m,n)} matrix with \eqn{i}-th column equal to the 
#' \eqn{i}-th column of the coefficient matrix \eqn{a_{p_i}}{a[p[i]]}. If a column of 
#' \eqn{a(z)} is zero, then
#' the elements of the corresponding column of the column end matrix are set to \code{NA}'s.
#'
#' @param x A polynomial matrix, i.e. an object of class \code{\link{polm}}.
#'
#' @return The column end matrix.
#' @export
#'
#' @examples
#' x = polm(array(c(0,1,1,0,
#'                  0,0,1,0,
#'                  0,0,0,1,
#'                  0,0,0,0), dim = c(2,2,4)))
#' x
#' degree(x)
#' degree(x, 'columns')
#' col_end_matrix(x)
col_end_matrix = function(x) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be a "polm" object!')
  }
  d = dim(x)
  x = unclass(x)
  NAvalue = ifelse(is.complex(x), NA_complex_, NA_real_)
  m = matrix(NAvalue, nrow = d[1], d[2])
  if (length(x) == 0) {
    return(m)
  }

  # degree of a univariate polynomial (= vector): 
  # length of vector - 1 - number of zero leading coefficients
  deg_scalar = function(x) {
    length(x) - sum(cumprod(rev(x) == 0)) - 1
  }
  deg = apply(x, MARGIN = c(1,2), FUN = deg_scalar)
  col_deg = apply(deg, MARGIN = 2, FUN = max)
  
  for (i in iseq(1, dim(x)[2])) {
    if (col_deg[i] >= 0) m[,i] = x[,i,col_deg[i]+1]
  }
  return(m)
}

#' Prune (Laurent) Matrix Polynomial
#'
#' Performs three steps to simplify a matrix polynomial.
#' \enumerate{
#' \item All leading coefficients where the absolute values of the real and the imaginary parts are
#'   less than or equal to \code{tol} are set to zero.
#'   For Laurent polynomials, the same is happening from the other direction.
#' \item The zero leading coefficient matrices are dropped. 
#'   For Laurent polynomials, the same is happening for the coefficient matrices pertaining to low powers.
#' \item If all the absolute values of the imaginary parts of the coefficients are less than or
#'       equal to \code{tol} then the coefficients are set to real values.
#' }
#' Empty polynomial matrices (i.e. matrices with zero rows or columns) are set to polynomial of zero degree.
#'
#' @param x \code{\link{polm}} or \code{\link{lpolm}} object.
#' @param tol Double. Tolerance parameter. Default set to \code{sqrt(.Machine$double.eps)}.
#' @param brutal Boolean. Default set to FALSE.
#'   If TRUE, all small elements are set to zero (irrespective of whether they are leading
#'   coefficients or not).
#'
#' @return A Laurent matrix polynomial, i.e. a \code{\link{polm}} or \code{\link{lpolm}} object.
#' @export
#'
#' @examples
#' x = polm(array(c(1,0,0,0,
#'                  0,1,0,0,
#'                  0,0,1,0,
#'                  0,0,0,1,
#'                  0,0,0,0), dim = c(2,2,5)) + 1e-4)
#' x
#' prune(x, tol = 1e-3)
#' prune(x, tol = 1e-3, brutal = TRUE)
#'
#' # Case of complex variables:
#' x = x + complex(imaginary = 1e-5)
#' x
#' prune(x, tol = 1e-3)
#'
#' # also works for constant matrix polynomials (i.e. matrices)
#' x = polm(array(0:3, dim = c(2,2,1))+1e-4)
#' x
#' prune(x, tol = 1e-3)
#' 
#' # empty polynomials are coerced to polynomials of degree zero 
#' x = polm(array(0, dim = c(0,2,5)))
#' x
#' prune(x)
#' 
#' # Laurent polynomials:
#' (lp = lpolm(c(0, 1:3, 0), min_deg = 2))
#' prune(lp)
prune = function(x, tol = sqrt(.Machine$double.eps), brutal = FALSE) {
  
  x_class = class(x)
  if (!inherits(x, c('polm', 'lpolm'))) {
    stop('argument "x" must be a polm or lpolm object!')
  }
  if ("lpolm" %in% x_class){
    min_deg = attr(x, "min_deg")
  }
  x = unclass(x)
  d = dim(x)
  if (min(d) <= 0) {
    # empty polynomial, or polynomial of degree (-1)
    return(polm(array(0, dim = c(d[1], d[2], 0))))
  }

  # step one: Set all small leading coefficients to zero
  issmall = ( (abs(Re(x)) <= tol) & (abs(Im(x)) <= tol) )
  
  issmall_polm = apply(issmall, MARGIN = c(1,2), FUN = function(x) { rev(cumprod(rev(x))) }) 
  
  # apply: returns an array of dim (d[3], d[1], d[2]) if d[3] > 1
  
  # make sure that issmall_polm is an array (also in the case where the matrix polynomial is constant)
  dim(issmall_polm) = d[c(3,1,2)]
  
  # permute the dimensions back to the polm form: 
  # necessary because apply returns an array of dim (d[3], d[1], d[2]) if d[3] > 1
  issmall_polm = aperm(issmall_polm, c(2,3,1))
  issmall_polm = (issmall_polm == 1)
  
  # finish step one
  x[issmall_polm] = 0
  
  # Same steps in the other direction for lpolm
  if ("lpolm" %in% x_class){
    issmall_lpolm = apply(issmall, MARGIN = c(1,2), FUN = function(x) { cumprod(x) }) 
    dim(issmall_lpolm) = d[c(3,1,2)]
    issmall_lpolm = aperm(issmall_lpolm, c(2,3,1))
    issmall_lpolm = (issmall_lpolm == 1)
    x[issmall_lpolm] = 0
  }
  
  # step two: drop leading zero matrix coefficients
  keep = apply(!issmall_polm, MARGIN = 3, FUN = any)
  
  if ("polm" %in% x_class){
    # keep[1] = TRUE # keep the constant
    keep
    x = x[,, keep, drop = FALSE]
  } else if("lpolm" %in% x_class){
    keep_lpolm = apply(!issmall_lpolm, MARGIN = 3, FUN = any)
    keep_lpolm = keep_lpolm * keep
    keep_lpolm = (keep_lpolm == 1)
    x = x[,, keep_lpolm, drop = FALSE]
    # Adjust min_deg
    min_deg = min_deg + sum(cumprod(!keep_lpolm))
  }
  
  # step three: drop imaginary part if all imaginary parts are small
  if (is.complex(x)) {
    if ( all(abs(Im(x)) <= tol) ) {
      x = Re(x)
    }
  }
  
  # This option is provided to see, e.g., the lower triangularity of the zero power coefficient 
  # matrix when using "transform_lower_triangular"
  if (brutal){
    issmall_brutal = ( (abs(Re(x)) <= tol) & (abs(Im(x)) <= tol) )
    x[issmall_brutal] = 0
  }
  
  if ("polm" %in% x_class){
    x = polm(x)
  } else if ("lpolm" %in% x_class){
    x = lpolm(x, min_deg = min_deg)
  }

  return(x)
}

# Small internal helpers for column reduction ####

# internal function
# l2 norm of a vector
l2_norm = function(x){
  return(sqrt(sum(x^2)))
} 

# internal function
# consider a vector x = c(x[1], ..., x[k],  x[k+1], ..., x[n])
# return a logical  i = c(FALSE, .., FALSE, TRUE, .., TRuE) 
# where k is the minimum integer such that | x[s] | <= tol for all s > k 
is_small = function(x, tol = sqrt(.Machine$double.eps), count = TRUE) {
  if (length(x) == 0) {
    i = logical(0)
  } else {
    i = rev( cumprod( rev(abs(x) <= tol) ) == 1 )
  }
  if (count) {
    return(sum(i))
  } else {
    return(i)
  }
}

# internal function
# consider a vector x = c(x[1], ..., x[k],  x[k+1], ..., x[n])
# return a logical  i = c(TRUE, .., TRUE, FALSE, .., FALSE) 
# where k is the minimum integer such that | x[s] | <= tol for all s > k 
is_large = function(x, tol = sqrt(.Machine$double.eps), count = TRUE) {
  if (length(x) == 0) {
    i = logical(0)
  } else {
    i = rev( cumprod( rev(abs(x) <= tol) ) == 0 )
  }
  if (count) {
    return(sum(i))
  } else {
    return(i)
  }
}

# Normal forms (Hermite, Smith, Wiener-Hopf) and essential helpers ####

#' Purge Rows or Columns of a Polynomial Matrix
#'
#' This helper function is the main work horse for computing the Hermite normal form 
#' (see \code{\link{hnf}}) and the Smith normal form (see \code{\link{snf}}) of 
#' polynomial matrices. It "purges" all elements below, above, to the right or to the 
#' left of a pivot element by elementary row- or column- operations. Here 
#' "purge" means that the elements are either reduced to zero or that the degree of the 
#' elements is made smaller than the degree of the pivot element. 
#' 
#' Suppose that the matrix \eqn{a(z)} has \eqn{m} rows, \eqn{n} columns and that the pivot is at 
#' position \eqn{(i,j)}. Furthermore, let us first consider the case \code{direction='down'} 
#' and \code{permute=FALSE}. In this case a suitable multiple - which is computed by the 
#' Euclidean polynomial division algorithm - of the \eqn{i}-th row is subtracted from all rows 
#' below the \eqn{i}-th row such that the respective degree of the elements below 
#' the pivot element have a degree which is smaller than the degree of the pivot. 
#' 
#' If the option \code{permute=TRUE} then first the rows \eqn{i:m} are permuted such that 
#' the \eqn{(i,j)}-th element has the smallest degree among all elements in the \eqn{j}-th 
#' column and the rows \eqn{i:m}. Next a suitable multiple of the \eqn{i}-th row is 
#' subtracted from all rows below the \eqn{i}-th row such that the respective degree 
#' of the elements below the pivot element have a degree which is smaller than the 
#' degree of the pivot. These two steps are repeated until all elements below the pivot element are zero. 
#' 
#' Quite analogously the cases \code{direction='up'},  \code{direction='left'} and  \code{direction='right'}
#' may be discussed. Note however, that for the cases  \code{direction='left'} and  \code{direction='right'} 
#' elementary column-operations are used. 
#' 
#' Finally, for \code{monic=TRUE} the pivot element is made monic, by multiplying the 
#' respective row (or column) by a suitable scalar. 
#'
#' @param a Polynomial matrix of dimension \eqn{(m,n)}, i.e. an object of class \code{\link{polm}}.
#' @param pivot Integer vector of length 2. Specifies the position of the "pivot" element.
#' @param direction Character string. 
#'   \describe{
#'   \item{down}{(default) "purge" all elements below the pivot element by elementary row-operations.}
#'   \item{up}{"purge" all elements above the pivot element by elementary row-operations}
#'   \item{right}{"purge" all elements to the right of the pivot element 
#'               by elementary column-operations}
#'   \item{left}{"purge" all elements to the left of the pivot element 
#'               by elementary column-operations}
#'   }
#' @param permute Logical, defaults to TRUE. See the details below.
#' @param tol Tolerance parameter, used for "pruning" the polynomial matrix (after each step). 
#'   See \code{\link{prune}}.
#' @param monic Logical, defaults to FALSE.
#'   If TRUE, the coefficient pertaining to the highest degree of the pivot element will 
#'   be normalized to 1.
#' @param debug Logical, default to FALSE. If TRUE, then some diagnostic messages are printed.
#'
#' @return List with three slots
#'   \item{\code{h}}{Polynomial matrix of dimension \eqn{(m,n)}, i.e. an object of class \code{\link{polm}}. 
#'         This matrix is the result of the "purging" operation(s).}
#'   \item{\code{u}}{Unimodular polynomial matrix, i.e. a class \code{\link{polm}} object. 
#'        For \code{direction = 'down'} or  \code{direction = 'up'} the matrix 
#'        \eqn{u(z)} is \eqn{(m,m)} dimensional and satisfies \eqn{a(z) = u(z) h(z)}.
#'        For \code{direction = 'left'} or  \code{direction = 'right'} the matrix 
#'        \eqn{u(z)} is \eqn{(n,n)} dimensional and satisfies \eqn{a(z) = h(z) u(z)}.}
#'   \item{\code{u_inv}}{Unimodular polynomial matrix, i.e. a class \code{\link{polm}} object.  
#'        This matrix is the inverse of \eqn{u(z)}.}
#'   
#' @export
#'
#' @examples
#' # Generate matrix polynomial
#' a = test_polm(dim = c(2,3), degree = 1)
#' print(a, format = 'c')
#'
#' ########################################################
#' # Purge first column downwards
#' out = purge_rc(a, pivot = c(1,1), monic = TRUE)
#'
#' # First col zero except for (1,1) element
#' print(out$h, digits = 2, format = 'c')
#'
#' # Check polynomial matrix products
#' all.equal(a, prune(out$u %r% out$h))
#' all.equal(out$h, prune(out$u_inv %r% a))
#' all.equal(polm(diag(2)), prune(out$u_inv %r% out$u))
#'
#' ########################################################
#' # Purge last column upwards
#' out = purge_rc(a, pivot = c(2,3), direction = "up", monic = TRUE)
#'
#' # Last col zero except for (2,3) element
#' print(out$h, digits = 2, format = 'c')
#'
#' # Check polynomial matrix products
#' all.equal(a, prune(out$u %r% out$h))
#' all.equal(out$h, prune(out$u_inv %r% a))
#' all.equal(polm(diag(2)), prune(out$u_inv %r% out$u))
#' 
#' ########################################################
#' # Purge first row right
#' out = purge_rc(a, pivot = c(1,1), direction = "right", monic = TRUE)
#'
#' # first row zero except for (1,1) element
#' print(out$h, digits = 2, format = 'c')
#'
#' # Check polynomial matrix products
#' all.equal(a, prune(out$h %r% out$u))
#' all.equal(out$h, prune(a %r% out$u_inv))
#' all.equal(polm(diag(3)), prune(out$u_inv %r% out$u))
#' 
#' ########################################################
#' # Purge last row left
#' out = purge_rc(a, pivot = c(2,3), direction = "left", monic = TRUE)
#'
#' # last row zero except for (2,3) element
#' print(out$h, digits = 2, format = 'c')
#'
#' # Check polynomial matrix products
#' all.equal(a, prune(out$h %r% out$u))
#' all.equal(out$h, prune(a %r% out$u_inv))
#' all.equal(polm(diag(3)), prune(out$u_inv %r% out$u))
purge_rc = function(a, pivot = c(1,1), direction = c('down','up','left','right'), 
                    permute = TRUE, tol = sqrt(.Machine$double.eps), 
                    monic = FALSE, debug = FALSE) {
  direction = match.arg(direction)

  # Check argument 'a'
  stopifnot("purge_rc(): Argument *a* is not a polm object!" = inherits(a, 'polm'))
  
  # Dimensions of input
  d = dim(unclass(a))
  m = d[1]
  n = d[2]
  p0 = d[3] - 1
  stopifnot("purge_rc(): Argument *a* must have more than zero inputs and outputs!" = m*n != 0)
  
  # check pivot
  pivot = as.integer(as.vector(pivot)) 
  stopifnot("purge_rc(): Argument *pivot* must be an integer vector of length 2, 1 <= pivot <= dim(a)" = (length(pivot) == 2) && (min(pivot) > 0) && (pivot[1] <= m) && (pivot[2] <= n))
  
  i = pivot[1]
  j = pivot[2]
  
  # If direction is not "down", transform (transpose etc) the polm object
  if (direction == 'up') {
    a = a[m:1, ]
    i = m - (i - 1)
  }
  if (direction == 'right') {
    a = t(a)
    
    junk = i
    i = j
    j = junk
    
    junk = m
    m = n
    n = junk
  }
  if (direction == 'left') {
    a = t(a)
    
    junk = i
    i = j
    j = junk
    
    junk = m
    m = n
    n = junk
    
    a = a[m:1, ]
    i = m - (i - 1)
  }
  
  # Initialization of unimodular matrix 
  u0 = polm(diag(m))
  u = u0
  u_inv = u0
  
  # (m x n) matrix of degrees of each element
  p = degree(a)  
  
  # degrees of entries in the j-th column
  p_col = p[, j]
  
  # no permutations allowed, but pivot element is zero!
  if ( (i < m) && (!permute) && (p_col[i] == -1) && any(p_col[(i+1):m] > -1) ) {
    stop("purge_rc(): Pivot element is zero but permutation is not allowed. Purging not possible.") 
  }
  
  # Main iteration ####
  iteration = 0
  
  # The column is not purged if 
  # (in the case where permutations are allowed) any element below the pivot is non-zero
  # (in the case where permutations are not allowed) any element below the is non-zero is of equal or larger degree than the pivot
  not_purged = (i < m) && ( ( permute && any(p_col[(i+1):m] > -1) ) || ( (!permute) && any(p_col[(i+1):m] >= p_col[i]) ) )

  while ( not_purged )  {
    
    iteration = iteration + 1

    if (debug) {
      message('purge_rc: iteration=', iteration)
      print(a, format = 'i|zj', digits = 2)
      # print(a)
      print(p)
    }
    
    if (permute){ 
      # Permutation step 

      # find (non zero) entry with smallest degree
      p_col[iseq(1, i-1)] = Inf   # ignore entries above the i-th row 
      p_col[p_col == -1]  = Inf   # ignore zero entries 
      k = which.min(p_col)

      # permute i-th row and k-th row
      perm = 1:m
      perm[c(k,i)] = c(i,k)
      
      a = a[perm, ]
      u = u[, perm]
      u_inv = u_inv[perm, ]
      p_col = p_col[perm]
    }
    
    # Division step 

    q = a[(i+1):m, j] %/% a[i, j] # polynomial divison
      
    M = u0
    Mi = u0
    
    M[(i+1):m, i] = -q
    Mi[(i+1):m, i] = q
    
    a = prune(M %r% a, tol = tol)
    u = u %r% Mi
    u_inv = M %r% u_inv
    
    p = degree(a)
    
    # degrees of entries in the j-th column
    p_col = p[, j]
    stopifnot("purge_rc(): Reduction of degree failed! This should not happen." = all(p_col[(i+1):m] < p_col[i]))
    
    # The column is not purged if 
    # (in the case where permutations are allowed) any element below the pivot is non-zero
    # (in the case where permutations are not allowed) any element below the is non-zero is of equal or larger degree than the pivot
    not_purged = (i < m) && ( ( permute && any(p_col[(i+1):m] > -1) ) || ( (!permute) && any(p_col[(i+1):m] >= p_col[i]) ) )
  }
  
  if (monic) {
    if (p[i,j] >= 0) {
      c = unclass(a)[i,j,p[i,j]+1]
      a[i, ] = a[i, ] %/% c # polynomial division, see ?Ops.ratm
      u[, i] = u[, i] * c
      u_inv[i, ] = u_inv[i, ] %/% c
    }
  }
  
  # If direction is not "down", transform (transpose etc) the polm object back
  if (direction == 'up') {
    a = a[m:1, ]
    u = u[m:1, m:1]
    u_inv = u_inv[m:1, m:1]
  }
  if (direction == 'right') {
    a = t(a)
    u = t(u)
    u_inv = t(u_inv)
  }
  if (direction == 'left') {
    a = t(a[m:1,])
    u = t(u[m:1,m:1])
    u_inv = t(u_inv[m:1,m:1])
  }
  
  return(list(h = a, u = u, u_inv = u_inv))
}

#' Construct a Column Reduced Polynomial Matrix
#' 
#' Let \eqn{a(z)} be a square (non singular) polynomial matrix. 
#' This helper function constructs a unimodular transformation matrix \eqn{v(z)} such that \eqn{a(z) v^{-1}(z)} is column reduced (i.e. the column end matrix has full rank).
#' Algorithmic implementation are described e.g. in \insertCite{Wolovich1974}{rationalmatrices} Theorem 2.5.7, page 27, \insertCite{KrishnaraoChen84_colred}{rationalmatrices}, and \insertCite{geurtspraagman96}{rationalmatrices} who show that the KC implementations fails when the degree of the unimodular matrix \eqn{v^{-1}(z)} exceeds the degree of \eqn{a(z)} (page 4 in GP).
#' While all these implementations use elementary column operations to obtain zero columns in the column-end-matrix (in order to reduce the degree of the matrix polynomial), this implementation uses the \link{svd}.
#' The examples below are taken from \insertCite{geurtspraagman96}{rationalmatrices} and \insertCite{KrishnaraoChen84_colred}{rationalmatrices}.
#' 
#' @section Possible "Improvements":
#' It is not clear whether the changes are improvements...
#' First: When a rank deficiency in the column-end-matrix is detected with the SVD, only one column is "reduced" to zero (even when the rank deficiency is larger than one). 
#' Fewer SVDs are calculated when both 
#' It might be better to reduce all columns to zero which pertain to (numerically) zero eigenvalues.
#' \cr
#' Second: The pivoting mechanism in the QR decomposition might be useful to single out the columns which should be set to zero.
#' It might be preferable to the SVD.
#' 
#' 
#' @param a \code{\link{polm}} object, which represents the square, polynomial matrix 
#'          \eqn{a(z)}.
#' @param tol Double. Tolerance parameter. Default set to sqrt(.Machine$double.eps).
#' @param debug Logical, default to FALSE. If TRUE, then some diagnostic messages are printed.
#'
#' @return List with components 
#' \item{a}{(\code{\link{polm}} object) is the transformed (column reduced) matrix.}
#' \item{v,v_inv}{(\code{\link{polm}} objects) are the the unimodular matrices \eqn{v(z)} and \eqn{v^{-1}(z)} such that \eqn{a(z) v^{-1}(z)} is column reduced}.
#' \item{col_degrees}{vector of column degrees of the transformed matrix. 
#'     Note that the columns are permuted such that the transformed matrix has *non-increasing* column degrees.}
#' \item{col_end_matrix}{the column end matrix of the transformed matrix.}
#' 
#' @export
#' 
#' @references    
#' \insertRef{Wolovich1974}{rationalmatrices}
#' \insertRef{KrishnaraoChen84_colred}{rationalmatrices}
#' \insertRef{geurtspraagman96}{rationalmatrices}
#' 
#' @seealso The column end matrix may be computed with \code{\link{col_end_matrix}}. 
#' The function \code{col_reduce} is mainly used to compute the Wiener-Hopf factorization of a polynomial matrix, see \code{\link{whf}}.
#'
#' @examples
#' # #############################################################################
#' # define a simple utility function for the computation of the rank of a matrix 
#' # compare e.g. Matrix::rankMatrix
#' rkM = function(x) {
#'   m = nrow(x)
#'   n = ncol(x)
#'   tol = max(m,n) * .Machine$double.eps
#'   
#'   if (min(m,n) == 0) return(0L)
#'   s = svd(x, 0, 0)$d
#'   return(sum(s >= tol*max(s)))
#' }
#' # #############################################################################
#' 
#' z = polm(c(0,1))
#' 
#' # Example 2.5.4 in W. A. Wolovich, Linear Multivariable Systems ###############
#' a = matrix(c(-3,2,0,1,2,3,0,0,2), nrow = 3, ncol = 3) +
#'     matrix(c(0,4,0,0,0,1,2,0,-3), nrow = 3, ncol = 3) * z + 
#'     matrix(c(1,0,-1,0,0,0,0,0,0), nrow = 3, ncol = 3) * z^2 
#' # Original Matrix:
#' print(a, format = 'c')
#' # Its column end matrix (and its rank) and column degrees
#' col_end_matrix(a)
#' col_end_matrix(a) %>% svd() %>% .$d
#' degree(a, "c")
#' 
#' # After column reduction:
#' out = col_reduce(a)
#' print(out$a, format = 'c', digits = 2)
#' print(out$col_degrees)
#' print(out$col_end_matrix)
#' print(out$col_end_matrix %>% svd() %>% .$d)
#' 
#' # Check correctness:
#' all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#' all.equal(prune(a), prune(out$a %r% out$v))
#'
#' 
#'   
#' # Random example: col degrees = (0,1,-1): throws an error #####################################
#' a = test_polm(dim = c(3,3), degree = c(0,1,-1), random = TRUE, digits = 2)
#' print(a, format = 'c')
#' 
#' \dontrun{
#' # this throws an error, since a(z) is singular
#' out = col_reduce(a)
#' }
#' 
#' 
#' # Random example: Generic matrices are row reduced ############################################
#' a = test_polm(dim = c(3,3), degree = c(2,1,0), random = TRUE, digits = 2)
#' print(a, format = 'c')
#' 
#' # Column reduction:
#' out = col_reduce(a) 
#' print(out$a, format = 'c', digits = 2)
#' print(out$col_degrees)
#' print(out$col_end_matrix)
#' print(out$col_end_matrix %>% svd() %>% .$d)
#' 
#' # Check:
#' all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#' all.equal(prune(a), prune(out$a %r% out$v))
#' 
#' 
#' # Random example: Column end matrix has rank 2 ################################
#' col_end_matrix = matrix(round(rnorm(2*3),1), nrow = 3, ncol = 2) %*% 
#'                  matrix(round(rnorm(2*3),1), nrow = 2, ncol = 3)
#' a = test_polm(dim = c(3,3), degree = c(2,1,0), random = TRUE, 
#'                digits = 2, col_end_matrix = col_end_matrix)
#' print(a, format = 'c')
#' a %>% degree("c")
#' print(a %>% col_end_matrix)
#' print(a %>% col_end_matrix %>% svd() %>% .$d)
#' 
#' # Column reduction:
#' out = col_reduce(a) 
#' print(out$a, format = 'c', digits = 2)
#' print(out$col_degrees)
#' print(out$col_end_matrix)
#' print(out$col_end_matrix %>% svd() %>% .$d)
#' 
#' # Check:
#' all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#' all.equal(prune(a), prune(out$a %r% out$v))
#' 
#' 
#' # Random example: Column end matrix has rank 1 ################################ 
#' col_end_matrix = matrix(round(rnorm(3),1), nrow = 3, ncol = 1) %*% 
#'                  matrix(round(rnorm(3),1), nrow = 1, ncol = 3)
#' a = test_polm(dim = c(3,3), degree = c(2,1,1), random = TRUE, 
#'                digits = 2, col_end_matrix = col_end_matrix)
#' print(a, format = 'c')
#' a %>% degree("c")
#' a %>% col_end_matrix() %>% svd() %>% .$d
#' 
#' # Column reduction:
#' out = col_reduce(a, debug = FALSE) 
#' print(out$a, format = 'c', digits = 2)
#' print(out$col_degrees)
#' print(out$col_end_matrix)
#' print(out$col_end_matrix %>% svd() %>% .$d)
#' 
#' # Check:
#' all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#' all.equal(prune(a), prune(out$a %r% out$v))
#' 
#' #################################################################################
#' 
#' # PG, Ex 1: Jordan Normal Form type (Unimodular matrix) ##########
#' # Result: Works fine here (does not work with KC implementation)
#' m0 = diag(3)
#' m1 = matrix(c(0,1,0,   
#'               0,0,1,   
#'               0,0,0), nrow = 3, byrow = TRUE)
#' 
#' # Polymat:
#' (a = polm(array(c(m0,m1), dim = c(3,3,2))))
#' a %>% print(format = "c")
#' a %>% degree("c")
#' 
#' # Column reduced
#' a_red = col_reduce(a)
#' a_red$a
#' degree(a_red$a, "c")
#' 
#' #################################################################################
#' 
#' # PG, Ex 2: Nothing special ##########
#' # Result: Same as PG
#' 
#' a = test_polm(dim = c(2,2), degree = -1)
#' a[1,1] = polm(c(4,12,13,6,1))
#' a[1,2] = polm(c(-2,-5,-4,-1))
#' a[2,2] = polm(c(2,1))
#' a
#' a %>% print(format = "c")
#' col_end_matrix(a)
#' degree(a, "c")
#' 
#' # Column reduction:
#' a_red = col_reduce(a)
#' a_red$a
#' a_red$a %>% print(format = "c")
#' degree(a_red$a, "c")
#' 
#' #################################################################################
#' 
#' # PG, Ex 3: Unimodular matrix ##########
#' # Result: Works, but different matrices, v_inv (U(s) in PG notation is slightly "smaller" here)
#' 
#' a = test_polm(dim = c(3,3), degree = -1)
#' a[1,1] = polm(c(0,0,0,0,1))
#' a[1,2] = polm(c(0,0,1))
#' a[1,3] = polm(c(1,0,0,0,0,0,1))
#' a[2,1] = polm(c(0,0,1))
#' a[2,2] = polm(1)
#' a[2,3] = polm(c(0,0,0,0,1))
#' a[3,1] = polm(1)
#' a[3,3] = polm(1)
#' a
#' a %>% print(format = "c")
#' col_end_matrix(a)
#' degree(a, "c")
#' 
#' # Column reduction:
#' a_red = col_reduce(a)
#' a_red$a
#' a_red$a %>% print(format = "c")
#' degree(a_red$a, "c")
#' 
#' # Verify the column-reduced matrix time v(z) is equal to original one:
#' with(a_red, a %r% v) %>% print(format = "c")
#' a %>% print(format = "c")
#' 
#' # Unimodular matrix transforming the column-reduced matrix to original one:
#' a_red$v %>% print(format = "c")
#' 
#' # Verify that original matrix a(z) times v^{-1}(z) is column-reduced: PRUNING NECESSARY!
#' (a %r% a_red$v_inv) %>% print(format = "c")
#' prune(a %r% a_red$v_inv) %>% print(format = "c")
#' a_red$v_inv %>% print(format = "c") 
#' 
#' #################################################################################
#' 
#' # PG, Ex 4: PG's "singular case" ##########
#' # Result: Different results!
#' #   If eps = 10^(-9) is chosen, 
#' #   then this algorithm breaks down because the a[0] is recognized as singular! 
#' #   If eps = 10^(-4) is chosen (as is done in PG), 
#' #   we obtain a result which is different from the one in PG
#' 
#' 
#' # Original matrix (could be argued to be numerically singular at a[0], depending on tolerance!)
#' a = test_polm(dim = c(3,3), degree = -1)
#' eps = 10^(-4)
#' a[1,1] = polm(c(0,0,1,1))
#' a[1,2] = polm(c(1,eps))
#' a[1,3] = polm(1)
#' a[2,1] = polm(c(0,0,2))
#' a[2,2] = polm(-1)
#' a[2,3] = polm(-1)
#' a[3,1] = polm(c(0,0,3))
#' a[3,2] = polm(1)
#' a[3,3] = polm(1)
#' a
#' a %>% print(format = "c")
#' col_end_matrix(a)
#' col_end_matrix(a) %>% svd() %>% .$d
#' degree(a, "c")
#' 
#' # Column reduction: Note that there is a column with "small length". 
#' #   It depends on the tolerance whether this column is considered to be zero.
#' #   Also, note the singular values of the column-end-matrix of the reduced polymat!     
#' a_red = col_reduce(a)
#' a_red$a
#' a_red$a %>% print(format = "c")
#' a_red$a %>% print(format = "c", digits = 3)
#' 
#' a_red$a %>% col_end_matrix()
#' a_red$a %>% col_end_matrix() %>% svd() %>% .$d
#' a_red$a %>% degree("c")
#' 
#' # Col-reduced matrix time v(z) = original: 
#' # It works up to a small numerical issue in the (1,1) element 
#' with(a_red, a %r% v) %>% print(format = "c")
#' a %>% print(format = "c")
#' 
#' # Check unimodular matrix taking the col-reduced matrix back to original:
#' # Small (i.e. unproblematic) numerical issue in the (2,1)-element) which is also reflected above.
#' a_red$v %>% print(format = "c")
#'
#' # Original times v^{-1}(z) = col-reduced:
#' # Works fine, up to small (non-problematic) issue in (1,3)-element
#' a_red$a %>% print(format = "c")
#' (a %r% a_red$v_inv) %>% print(format = "c")
#' 
#' # Matrix transforming the original to column-reduced:
#' # Small (unproblematic) issue in (2,3)-element
#' a_red$v_inv %>% print(format = "c")
#' 
#' # Algebraic result given in PG: Different from the one obtained here! ####
#' # The element (3,1) of matrix R(s) in PG seems to be incorrect. 
#' # Changing this element below results in a column-reduced matrix.
#' eta = 1/eps
#' 
#' r = test_polm(dim = c(3,3), degree = -1)
#' r[1,1] = polm(c(0,-3*eta))
#' r[1,2] = polm(c(1))
#' r[1,3] = polm(1)
#' r[2,1] = polm(c(0,3*eta))
#' r[2,2] = polm(c(-1, eps))
#' r[2,3] = polm(-1)
#' r[3,1] = polm(c(1,-3*eta,5))
#' r[3,2] = polm(c(1,-eps))
#' r[3,3] = polm(1)
#' r
#' # PG: Col-reduced
#' r %>% print(format = "c")
#' # Here: Col-reduced
#' a_red$a %>% print(format = "c")
#' 
#' u = test_polm(dim = c(3,3), degree = -1)
#' u[1,1] = polm(c(1))
#' u[1,2] = polm(c(0))
#' u[1,3] = polm(0)
#' u[2,1] = polm(c(0,-3*eta,-eta))
#' u[2,2] = polm(c(1))
#' u[2,3] = polm(0)
#' u[3,1] = polm(c(0,0,eta+2))
#' u[3,2] = polm(c(0,-eps))
#' u[3,3] = polm(1)
#' u
#' u %>% print(format = "c")
#' 
#' # Check whether their result makes sense:
#' # Column-reduced matrix r is different from a %r% u!!! 
#' # This mistake can be corrected by adjusting element (3,1) of r
#' a %>% print(format = "c")
#' 
#' a %r% u
#' (a %r% u) %>% print(format = "c")
#' (a %r% u) %>% col_end_matrix()
#' (a %r% u) %>% col_end_matrix() %>% svd() %>% .$d
#' (a %r% u) %>% degree("c")
#' 
#' r
#' r %>% print(format = "c")
#' 
#' #########################
#' # Change an element in r
#' r2 = r
#' r2[3,1] = polm(c(0,-3*eta,5))
#' r2 %>% print(format = "c")
#' r2 %>% col_end_matrix()
#' r2 %>% col_end_matrix() %>% svd() %>% .$d
#' r2 %>% col_end_matrix() %>% rkM()
#' r2 %>% degree("c")
#' 
#' (a %r% u) %>% print(format = "c")
#'  
#' 
#' #################################################################################
#' 
#' # PG, Ex 5: Breaks down in PG, but they give an algebraic solution ##########
#' # Result: Works here! 
#' #  Different result as algebraic solution indicated in PG (up to column permutation).
#' #  Of course, the obtained result is also column-reduced!
#' #  The unimodular matrix (which column-reduces the original polynomial matrix) is different 
#' 
#' a = test_polm(dim = c(4,4), degree = -1)
#' eps = 10^(-9)
#' a[1,1] = polm(c(1,1,1))
#' a[1,2] = polm(c(0,eps))
#' a[1,3] = polm(c(0,0,0,1))
#' a[1,4] = polm(c(1,0,0,1))
#' a[2,1] = polm(c(0,2))
#' a[2,2] = polm(0)
#' a[2,3] = polm(c(1,0,0,2))
#' a[2,4] = a[2,3]
#' a[3,1] = polm(c(1,3))
#' a[3,2] = polm(3)
#' a[3,3] = polm(c(3,0,3))
#' a[3,4] = a[3,3]
#' a[4,1] = polm(c(0,4))
#' a[4,2] = polm(c(0))
#' a[4,3] = polm(c(1,0,0,4))
#' a[4,4] = a[4,3]
#' a
#' 
#' # Original matrix, its column-end-matrix with its singular values, and its column degrees
#' a %>% print(format = "c")
#' col_end_matrix(a)
#' col_end_matrix(a) %>% svd() %>% .$d
#' col_end_matrix(a) %>% rkM()
#' degree(a, "c")
#' 
#' # Column-reduction:
#' a_red = col_reduce(a)
#' a_red$a
#' a_red$a %>% print(format = "c")
#' col_end_matrix(a_red$a)
#' col_end_matrix(a_red$a) %>% svd() %>% .$d
#' degree(a_red$a, "c")
#' 
#' # Check whether col-reduced matrix times unimodular v(z) = original:
#' # Works with small (non-problematic) numerical mistakes.
#' # Brutal pruning "solves" it
#' # Similar for the unimodular matrix v(z) itself
#' with(a_red, a %r% v) %>% print(format = "c")
#' with(a_red, a %r% v) %>% prune(brutal = TRUE) %>% print(format = "c")
#' a_red$v %>% print(format = "c")
#' a_red$v %>% prune(brutal = TRUE) %>% print(format = "c")
#' 
#' # Check whether original times v^{-1}(z) = col-reduced:
#' # Same as above: non-problematic numerical mistakes, "solved" by brutally pruning
#' # Same for unimodular v^{-1}(z)
#' (a %r% a_red$v_inv) %>% print(format = "c")
#' (a %r% a_red$v_inv) %>% prune(brutal = TRUE) %>% print(format = "c")
#' a_red$a %>% print(format = "c")
#' a_red$v_inv %>% print(format = "c")
#' a_red$v_inv %>% prune(brutal = TRUE) %>% print(format = "c")
#' 
#' # Krishnarao and Chen example #####################
#' 
#' (a = polm(array(c(4,0,-2,2,  
#'                   12, 0, -5, 1,    
#'                   13,0,1,0,   
#'                   2,0,0,0), dim = c(2,2,4))))
#' out = col_reduce(a)
#' out$a
#' a %r% out$v_inv
#' with(out, a %r% v)
#' 
#' # Majid: Ex 3 ################
#' 
#' (m0 = matrix(c(-3,1,0, 2,2,0,   0,3,2), nrow = 3, byrow = TRUE))
#' (m1 = matrix(c(0,0,2,   4,0,0,   0,1,-3), nrow = 3, byrow = TRUE))
#' (m2 = matrix(c(1,0,0,   0,0,0,   -1,0,0), nrow = 3, byrow = TRUE))
#' m = polm(array(c(m0,m1,m2), dim = c(3,3,3)))
#' m %>% print(format = "c")
#' m %>% col_end_matrix()
#' m %>% col_end_matrix() %>% svd() %>% .$d
#' m %>% degree("c")
#' 
#' (out = col_reduce(m))
#' 
#'  # Majid: Ex 4 ################3
#' (m0 = matrix(c(1,0,   0,1), nrow = 2, byrow = TRUE))
#' (m1 = matrix(c(0,0,   2,0), nrow = 2, byrow = TRUE))
#' (m2 = matrix(c(1,1,   0,0), nrow = 2, byrow = TRUE))
#' (m3 = matrix(c(2,0,   0,0), nrow = 2, byrow = TRUE))
#' m = polm(array(c(m0,m1,m2,m3), dim = c(2,2,4)))
#' m %>% print(format = "c")
#' m %>% col_end_matrix()
#' m %>% col_end_matrix() %>% svd() %>% .$d
#' m %>% degree("c")
#' 
#' (out = col_reduce(m))
col_reduce = function(a, tol = sqrt(.Machine$double.eps), debug = FALSE) {
  
  # Check inputs
  stopifnot("col_reduce(): Input *a* must be a polm object!" = inherits(a, 'polm'))
  
  # Integer-valued parameters
  x = unclass(a)
  d = dim(x)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  stopifnot("col_reduce(): Input *a* must be a square, non empty, non zero polynomial matrix." = (m*n != 0) && (m == n) && (p >= 0))
  
  # Initialize unimodular matrices 
  v0 = polm(diag(n))
  v = v0
  v_inv = v0
  
  {# # balance the 'column norms'
  # col_norm = apply(x, MARGIN = 2, FUN = l2_norm)
  # a = a %r% diag(sqrt(mean(col_norm^2))/col_norm)
  # v_inv = v_inv %r% diag(sqrt(mean(col_norm^2))/col_norm)
  # v = diag(col_norm / sqrt(mean(col_norm^2))) %r% v
  }
  # Column degrees (taking rounding error into account)
  # output is a matrix! rows correspond to the norm of the respective column, columns correspond to degrees!
  col_norms = apply(x, MARGIN = c(2,3), FUN = l2_norm) 
  col_degrees = apply(col_norms, MARGIN = 1, FUN = is_large, tol = tol, count = TRUE) - 1
  stopifnot("col_reduce(): The input *a* has (close to) zero columns." = min(col_degrees) >= 0)
  
  # set "small columns" to zero and retrieve column end matrix
  col_end_matrix = matrix(0, nrow = m, ncol = n)
  for (i in (1:n)) {
    x[, i, iseq(col_degrees[i]+2, p+1)] = 0
    col_end_matrix[,i] = x[, i, col_degrees[i]+1]
  }
  
  # reduce order of polynomial
  p = max(col_degrees)
  x = x[, , 1:(p+1), drop = FALSE]
  a = polm(x)
  
  # sort by column degrees 
  o = order(col_degrees, decreasing = FALSE)
  col_degrees = col_degrees[o]
  col_end_matrix = col_end_matrix[,o,drop = FALSE]
  x = x[,o,,drop = FALSE]
  a = a[,o]
  v = v[o, ]
  v_inv = v_inv[, o]
  
  # SVD of column end matrix 
  svd_x = svd(col_end_matrix, nv = n, nu = 0)
  
  if (debug) {
    message('col_reduce: initial matrix')
    cat('column degrees:', col_degrees,'\n')
    print(col_end_matrix)
    cat('singular values of column end matrix:', svd_x$d,'\n')
    print(svd_x$v)
  }
  
  z = polm(c(0,1))
  iter = 0
  while (min(svd_x$d) < tol) {
    iter = iter + 1
    
    # Skip small entries at the end of svd_x$v[,n] (last singular value)
    k = is_large(svd_x$v[,n])
    
    if (debug) {
      message('col_reduce: reduce degree of column ',k)
    }
    
    v_step = v0
    v_step_inv = v0

    b = numeric(n)
    b[1:k] = svd_x$v[1:k,n]/svd_x$v[k,n]
    for (i in (1:k)) {
      v_step_inv[i, k] = b[i] * z^(col_degrees[k] - col_degrees[i])
      v_step[i, k] = -v_step_inv[i, k]
    }
    v_step[k,k] = 1
    
    a[, k] = a %r% v_step_inv[, k]
    v_inv[, k] = v_inv %r% v_step_inv[, k]
    v = v_step %r% v 
    
    {# balance the 'column norms'
    #     col_norm = apply(unclass(a), MARGIN = 2, FUN = l2_norm)
    # print(col_norm)
    #     a = a %r% diag(sqrt(mean(col_norm^2))/col_norm, nrow = n, ncol = n)
    #     v_inv = v_inv %r% diag(sqrt(mean(col_norm^2))/col_norm, nrow = n, ncol = n)
    #     v = diag(col_norm / sqrt(mean(col_norm^2)), nrow = n, ncol = n) %r% v
    }
    x = unclass(a)
    # eventually the degree of a has been reduced !?
    p = dim(x)[3] - 1
    
    # recompute col_degrees and col_end_matrix
    x[ , k, iseq(col_degrees[k]+1, p+1)] = 0 # this column has been purged!
    a = polm(x)
    
    col_norm = apply(matrix(x[, k, ], nrow = m, ncol = p+1), MARGIN = 2, FUN = l2_norm)
    col_degrees[k] = is_large(col_norm, tol = tol, count = TRUE) - 1
    if (col_degrees[k] < 0) {
      print(col_degrees)
      print(col_norm)
      print(x)
      stop('input "a" is (close to) singular')
    }
    x[, k, iseq(col_degrees[k]+2,p+1)] = 0
    col_end_matrix[,k] = x[, k, col_degrees[k] + 1]
    if (max(col_degrees) < p) {
      p = max(col_degrees)
      x = x[,,1:(p+1), drop = FALSE]
    }
    
    # (re) sort by column degrees 
    o = order(col_degrees, decreasing = FALSE)
    col_degrees = col_degrees[o]
    col_end_matrix = col_end_matrix[,o,drop = FALSE]
    x = x[,o,,drop = FALSE]
    a = polm(x)
    v = v[o, ]
    v_inv = v_inv[, o]
    
    # iterate 
    # SVD of column end matrix 
    svd_x = svd(col_end_matrix, nv = n, nu = 0)
    
    if (debug) {
      message('col_reduce: step=', iter)
      cat('column degrees:', col_degrees,'\n')
      print(col_end_matrix)
      cat('singular values of column end matrix:', svd_x$d,'\n')
      print(svd_x$v)
    }
  }
  
  # Resort column degrees in non-increasing direction
  o = order(col_degrees, decreasing = TRUE)
  col_degrees = col_degrees[o]
  col_end_matrix = col_end_matrix[,o,drop = FALSE]
  x = x[,o,,drop = FALSE]
  a = polm(x)
  v = v[o, ]
  v_inv = v_inv[, o]

  return(list(a = a, 
              v = v, v_inv = v_inv, 
              col_degrees = col_degrees, col_end_matrix = col_end_matrix))
}

#' Hermite Normal Form
#'
#' Calculate the \emph{column Hermite} (default) or \emph{row Hermite form} 
#' of a polynomial matrix \eqn{a(z)}, by using either (elementary) row operations (default) 
#' or column operations.
#'
#' For any \eqn{(m,n)} dimensional polynomial matrix \eqn{a(z)} with rank \eqn{r} 
#' (when considered as rational matrix) there exists a unimodular matrix \eqn{u(z)} 
#' and indices \eqn{j(1)<j(2)<\cdots<j(r)}{j(1)<j(2)<\dots<j(r)} such that 
#' \eqn{h(z)= u^{-1}(z) a(z)} is "quasi-upper-triangular" in the sense that 
#' \itemize{
#' \item \eqn{h_{ij(i)}}{h[i,j(i)]} is monic (the coefficient pertaining to the highest 
#'       degree is equal to one), 
#' \item the elements above \eqn{h_{ij(i)}}{h[i,j(i)]} have lower polynomial degree 
#'       than \eqn{h_{ij(i)}}{h[i,j(i)]} and 
#' \item \eqn{h_{i,j}}{h[i,j]} is zero for \eqn{i >r} or \eqn{j < j(i)}.
#' }
#' The matrix \eqn{h(z)} is called the row Hermite form of \eqn{a(z)}. The matrix \eqn{u^{-1}(z)} 
#' corresponds to the sequence of elementary row operations which renders \eqn{a(z)} into the 
#' desired upper-triangular form.
#' 
#' Quite analogously one may transform the matrix \eqn{a(z)} by elementary column operations into 
#' "quasi-lower-triangular" form \eqn{h(z) = a(z)u^{-1}(z)}. The corresponding normal form is 
#' called \emph{row Hermite form}. 
#' 
#' For a more detailed description, see e.g., 
#' \insertCite{Kailath80}{rationalmatrices} (page 375, Theorem 375) or 
#' the package vignette "Rational Matrices". 
#'
#' @param a Matrix polynomial, i.e. an object of class \code{\link{polm}}.
#' @param from_left Logical.
#'   Default set to TRUE, in which case unimodular row-transformations are used to obtain 
#'   the column Hermite normal form, i.e. \eqn{a(z) = u(z) h(z)}.
#'   If FALSE, unimodular column-transformations are used to obtain the row Hermite 
#'   normal form, i.e. \eqn{a(z) = h(z) u(z)}.
#' @param tol Tolerance parameter. Default set to \code{sqrt(.Machine$double.eps)}.
#' @param debug Logical. If TRUE, some diagnostic messages are printed.
#'
#' @return A list with the following slots.
#' \item{\code{h}}{\code{\link{polm}} object which represents the triangular matrix \eqn{h(z)}. 
#'                 Depending on \code{from_left} the matrix \eqn{h(z)} is either quasi-upper- 
#'                 or quasi-lower-triangular.}
#' \item{\code{u_inv}}{\code{\link{polm}} object, which represents the unimodular matrix 
#'                 which transform \eqn{a(z)} into the desired normal form, 
#'                 i.e. \eqn{h(z) = u^{-1}(z) a(z)}  or \eqn{h(z) = a(z)u^{-1}(z)}.}
#' \item{\code{u}}{\code{\link{polm}} object, which represents the unimodular matrix \eqn{u(z)} 
#'                 such that \eqn{a(z) = u(z) h(z)} or  \eqn{a(z) = h(z)u(z)}.}
#'
#'
#' @references    
#' \insertRef{Kailath80}{rationalmatrices}
#' 
#' @export
#'
#' @examples
#' #####################################################################
#' # Generate polynomial matrix
#' square = test_polm(dim = c(2,2), degree = 3)
#' # wide matrix, where all elements have a common factor (1-z)
#' wide = test_polm(dim = c(2,3), degree = 2) * polm(c(1,-1))
#' # tall matrix with a "right factor" ((2 x2) random polynomial matrix)
#' tall = test_polm(dim = c(3,2), degree = 1) %r% 
#'           test_polm(dim = c(2,2), degree = 1, random = TRUE, digits = 1)
#' 
#' a = tall  # choose one of the above cases
#'
#' ############################
#' # column Hermite form
#' out = hnf(a)
#' print(out$h, digits = 2, format = 'c')
#' 
#' # check result(s)
#' all.equal(a, prune(out$u %r% out$h))
#' all.equal(polm(diag(dim(a)[1])), prune(out$u_inv %r% out$u))
#' if (dim(a)[1] == dim(a)[2]) {
#'   rbind(sort(zeroes(a)), sort(zeroes(out$h)))
#' }
#'
#' ############################
#' # row Hermite form
#' out = hnf(a, from_left = FALSE)
#' print(out$h, digits = 2, format = 'c')
#'
#' # check result(s)
#' all.equal(a, prune(out$h %r% out$u))
#' all.equal(polm(diag(dim(a)[2])), prune(out$u_inv %r% out$u))
#' if (dim(a)[1] == dim(a)[2]) {
#'   rbind(sort(zeroes(a)), sort(zeroes(out$h)))
#' }
hnf = function(a, from_left = TRUE, tol = sqrt(.Machine$double.eps), debug = FALSE) {
  
  # Check inputs
  if (!inherits(a, 'polm')) {
    stop('input "a" is not a "polm" object!')
  }
  
  # skip zero leading coefficients
  a = prune(a, tol = tol)
  
  # Dimensions
  d = unname(dim(a))
  m = d[1]
  n = d[2]
  if (m*n == 0) {
    stop('input "a" is an empty polynomial matrix!')
  }
  
  
  # from_right
  if (!from_left) {
    a = t(a)
    
    # recompute dimensions
    d = unname(dim(a))
    m = d[1]
    n = d[2]
  }
  
  # Init of unimodular matrices 
  u0 = polm(diag(m))
  u = u0
  u_inv = u0
  
  i = 0
  j = 0
  pivots = integer(0)
  while ((i<m) && (j<n))
  {
    if (debug) {
      message('hnf pivot: i=', i,' j=',j,'\n')
    }
    
    p = degree(a)
    # code zero elelements with Inf
    p[p == -1] = Inf 
    
    if (all(is.infinite(p[(i+1):m,j+1]))) {
      # all remaining elements in (j+1)-th column are zero
      j = j+1 
      next
    }
    
    # purge (j+1)-th column
    out = purge_rc(a, pivot = c(i+1, j+1), direction = "down", permute = TRUE, 
                   tol = tol, monic = TRUE, debug = debug)
    a = out$h
    u = u %r% out$u
    u_inv = out$u_inv %r% u_inv
    
    # make sure that elements ABOVE the diagonal are smaller in degree than the diagonal element
    if (i > 0) {
      out = purge_rc(a, pivot = c(i+1,j+1), direction = "up", permute = FALSE, 
                     tol = tol, monic = TRUE, debug = debug)
      a = out$h
      u = u %r% out$u
      u_inv = out$u_inv %r% u_inv
    }
    
    pivots = c(pivots, j+1)
    i = i+1
    j = j+1
    
  }
  
  if (from_left){
    return(list(h = a, u = u, u_inv = u_inv, pivots = pivots, rank = length(pivots)))
  } else {
    return(list(h = t(a), u = t(u), u_inv = t(u_inv), pivots = pivots, rank = length(pivots)))
  }
}

#' Smith Normal Form
#'
#' Calculates the \emph{Smith normal form} of an \eqn{(m,n)}-dimensional polynomial matrix \eqn{a(z)},
#' i.e. a factorization of \eqn{a(z)} of the form 
#' \deqn{a(z) = u(z) s(z) v(z),}
#' where \eqn{u(z)} and \eqn{v(z)} are unimodular polynomial matrices of dimensions \eqn{(m,m)} and 
#' \eqn{(n,n)} respectively, and \eqn{s(z)} is a quasi-diagonal polynomial matrix of dimension \eqn{(m,n)}
#' with diagonal elements \eqn{d_i(z)}{d[i,i]} which satisfy
#' \itemize{
#' \item \eqn{d_{ii}}{d[i,i]} is monic for \eqn{i\leq r}{i \le r} and 
#'       \eqn{d_{ii}}{d[i,i]} is zero for \eqn{i>r}, and 
#' \item \eqn{d_{ii}}{d[i,i]} divides \eqn{d_{i+1,i+1}}{d[i+1,i+1]} for \eqn{i < r}.
#' }
#' Here \eqn{0\leq r \leq \min(m,n)}{0\le r \le min(m,n)} is the rank of \eqn{a(z)} 
#' when considered as a rational matrix.   
#' See *Gohberg, Lancaster, Rodman 09 - Matrix Polynomials* page 318 or \insertCite{Hannan.Deistler12}{rationalmatrices} page 42, Lemma 2.2.3.
#'
#' @param a Polynomial matrix, i.e. an object  of class \code{\link{polm}}.
#' @param tol Tolerance parameter, used for "pruning" the polynomial matrix (after each step). 
#'   See \code{\link{prune}}.
#' @param debug Logical, default to FALSE. If TRUE, then some diagnostic messages are printed.
#'
#' @return A list with five elements
#' \item{\code{s}}{Matrix polynomial of class \code{\link{polm}} and dimensions \eqn{(m,n)},
#'     representing \eqn{s(z)} in the Smith-form (whose only non-zero elements are on the diagonal).}
#' \item{\code{u}}{Unimodular matrix polynomial of class \code{\link{polm}} and dimensions \eqn{(m,m)},
#'     Represents \eqn{u(z)} in the Smith form \eqn{a(z) = u(z) s(z) v(z)}.}
#' \item{\code{u_inv}}{Unimodular matrix polynomial of class \code{\link{polm}} and dimensions \eqn{(m,m)},
#'     Represents the inverse of \eqn{u(z)}.}
#' \item{\code{v}}{Unimodular matrix polynomial of class \code{\link{polm}} and dimensions \eqn{(n,n)},
#'     Represents \eqn{v(z)} in the Smith form \eqn{a(z) = u(z) s(z) v(z)}.}
#' \item{\code{v_inv}}{Unimodular matrix polynomial of class \code{\link{polm}} and dimensions \eqn{(n,n)},
#'     Represents the inverse of \eqn{v(z)}.}
#'
#' @references    
#' \insertRef{Hannan.Deistler12}{rationalmatrices}
#' 
#' @name snf
#' @export
#' 
#' @examples
#' ##############
#' # Quadratic case
#'
#' a = test_polm(dim = c(2,2), degree = 1)
#' out = snf(a)
#'
#' print(out$s, digits = 2, format = 'c')
#'
#' all.equal(a, prune(out$u %r% out$s %r% out$v))
#'
#' ##############
#' # Tall case
#'
#' a = test_polm(dim = c(3,2), degree = 1)
#' out = snf(a)
#'
#' print(out$s, digits = 2, format = 'c')
#'
#' all.equal(a, prune(out$u %r% out$s %r% out$v))
#'
#' ##############
#' # Wide case
#'
#' a = test_polm(dim = c(2,3), degree = 1)
#' out = snf(a)
#'
#' print(out$s, digits = 2, format = 'c')
#'
#' all.equal(a, prune(out$u %r% out$s %r% out$v))
#'
#' ##############
#' # Diagonal case 
#' z = polm(c(0,1))
#' a = polm(diag(3))
#' a[3,3] = 1+2*z 
#' a[2,2] = a[3,3] * (1-z)
#' a[1,1] = a[2,2] * (1+z)
#' print(a, format = 'c')
#'
#' out = snf(a)
#'
#' print(out$s, digits = 2, format = 'c')
#'
#' all.equal(a, prune(out$u %r% out$s %r% out$v))
#' 
#' ##############
#' # Common factor(s) 
#' a = test_polm(dim = c(3,3), degree = 1, random = TRUE, digits = 1)
#' a = a * z
#' a[,2] = a[,2] * (1+z)
#' a[,1] = a[,2] * (1-z)
#' 
#' print(a, format = 'c')
#'
#' out = snf(a)
#'
#' print(out$s, digits = 2, format = 'c')
#'
#' all.equal(a, prune(out$u %r% out$s %r% out$v))
snf = function(a, tol = sqrt(.Machine$double.eps), debug = FALSE) {
  # Check inputs
  stopifnot("snf(): Input argument *a* must be a polm object" =  inherits(a, 'polm'))
  
  # skip "zero" leading coefficients 
  a = prune(a, tol = tol)
  # a0 = a      
  
  # Dimensions of a
  d = unname(dim(a))
  m = d[1]
  n = d[2]
  
  stopifnot("snf(): Input argument *a* must be a non-empty polynomial matrix!" = m*n != 0)

  # initialize unimodular matrices. a = u s v, u_inv is 
  u = polm(diag(m))
  u_inv = u
  v = polm(diag(n))
  v_inv = v
  
  i = 1
  iteration = 0
  while (i <= min(m,n)) {
    if (iteration > 100) stop('iteration maximum reached')
    
    # a is block diagonal
    # a[1:(i-1),1:(i-1)] is already in SNF form
    # handle the lower, right block: a[i:m, i:n]
    iteration = iteration + 1
    
    p = degree(a) # degree of each element (i,j) of the polynomial matrix
    if (debug) {
      message('snf: i=', i, ', iteration=', iteration)
      print(a, format = 'i|zj', digits = 2)
      print(p)
    }
    {# print(all.equal(a0, u %r% a %r% v))      
    # print(prune(u %r% u_inv, tol = tol))      
    # print(prune(v %r% v_inv, tol = tol))      
    }
    # code zero entries as Inf
    p[p == -1] = Inf
    
    # a[i:m, i:n] is zero => a is in SNF form!
    if (all(is.infinite(p[i:m, i:n]))) { 
      return(list(s = a, u = u, u_inv = u_inv, v = v, v_inv = v_inv))
    }
    
    # a[i,i] is non zero and all entries to the right and below the (i,i)-th element are zero 
    if (is.finite(p[i,i]) && 
        all(is.infinite(p[iseq(i+1,m), i])) && 
        all(is.infinite(p[i, iseq(i+1,n)])) ) {
      
      # At most bottom-right element: Make monic and return
      if (i == min(m,n)) {
        c = unclass(a)[i, i, p[i,i]+1]
        a[i, ] = a[i, ] %/% c
        u[, i] = u[, i] * c
        u_inv[i, ] = u_inv[i, ] %/% c
        return(list(s = a, u = u, u_inv = u_inv, v = v, v_inv = v_inv))
      }
      
      # Check remainder a[,] %% a[i,i] ####
      # (%% is polynomial remainder, %/%  division, see ?Ops.ratm)
      ra = test_polm(dim = c(m,n), degree = -1)
      ra[(i+1):m, (i+1):n] = prune(a[(i+1):m, (i+1):n] %% a[i,i], tol = tol)
      rp = degree(ra)
      rp[rp == -1] = Inf
      
      # all remainders are zero => next step i -> i+1
      if (all(is.infinite(rp))) {
        
        # first make diagonal element monic
        c = unclass(a)[i, i, p[i,i]+1]
        a[i, ] = a[i, ] %/% c
        u[, i] = u[, i] * c
        u_inv[i, ] = u_inv[i, ] %/% c
        
        i = i+1
        next
      }
      
      # find element with minimal (remainder) degree 
      c = which.min(apply(rp, MARGIN = 2, FUN = min))
      r = which.min(rp[, c])
      f = a[r,c] %/% a[i,i] # element-wise polynomial division
      
      # a[r,c] <- (a[r,c] %% a[i,i]) = (a[r,c] - f * a[i,i]), (%% gives remainder, %/% divides polynomial, and discards remainder)
      if (debug) {
        message('snf: a[', r, ',', c, '] <- a[', r, ',', c, '] %% a[i,i]\n')
      }

      {# Go through th code below with a diagonal matrix containing two different factors. First, add the (i,i) element to the zero element in row r, then use a column transformation (Euclidean algorithm) to subtract a factor times the (i,i)-element from the (r,c)-element
      
      # add i-th row to r-th row
      a[r,] = a[r,] + a[i,]
      u_inv[r,] = u_inv[r,] + u_inv[i,]
      u[,i] = u[,i] - u[,r]
      # substract f*(i-th column) from c-th column
      a[,c] = a[,c] - f * a[,i]
      v_inv[,c] = v_inv[,c] - f * v_inv[,i]
      v[i,] = v[i,] + f * v[c,]
      a = prune(a, tol = tol)
      }

      # next iteration
      next
    } 
    
    
    if (i > 1) {
      diag(p)[1:(i-1)] = Inf
    }
    
    # find element with minimal degree 
    c = which.min(apply(p, MARGIN = 2, FUN = min))
    r = which.min(p[, c])
    
    # bring this element to position (i,i)
    if (debug) {
      message('snf: a[i,i] <- a[', r, ',', c, ']\n')
    }
    
    rperm = 1:m
    rperm[c(i, r)] = c(r, i)
    cperm = 1:n
    cperm[c(i, c)] = c(c, i)
    a = a[rperm, cperm]
    p = p[rperm, cperm]
    u = u[, rperm]
    u_inv = u_inv[rperm, ]
    v = v[cperm, ]
    v_inv = v_inv[, cperm]
    
    # apply column purge
    out = purge_rc(a, pivot = c(i,i), direction = 'right', 
                   monic = FALSE, permute = TRUE, tol = tol, debug = debug)
    a = out$h 
    v = out$u %r% v
    v_inv = v_inv %r% out$u_inv
    
    # apply row purge 
    # note: this may generate non zero elements in the i-th column
    out = purge_rc(a, pivot = c(i,i), direction = 'down', 
                   monic = FALSE, permute = TRUE, tol = tol, debug = debug)
    a = out$h 
    u = u %r% out$u
    u_inv = out$u_inv %r% u_inv
    
    
    # next iteration
  }
  
  return(list(s = a, u = u, u_inv = u_inv, v = v, v_inv = v_inv))
}

# internal function 
# factorize a scalar polynomial into an stable and an unstable part
whf_scalar = function(a, tol = sqrt(.Machine$double.eps)) {
  a = prune(a, tol = tol)
  z = polyroot(unclass(a))
  a_f = a       # forward part (zeroes |z| < 1))
  a_b = polm(1) # backward part (zeroes |z| > 1))
  for (i in iseq(1,length(z))) {
    if (abs(z[i]) == 1) stop('unit roots are not allowed')
    if (abs(z[i]) > 1) {
      a_i = polm(c(-z[i], 1))
      a_f = a_f %/% a_i
      a_b = a_b * a_i
    }
  }
  a_f = prune(a_f, tol = tol)
  a_b = prune(a_b, tol = tol)
  if (is.complex(c(unclass(a_f), unclass(a_b)))) {
    stop('factors "a_f", "a_b" are complex')
  }
  return(list(a_f = a_f, a_b = a_b))
}

#' Wiener-Hopf Factorization
#'
#' A (Right-) Wiener-Hopf factorization (R-WHF) of a (square \eqn{(m,m)}-dimensional, non singular) polynomial matrix \eqn{A(z)} is a factorization of the form 
#' \deqn{A(z) = A_f(z) A_0(z) A_b(z) = A_r(z) A_b(z),}{A(z) = Af(z) A0(z) Ab(z) = Ar(z) Ab(z),}
#' where 
#' \describe{
#' \item{\eqn{A_b(z)}{Ab(z)}}{is a polynomial matrix which has only zeros outside the unit circle.}
#' \item{\eqn{A_r(z)}{Ar(z)}}{is a column reduced polynomial matrix with column degrees 
#'                            \eqn{\kappa_i}{\kappa[i]} and all zeroes inside the unit circle.}
#' \item{\eqn{A_0(z)}{A0(z)}}{is a diagonal matrix with diagonal entries 
#'                    \eqn{z^{\kappa_i}}{z^{\kappa[i]}} where \eqn{\kappa_i \geq \kappa_{i+1}}{\kappa[i] \ge \kappa[i+1]}}
#' \item{\eqn{A_f(z)}{Af(z)}}{is a polynomial in \eqn{z^{-1}}. 
#'                    Note that \eqn{A_f(z) = Ar(z) A_0^{-1}(z)}{Af(z) = Ar(z) A0^{-1}(z)}.}
#' }                                            
#' The factors \eqn{A_f(z), A_0(z), A_b(z)}{Af(z), A0(z), Ab(z)} are called 
#' \emph{forward}, \emph{null} and \emph{backward} components of \eqn{A(z)} and the integers 
#' \eqn{(\kappa_1,\ldots,\kappa_n)}{(\kappa[1],\dots,\kappa[n])} are the 
#' \emph{partial indices} of \eqn{A(z)}.
#' \cr
#' Similarly, the Left-WHF is defined as
#' \deqn{A(z) = A_b(z) A_0(z) A_g(z) = A_b(z) A_r(z),}{A(z) = Ab(z) A0(z) Af(z) = Ab(z) Ar(z),}
#' where \eqn{A_r(z)}{Ar(z)} is now row-reduced.
#' \cr
#' Note that zeroes on the unit circle are not allowed. 
#' In this case the procedure \code{whf()} throws an error.
#' \cr
#' The Wiener-Hopf factorization plays an important role for the analysis of linear, 
#' rational expectation models. See e.g. \insertCite{Al-Sadoon2017}{rationalmatrices}.
#' 
#' The algorithm is based on the Smith normal form (SNF), \link{snf}, and a column reduction step, see \link{col_reduce}.
#' An alternative is described in \insertCite{gohkaaspit03_summerschool}{rationalmatrices} pages 7ff.
#' 
#' @references
#' \insertRef{Al-Sadoon2017}{rationalmatrices}
#' \insertCite{gohkaaspit03_summerschool}{rationalmatrices}
#'
#' @param a \code{\link{polm}} object, which represents the polynomial matrix \eqn{A(z)}. 
#' @param right_whf Boolean. Default set to TRUE. 
#'                  If FALSE, then the left WHF \deqn{A(z) = A_b(z) A_0(z) A_f(z) = A_b(z) A_r(z),}{A(z) = Ab(z) A0(z) Af(z) = Ab(z) Ar(z),}
#'                  where the matrix \eqn{A_r(z)}{Ar(z)} is row-reduced.
#' @param tol Tolerance parameter, used for "pruning" the polynomial matrix (after each step). 
#'            See \code{\link{prune}}.
#' @param debug Logical, default set to FALSE. 
#'              If TRUE, then some diagnostic messages are printed.
#'
#' @return List with components 
#'         \itemize{
#'            \item{af: }{A \code{\link{lpolm}} object. Forward component \eqn{A_f}{Af}, a Laurent polynomial whose coefficients pertaining to positive powers are zero.} 
#'            \item{ab: }{A \code{\link{polm}} object. Backward component \eqn{A_b}{Ab}}
#'            \item{a0: }{A \code{\link{polm}} object. Diagonal matrix with monomials of degrees equal to the partial indices \eqn{A_0}{A0}}
#'            \item{ar: }{A \code{\link{polm}} object) the column reduced polynomial \eqn{A_r}{Ar}}
#'            \item{idx: }{ A vector of integers. partial indices \eqn{(\kappa_1,\ldots,\kappa_n)}{(\kappa[1],\dots,\kappa[n])}}
#'          }
#' @export
#' 
#' @examples
#' set.seed(1234) 
#' 
#' # create test polynomial
#' a = test_polm(dim = c(3,3), deg = 2, digits = 2, random = TRUE)
#' 
#' # compute WHF and print the result
#' out = whf(a)
#' print(out$af, digits = 2, format = 'c')
#' print(out$a0, digits = 2, format = 'c')
#' print(out$ab, digits = 2, format = 'c')
#' 
#' # check the result
#' all.equal(a, prune(out$ar %r% out$ab))           # A = Ar * Ab
#' 
#' # check A(z) = Ab(z^{-1}) A0(z) Ab(z)
#' # generate random complex z's
#' z = complex(real = rnorm(10), imaginary = rnorm(10))
#' a_z  = zvalues(a, z)         # A(z)
#' ab_z = zvalues(out$ab, z)    # Ab(z)
#' a0_z = zvalues(out$a0, z)    # A0(z)
#' af_z = zvalues(out$af, 1/z)  # Af(z^{-1})  
#' attr(af_z, 'z') = z           # in order to combine the 'zvalues' objects, 
#'                               # the attribute 'z' must be identical
#' all.equal(a_z, af_z %r% a0_z %r% ab_z)
#' 
#' all.equal(out$idx, degree(out$ar, 'columns'))    # idx = column degrees of Ar
#' all(svd(col_end_matrix(out$ar))$d > 1e-7)     # Ar is column reduced
#' abs(zeroes(out$ar, print_message = FALSE))       # Ar has zeroes inside the unit circle
#' abs(zeroes(out$ab, print_message = FALSE))       # Ab zeroes outside the unit circle
#' 
#' set.seed(NULL)
whf = function(a, right_whf = TRUE, tol = sqrt(.Machine$double.eps), debug = FALSE) {

  # Check inputs ####
  d = dim(a)
  stopifnot("whf(): Input *a* must be a polynomial matrix" = inherits(a, 'polm'),
            "whf(): Input *a* must be a square, non singular, polynomial matrix" = (d[1]*d[2]*d[3] != 0) && (d[1] == d[2]))
  m = d[1]
  
  # Left or right WHF?
  if (!right_whf){
    a = t(a)
  }
  
  # Obtain Smith normal form (in particular the diagonal matrix)
  snf = snf(a, tol = tol, debug = debug)
  
  {# cat('**************** SNF\n')
  # print(snf$s, digits = 3, format = 'i|zj')
  # print(snf$u, digits = 3, format = 'i|zj')
  # print(snf$v, digits = 3, format = 'i|zj')
  }
  
  # Factorize the diagonal matrix into stable and unstable part
  s_f = polm(diag(m))
  s_b = s_f 
  for (i in (1:m)) {
    out = whf_scalar(snf$s[i,i])
    s_f[i,i] = out$a_f
    s_b[i,i] = out$a_b
  }
  
  {# cat('************** diagonal\n')
  # print(snf$s, digits = 2, format = 'i|zj')
  # print(s_b*s_f, digits = 2, format = 'i|zj')
  }
  
  ar = snf$u %r% s_f
  ab = s_b %r% snf$v
  
  if (debug) {
    message('whf:')
    cat('col degree Ar', degree(ar, 'columns'),'\n')
    cat('row degree Ab', degree(ar, 'rows'),'\n')
  }
    
  {# cat('************** Ar Ab\n')
  # print(a, digits = 2, format = 'i|zj')
  # print(ar * ab, digits = 2, format = 'i|zj')
  }

  # column reduction of ar
  out = col_reduce(ar, tol = tol, debug = debug)
  # print(out)
  ar = out$a
  ab = out$v %r% ab
  idx = out$col_degrees

  # try to 'balance' Ar and Ab
  l2norm = function(x) sqrt(sum(x^2))
  nr = apply(unclass(ar), MARGIN = 2, FUN = l2norm) 
  nb = apply(unclass(ab), MARGIN = 1, FUN = l2norm) 
  nn = sqrt(nb/nr)
# print(rbind(nr,nb,nn))
  ar = prune(ar %r% diag(nn, nrow = m, ncol = m), tol = tol)
  ab = prune(diag(nn^{-1}, nrow = m, ncol = m) %r% ab, tol = tol)
{# nr = apply(unclass(ar), MARGIN = 2, FUN = l2norm) 
# nb = apply(unclass(ab), MARGIN = 1, FUN = l2norm) 
# nn = sqrt(nb/nr)
# print(rbind(nr,nb,nn))
  } 
  a0 = polm(diag(m))
  z = polm(c(0,1))
  for (i in (1:m)) {
    if (idx[i] > 0) a0[i,i] = z^idx[i]
  }
  
  # create the forward part Af 
  # multiply the j-th column with z^(-idx[j]) => polynomial in z^(-1)
  ar_tmp = unclass(ar)
  af = array(0, dim = dim(ar_tmp))
  idx_max = max(idx)
  for (j in (1:m)) {
    af[,j,(1+idx_max-idx[j]):(1+idx_max)] = ar_tmp[,j,1:(idx[j]+1), drop = FALSE]
  }
  af = lpolm(af, min_deg = -idx_max)
  af = prune(af, tol = tol)
  
  # if left- WHF transform elements back
  if (!right_whf){
    af = t(af)
    ab = t(ab)
    ar = t(ar)
  }
  
  return(list(af = af, a0 = a0, ab = ab, ar = ar, idx = idx))
}



# Laurent polynomial transformations ####

#' Transforms to Polynomial in Forward Shift
#' 
#' Transform \code{\link{polm}} object to polynomial in forward shift (represented as \code{\link{lpolm}} object), i.e. transform
#' \deqn{a(z) = a_0 + a_1 z^1 + \cdots + a_p z^p}{
#'       a(z) = a[0] + a[1] z^1 + \dots + a[p] z^p}
#' to       
#' \deqn{a(z) = a_0 + a_1 z^{-1} + \cdots + a_p z^{-p}}{
#'       a(z) = a[0] + a[1] z^(-1) + \dots + a[p] z^(-p)}
#' 
#' @param polm_obj \code{\link{polm}} object
#'
#' @return \code{\link{lpolm}} object
#' @export
#'
#' @examples
#' (p = test_polm(degree = 3))
#' polm2fwd(p)
polm2fwd = function(polm_obj){
  x = unclass(polm_obj)
  d = dim(x)
  if (d[3] > 0){
    return(lpolm(x[,,d[3]:1], min_deg = -(d[3]-1)))
  } else {
    return(lpolm(x, min_deg = 0))
  }
}

#' Forward and Backward Bracket
#' 
#' For an \code{\link{lpolm}} object as input, \code{get_bwd} discards all coefficient matrices pertaining to \strong{negative} powers and returns a \code{\link{lpolm}} object with \code{min_deg = 0}.
#' Similarly, \code{get_fwd} discards all coefficient matrices pertaining to \strong{non-negative} powers, and also returns an \code{\link{lpolm}} object.
#' 
#' Obtain the forward or backward part of a Laurent polynomial, i.e. apply \eqn{[.]_-} or \eqn{[.]_+} to 
#' \deqn{a(z) = a_{-q} z^{-q} + \cdots + a_{-1} z^{-1} + a_0 + a_1 z^1 + \cdots + a_p z^p}{
#'       a(z) = a[-q] z^(-q) + \dots + a[-1] z^{-1} + a[0] + a[1] z^1 + \dots + a[p] z^p}
#' and obtain for \code{get_fwd}
#' \deqn{[a(z)]_- = a_{-q} z^{-q} + \cdots + a_{-1} z^{-1}}{
#'       [a(z)]_- = a[-q] z^(-q) + \dots + a[-1] z^{-1}}
#' or for \code{get_bwd}
#' \deqn{[a(z)]_+ = a_0 + a_1 z^1 + \cdots + a_p z^p}{
#'       [a(z)]_+ = a[0] + a[1] z^1 + \dots + a[p] z^p}
#'       
#' @param lpolm_obj Laurent polynomial object \code{\link{lpolm}}
#'
#' @return Laurent polynomial object \code{\link{lpolm}} without non-negative coefficients or without negative
#' @export
#'
#' @examples
#' (lp = test_lpolm(degree_max = 2, degree_min = -2))
#' get_fwd(lp)
get_fwd = function(lpolm_obj){
  md = attr(lpolm_obj, which = "min_deg")
  if (md >=0){
    return(lpolm_obj)
  } else {
    x = unclass(lpolm_obj)[,,1:(-md), drop = FALSE]
    return(lpolm(x, min_deg = md))
  }
}

#' @rdname get_fwd
#' @export
#' @examples 
#' (lp = test_lpolm(degree_max = 2, degree_min = -2))
#' get_bwd(lp)
#' 
#' (lp = lpolm(1:3, min_deg = 2))
#' get_bwd(lp)
#' 
#' (lp = lpolm(1:3, min_deg = -1))
#' get_bwd(lp)
#' 
#' (lp = lpolm(1:3, min_deg = -5))
#' get_bwd(lp)
get_bwd = function(lpolm_obj){
  attr_l = attributes(lpolm_obj)
  min_deg = attr_l$min_deg
  d = attr_l$dim
  attributes(lpolm_obj) = NULL
  dim(lpolm_obj) = d
  if (d[3]+1 < -min_deg){
    # No non-negative coefficients: Empty lpolm
    return(lpolm(array(0, c(d[1],d[2],0)), min_deg = 0))
  } else if (min_deg >= 0){
    polm_offset = array(0, dim = c(dim(lpolm_obj)[1:2], min_deg))
    return(lpolm(dbind(d = 3, polm_offset, lpolm_obj), min_deg = 0))
  } else {
    return(lpolm(lpolm_obj[,,(-min_deg+1):d[3], drop = FALSE], min_deg = 0))
  }
}



