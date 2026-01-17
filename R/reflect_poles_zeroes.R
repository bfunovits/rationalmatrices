#' Division Algorithm for Polynomial Matrices
#'
#' For given polynomial matrices \eqn{a(z), b(z)} compute two matrices \eqn{c(z), d(z)} such that 
#' \deqn{a(z) = c(z) b(z) + d(z)}
#' where the degree of \eqn{d(z)} is smaller than the degree of \eqn{b(z)}. 
#' The matrix \eqn{b(z)} must be square with a non singular leading coefficient matrix! 
#' The matrices must be compatible, i.e. the number of columns of \eqn{a(z)} must equal the number of rows (and columns) of \eqn{b(z)}.  
#'
#' @param a,b Two compatible polynomial matrices.  
#'
#' @return List with two slots 
#' \itemize{
#'   \item{qu}{contains the polynomial \eqn{c(z)}}
#'   \item{rem}{contains the polynomial \eqn{d(z)}.}
#' }
#' 
#' @export
#'
#' @examples
#' a = test_polm(dim = c(3,2), degree = 4, random = TRUE)
#' b = test_polm(dim = c(2,2), degree = 2, random = TRUE)
#' 
#' (out = polm_div(a, b))
#' 
#' all.equal(a, out$qu %r% b + out$rem)
polm_div = function(a, b) {
  
  # a,b must be 'polm'-objects
  if ((!inherits(a, 'polm')) || (!inherits(b, 'polm'))) {
    stop('The input parameter "a", "b" must be "polm" objects!')
  }
  a = unclass(a)
  dim.a = dim(a)
  b = unclass(b)
  dim.b = dim(b)
  
  if ((length(dim.a) !=3) || (length(dim.b) != 3) || 
      (dim.a[2] != dim.b[1]) || (dim.b[1] != dim.b[2]) || (dim.b[1] == 0)) {
    stop('The polynomial matrices a, b are not compatible!')
  }
  
  m = dim.a[1]
  n = dim.a[2]
  p = dim.a[3]-1
  q = dim.b[3]-1
  
  if (q > p) {
    c = polm(array(0, dim = c(m, n, 1)))
    d = polm(a)
    return(list(qu = c, rem = d))
  }
  if (q == -1) {
    stop('The polynomial matrix b is zero!')
  }
  a = matrix(a, nrow = m, ncol = n*(p+1))
  b = matrix(b, nrow = n, ncol = n*(q+1))  
  c = matrix(0, nrow = m, ncol = n*((p-q)+1))
  
  inv_bq = try(solve(b[ , (q*n+1):((q+1)*n), drop = FALSE])) 
  if (inherits(inv_bq, 'try-error')) {
    # print(b[ , (q*n+1):((q+1)*n), drop = FALSE])
    stop('The leading coefficient of the polynomial matrix b is singular!')
  }
  
  for (i in ((p-q):0)) {
    # print(i)
    c[ , (i*n+1):((i+1)*n)] = 
      a[ , ((q+i)*n+1):((q+i+1)*n), drop = FALSE] %*% inv_bq
    if (q > 0) {
      # print(a[ , (i*n+1):((i+q)*n), drop = FALSE])
      # print(c[ , (i*n+1):((i+1)*n), drop = FALSE])
      # print(b[, 1:(q*n), drop = FALSE])
      a[ , (i*n+1):((i+q)*n)] = a[ , (i*n+1):((i+q)*n), drop = FALSE] - 
        c[ , (i*n+1):((i+1)*n), drop = FALSE] %*% b[, 1:(q*n), drop = FALSE]
    }
  }
  
  c = polm(array(c, dim = c(m,n,p-q+1)))
  if (q > 0) {
    d = polm(array(a[,1:(q*n)], dim = c(m,n,q)))
  } else {
    d = polm(array(0, dim = c(m,n,1)))
  }
  
  return(list(qu = c, rem = d))
}


#' Blaschke Factors
#'
#' The Blaschke factor at \eqn{\alpha}{\alpha} is the rational function
#' \deqn{B(z) := \frac{1-\bar{\alpha}z}{-\alpha + z}}{
#'       B(z) := (1-Conj(\alpha)z) / (-\alpha+z)}
#' This is an all-pass function with a pole at \eqn{z=\alpha} and a zero at \eqn{z=1/\bar{\alpha}}{z=1/Conj(\alpha)}. 
#' The function \code{blaschke(alpha)} returns this rational \eqn{(1 \times 1)}{(1 x 1)} matrix in \code{\link{lmfd}} form. 
#' Clearly \eqn{B(z)} has complex coefficients, if \eqn{\alpha} is complex.
#' \cr
#' The call \code{blaschke2(alpha, row=NULL)} computes the product of the Blaschke factors at \eqn{\alpha} and at \eqn{\bar{\alpha}}{Conj(\alpha)}, i.e. the rational function
#' \deqn{B_{s}(z) := 
#' \frac{1-2\Re(\alpha)z + |\alpha|^2 z^2}{|\alpha|^2 -2\Re(\alpha) z + z^2}}{
#'       Bs(z) := 
#' (1-2Re(\alpha)z + Mod(\alpha)^2 z^2) / (Mod(\alpha)^2 -2Re(\alpha) z + z^2)}
#' \cr
#' If \code{blaschke2} is called with an optional argument \code{w} (a non zero complex vector of length 2) then  \code{blaschke2} constructs a \eqn{(2 \times 2)}{(2 x 2)} rational, all-pass matrix of the form 
#' \deqn{B_{2}(z) := a^{-1}(z) b(z)}{B2(z) := a^{-1}(z) b(z)} 
#' where \eqn{a(z), b(z)} are  two \eqn{(2 \times 2)}{(2 x 2)} polynomial matrices (with real coefficients) of degree one. 
#' This matrix is constructed such that the column space of \eqn{a(\alpha)} is spanned by  \eqn{\bar{w}}{Conj(w)} and the column space of \eqn{a(\bar{\alpha})}{a(Conj(\alpha))} is spanned by the vector \eqn{w}. 
#'  
#' @note The routine \code{blaschke2} throws an error if \eqn{\alpha} is not complex (i.e. the imaginary part is zero). 
#'   If \eqn{\alpha} is close to the unit circle then \code{blaschke2(alpha, w)} simply returns an \code{\link{lmfd}} representation of the bivariate identity matrix. 
#'   If \eqn{w} and \eqn{\bar{w}}{Conj(w)} are almost linearly dependent, then an error is thrown.
#' 
#' @param alpha complex or real scalar, represents \eqn{\alpha}.
#' @param w \code{NULL} or a (complex) vector of length 2. 
#' @param tol Tolerance (used to decide whether \code{alpha} 
#'        has modulus equal to one).
#'
#' @return \code{\link{lmfd}} object, which represents the constructed "Blaschke factors" \eqn{B(z)}, \eqn{B_s(z)}{Bs(z)} or \eqn{B_2(z)}{B2(z)}.
#' @name blaschke
#' @export
#'
#' @examples
#' # Blaschke factor with a real alpha
#' (B = blaschke(1.5))
#' zvalues(B) %>% abs()
#' 
#' # Blaschke factor with a complex alpha
#' (B = blaschke(complex(real = 1.5, imaginary = 0.5)))
#' zvalues(B) %>% abs()
#' 
#' # product of the Blaschke factors at alpha and Conj(alpha) 
#' # this gives a scalar, rational, all-pass matrix with real coefficients
#' (B = blaschke2(complex(real = 1.5, imaginary = 0.5)))
#' zvalues(B) %>% abs()
#' 
#' #############################################################
#' # a "bivariate" Blaschke factor 
#' 
#' # case 1: alpha is "outside the unit circle" ################
#' (alpha = complex(real = 1.5, imaginary = 0.5))
#' (w = complex(real = c(0.1,0.9), imaginary = c(0.75,-0.5)))
#' 
#' (B = blaschke2(alpha, w = w))
#' # B(z) is all-pass 
#' print(zvalues(B) %r% Ht(zvalues(B)), digits = 3)
#' 
#' # B(z) has poles at z=alpha, z=Conj(alpha) and 
#' # zeroes at z=1/alpha and z=1/Conj(alpha)
#' poles(B)
#' zeroes(B)
#' 
#' # The column space of a(alpha) is spanned by the vector Conj(w).
#' max(abs( Conj(c(-w[2], w[1])) %*% zvalue(B$a, alpha) ))
#' 
#' # case 2: alpha is "inside the unit circle" #################
#' (alpha = 1 / complex(real = 1.5, imaginary = 0.5))
#' (w = complex(real = c(0.1,0.9), imaginary = c(0.75,-0.5)))
#' 
#' (B = blaschke2(alpha, w = w))
#' # B(z) is all-pass 
#' print(zvalues(B) %r% Ht(zvalues(B)), digits = 3)
#' 
#' # B(z) has poles at z=alpha, z=Conj(alpha) and 
#' # zeroes at z=1/alpha and z=1/Conj(alpha)
#' poles(B)
#' zeroes(B)
#' 
#' # The column space of a(alpha) is spanned by the vector Conj(w).
#' max(abs( Conj(c(-w[2], w[1])) %*% zvalue(B$a, alpha) ))
#' 
#' # case 3: alpha is "on the unit circle" #####################
#' alpha = alpha / Mod(alpha)
#' blaschke2(alpha) %>% print(digits = 2)
#' blaschke2(alpha, w = w) %>% print(digits = 2)
blaschke = function(alpha) {
  alpha = as.vector(alpha)[1]
  BM = lmfd(polm(c(-alpha,1)), polm(c(1, -Conj(alpha))))
  return(BM)
}

#' @name blaschke 
#' @export
blaschke2 = function(alpha, w = NULL, tol = 100*.Machine$double.eps) {
  alpha = as.vector(alpha)[1]
  if (Im(alpha) == 0) {
    stop('"alpha" must be complex with a non zero imaginary part!')
  }
  
  # Univariate case
  if (is.null(w)) { 
    BM = lmfd(polm(c(Mod(alpha)^2, -2*Re(alpha), 1)), 
              polm(c(1, -2*Re(alpha), Mod(alpha)^2)))
    return(BM)
  }

  # root 'alpha' is on the unit circle: simply return the identity!
  if (abs( Mod(alpha) - 1 ) < tol) {
    return(lmfd(polm(diag(2)), polm(diag(2))))
  } 
  
  # w must be two dimensional complex vector and 
  # w and Conj(w) must be linearly independent
  if (Mod(alpha) < 1) {
    # root alpha is inside the unit circle
    # a(z) = (-A + Iz)
    A = try(t(solve( t( cbind(w, Conj(w)) ), 
                     t( cbind(alpha*w, Conj(alpha*w))) )), silent = TRUE)
    if (inherits(A, 'try-error')) {
      stop('w and Conj(w) are linearly dependent!')
    }
    A = Re(A)
    Gamma0 = lyapunov(t(A), diag(2), non_stable = 'ignore')
    # b(z) = (I - Bz)T^{-1}
    B = solve(Gamma0, t(A) %*% Gamma0)
    T = chol(Gamma0 - t(B) %*% Gamma0 %*% B)
    
    BM = lmfd(polm(array(cbind(-A, diag(2)), dim = c(2,2,2))), 
              polm(array(cbind(diag(2), -B), dim = c(2,2,2))) %r% solve(T))
    return(BM)
    
  } else {
    # root alpha is outside the unit circle
    # a(z) = (I - Az)
    A = try( t(solve( t( cbind(alpha*w, Conj(alpha*w)) ), 
                      t( cbind(w, Conj(w))) )), silent = TRUE)
    if (inherits(A, 'try-error')) {
      stop('w and Conj(w) are linearly dependent!')
    }
    A = Re(A)
    Gamma0 = lyapunov(t(A), diag(2), non_stable = 'ignore')
    # b(z) = (B - Iz)T^{-1}
    B = solve(Gamma0, t(A) %*% Gamma0)
    T = chol(Gamma0 - t(B) %*% Gamma0 %*% B)
    
    BM = lmfd(polm(array(cbind(diag(2), -A), dim = c(2,2,2))), 
              polm(array(cbind(B, -diag(2)), dim = c(2,2,2))) %r% solve(T))
    return(BM)
  }
} 

#' Create an All-Pass Rational Matrix
#' 
#' Create a square, all-pass rational matrix in statespace form with given \eqn{A,B}. 
#'
#' @param A square \eqn{(s,s)} matrix 
#' @param B \eqn{(s,m)} matrix
#'
#' @return A \code{\link{stsp}} object which represents the all-pass 
#'         rational matrix.
#' @export
#' @keywords internal
#' 
#' @section Notes: 
#' The function is intended as an internal helper function and thus does not check the inputs. See also \code{\link{reflect_poles}} and \code{\link{reflect_zeroes}}. 
#' 
#' If \eqn{A,C} is given we can use \code{t(make_allpass(t(A), t(C)))}. 
#'
#' @examples
#' m = 2
#' s = 6
#' A = matrix(rnorm(s*s), nrow = s, ncol = s)
#' B = matrix(rnorm(s*m), nrow = s, ncol = m)
#' 
#' K = make_allpass(A, B)
#' all.equal(cbind(A,B), cbind(K$A, K$B))
#' # check that K(z) is all-pass
#' Kf = zvalues(K)
#' all.equal(zvalues(polm(diag(m))) + complex(imaginary = 0), 
#'           Kf %r% Ht(Kf) + complex(imaginary = 0))
#' 
#' # if (A,C) is given, we proceed as follows:
#' C = matrix(stats::rnorm(m*s), nrow = m, ncol = s)
#' K = t(make_allpass(t(A), t(C)))
#' all.equal(rbind(A,C), rbind(K$A, K$C))
#' # check that K(z) is all-pass
#' Kf = zvalues(K)
#' all.equal(zvalues(polm(diag(m))) + complex(imaginary = 0), 
#'           Kf %r% Ht(Kf) + complex(imaginary = 0))
make_allpass = function(A, B) {
  m = ncol(B)
  s = nrow(B)
  
  P = lyapunov(A, B %*% t(B))
  C = -t(solve(A %*% P, B))
  # 
  M = diag(m) - t(solve(A, B)) %*% t(C)
  # make sure that M is symmetric 
  M = (M + t(M))/2
  D = solve(t(chol(M)))
  C = D %*% C
  # print(D %*% M %*% t(D))
  
  x = stsp(A, B, C, D)
  return(x)
}


#' Reflect Zeroes of Rational Matrices
#' 
#' Given a square, rational matrix \eqn{x(z)} and a set of zeroes of \eqn{x(z)}, 
#' the function \code{reflect_zeroes} constructs a rational, allpass matrix 
#' \eqn{U(z)} such that \eqn{x(z)U(z)} has zeroes which are the reciprocals 
#' of the selected zeroes. The other zeroes of \eqn{x(z)} and the poles of 
#' \eqn{x(z)} are not changed. 
#'  
#' @param x Square, rational matrix object (with real coefficients).
#' @param zeroes Complex or real vector of zeroes to be reflected. 
#'        It is assumed that for complex conjugated pairs of roots, 
#'        only \bold{one} is contained in this vector.
#' @param check_zeroes If \code{TRUE} then the procedure checks that the given 
#'        zeroes are zeroes of the rational matrix \eqn{x(z)}.
#' @param tol (double) Tolerance parameter used for the checks. 
#' @param ... Not used.
#' 
#' @seealso The procedure uses the helper functions \code{\link{blaschke}} 
#'          and \code{\link{make_allpass}}. For reflecting poles, see 
#'          \code{\link{reflect_poles}}. 
#'          
#' @note 
#' The procedures are only implemented for rational matrices with real coefficients. 
#' Therefore complex zeroes occur in complex conjugated pairs and such pairs are **jointly** reflected to ensure that the result again has real coefficients. 
#' Note however, that the argument \code{zeroes} must only contain \bold{one} 
#' element of the pairs to be reflected.
#' 
#' In some degnerate cases the procedure(s) may run into numerical problems: 
#' For polynomial matrices and matrices in \code{lmfd} form, 
#' the procedure assumes that the value of the matrix evaluated at 
#' \eqn{z=0} is non singular. For the matrices in statespace form, zeroes 
#' on (or close to) the unit circle are not allowed. In addition multiple 
#' zeroes may lead to unreliable results. 
#'
#' @rdname reflect_zeroes
#' 
#' @return Rational matrix object which represents the rational matrix 
#'         \eqn{x(z)U(z)} with the reflected zeroes. The object is of 
#'         the same class as the input object \code{x}.
#' @export
reflect_zeroes = function(x, zeroes, ...) {
  UseMethod("reflect_zeroes")
}

#' @rdname reflect_zeroes
#' @export
#'
#' @examples
#' # ###################################################################
#' # polynomial matrix 
#' 
#' # construct a (2 x 2) polynomial matrix with degree 2
#' x = polm(array(c(1,0,
#'                  0,1,
#'                  -1,0,
#'                  0,0,
#'                  1,0.25,
#'                  0.75,-1), dim = c(2,2,3)))
#' x
#' 
#' # determine zeroes of x(z)
#' ( x_zeroes = zeroes(x) )
#' 
#' # evaluate the polynomial at the (real) zeroe x_zeroes[1] 
#' ( x_value = zvalue(x, x_zeroes[1]) )
#'
#' # Verify that the matrix evaluated at z = x_zeroes[1] is rank deficient
#' d = svd(x_value)$d
#' (min(d) < (max(d) * sqrt(.Machine$double.eps)))
#'
#' # Reflect this zeroe at the unit circle and check the zeroes of the result
#' x1 = reflect_zeroes(x, x_zeroes[1])
#' 
#' r_zeroes = x_zeroes
#' r_zeroes[1] = 1 / r_zeroes[1]
#' x1_zeroes = zeroes(x1)
#' j = match_vectors(r_zeroes, x1_zeroes)
#' all.equal(r_zeroes, x1_zeroes[j]) 
#' 
#' # x_zeroes[2] and x_zeroes[3] form a pair of complex conjugated zeroes
#' ( x_value = zvalue(x, x_zeroes[2]) )
#'
#' # Verify that the matrix evaluated at x_zeroes[2] is rank deficient
#' d = svd(x_value)$d
#' (min(d) < (max(d) * sqrt(.Machine$double.eps)))
#'
#' # reflect the real zeroe x_zeroes[1] and this pair of 
#' # complex conjugated zeroes at the unit circle and verify the result
#' x1 = reflect_zeroes(x, x_zeroes[c(1,2)])
#' 
#' r_zeroes = x_zeroes
#' r_zeroes[1:3] = 1 / r_zeroes[1:3]
#' x1_zeroes = zeroes(x1)
#' j = match_vectors(r_zeroes, x1_zeroes)
#' all.equal(r_zeroes, x1_zeroes[j]) 
#' 
#' # Check that the transformation matrix U (x1 = x %r% U) is all-pass
#' all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
reflect_zeroes.polm = function(x, zeroes, tol = sqrt(.Machine$double.eps), 
                               check_zeroes = TRUE, ...) {
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('The argument "x" must be a square polynomial matrix (in polm form).')
  }
  m = d[1]
  
  zeroes = as.vector(zeroes)
  k = length(zeroes)
  if (k == 0) {
    # nothing to do
    return(x)
  }
  
  if (check_zeroes) {
    z = zeroes(x, tol = tol, print_message = FALSE)
    zz = c(zeroes, Conj(zeroes[Im(zeroes) != 0]))
    if (length(z) < length(zz)) {
      stop(paste('The polynomial "x" has less zeroes than specified',
                 'in the argument "zeroes"!'))
    }
    j = match_vectors(zz, z)
    if (max(abs(zz - z[j])) > tol) {
      print(cbind(zz, z[j], zz - z[j]))
      stop(paste('The specified zeroes do not match',
                 'with the zeroes of the polynomial matrix "x"!'))
    }
  }
  
  # should we sort/order the zeroes?
  for (i in (1:k)) {
    z0 = zeroes[i]
    # cat('zero:', z0, '\n')
    
    if (Im(z0) == 0) {
      # real zero
      z0 = Re(z0)
      x0 = zvalue(x, z0)
      
      U = svd(x0, nu = 0, nv = m)$v
      x = x %r% U
      B = blaschke(z0)
      x[, m] = (polm_div(x[, m], B$a)$qu ) %r% B$b
    } else {
      # complex zero
      x0 = zvalue(x, z0)
      
      # flip z0
      U = svd(x0, nu = 0, nv = m)$v
      x = x %r% U
      B = blaschke(z0)
      # print(zvalue(x, z0)[, m])
      x[, m] = (polm_div(x[, m], B$a)$qu) %r% B$b
      
      # flip conjugate zeroe Conj(z0) and the corresponding null vector Conj(w)
      w0 = t(Conj(U)) %*% Conj(U[, m])
      w0[m] = w0[m] / zvalue(B, Conj(z0))
      
      U = svd(w0, nu = m, nv = 0)$u
      x = x %r% U
      B = blaschke(Conj(z0))
      # print(zvalue(x, Conj(z0))[, 1])
      x[, 1] = (polm_div(x[, 1], B$a)$qu) %r% B$b
      
      # make real!
      # since we use x(z0) w0 = 0 <=> x(Conj(z0)) Conj(w0) = 0
      x0 = zvalue(x, 0)
      L = t(chol(Re(x0 %*% t(Conj(x0)))))
      x = x %r% solve(x0, L)
      x = Re(x)
    }
  }
  return(x)
}


#' @rdname reflect_zeroes
#' @export
#' 
#' @examples 
#' 
#' 
#' # ###################################################################
#' # rational matrix in LMFD form
#' 
#' set.seed(12345)
#' (x = test_lmfd(dim = c(3,3), degrees = c(2,2), digits = 2))
#' (x_zeroes = zeroes(x))
#' 
#' # reflect all zeroes inside the unit circle
#' # note: for complex zeroes, select only one of the complex conjugated pair!
#' x1 = reflect_zeroes(x, x_zeroes[(abs(x_zeroes) <1) & (Im(x_zeroes) >=0 )])
#' 
#' r_zeroes = x_zeroes
#' r_zeroes[abs(r_zeroes) < 1] = 1 / r_zeroes[abs(r_zeroes) < 1]
#' (x1_zeroes = zeroes(x1))
#' j = match_vectors(r_zeroes, x1_zeroes)
#' all.equal(r_zeroes, x1_zeroes[j]) 
#' 
#' # Check that the transformation matrix U (x1 = x %r% U) is all-pass
#' all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#' 
#' set.seed(NULL)
reflect_zeroes.lmfd = function(x, zeroes, tol = sqrt(.Machine$double.eps), 
                               check_zeroes = TRUE, ...) {
  if (!is.numeric(x)) {
    stop('The argument "x" must be a rational matrix with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('The argument "x" must be a square polynomial matrix (in polm form).')
  }
  m = d[1]
  
  zeroes = as.vector(zeroes)
  k = length(zeroes)
  
  if (k == 0) {
    # nothing to do
    return(x)
  }
  
  x = lmfd(a = x$a, 
           b = reflect_zeroes(x$b, zeroes, check_zeroes = check_zeroes, 
                              tol = tol))
  return(x)
}


# Internal helper for reflect_poles.stsp and reflect_zeroes.stsp
# Validates input and preprocesses roots (poles/zeroes)
# Returns NULL if roots is empty (caller should return x early)
# Otherwise returns the roots with complex conjugates appended
validate_stsp_reflect_input <- function(x, roots, tol, root_name) {
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('argument "x" must be a square rational matrix (in stsp form).')
  }

  roots = as.vector(roots)
  if (length(roots) == 0) {
    return(NULL)
  }

  if (min(abs(abs(roots) - 1)) < tol) {
    stop(paste0('one of the selected ', root_name, ' has modulus close to one.'))
  }

  # append complex conjugates
  c(roots, Conj(roots[Im(roots) != 0]))
}

#' @rdname reflect_zeroes
#' @export
#'
#'
#' @examples
#'
#'
#' # ###################################################################
#' # rational matrix in statespace form (stsp object)
#'
#' # create a random (2,2) rational matrix in state space form with
#' # state dimension s=5
#' set.seed(12345)
#' (x = test_stsp(dim = c(2,2), s = 5))
#' # zeroes of x(z)
#' (x_zeroes = zeroes(x))
#'
#' # reflect all unstable zeroes (inside the unit circle)
#' # note: for complex zeroes, select only one of the complex conjugated pair!
#' x1 = reflect_zeroes(x, x_zeroes[(abs(x_zeroes) <1) & (Im(x_zeroes) >=0 )])
#'
#' r_zeroes = x_zeroes
#' r_zeroes[abs(r_zeroes) < 1] = 1 / r_zeroes[abs(r_zeroes) < 1]
#' (x1_zeroes = zeroes(x1))
#' j = match_vectors(r_zeroes, x1_zeroes)
#' all.equal(r_zeroes, x1_zeroes[j])
#'
#' # Check that the transformation matrix U (x1 = x %r% U) is all-pass
#' all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#'
#' set.seed(NULL)
reflect_zeroes.stsp = function(x, zeroes, tol = sqrt(.Machine$double.eps), ...) {
  zeroes = validate_stsp_reflect_input(x, zeroes, tol, "zeroes")
  if (is.null(zeroes)) return(x)

  A = x$A
  B = x$B
  C = x$C
  D = x$D
  s = nrow(B)
  m = ncol(B)
  
  # transform (A - BD^{-1} C) to upper block diagonal matrix, 
  # where the top block corresponds to the selected zeroes 
  # schur decomposition of (A - BD^{-1} C)
  out = try(schur(A - B %*% solve(D, C), 1/zeroes))
  if (inherits(out, 'try-error')) stop('ordered schur decomposition failed.')
  
  # (A - B D^{-1} C) = U S U' 
  A = t(out$U) %*% A %*% out$U
  B = t(out$U) %*% B
  C = C %*% out$U
  k = out$k
  
  i1 = (1:k)
  i2 = iseq((k+1), s)
  
  # create all-pass function with given A,C
  U = t(make_allpass(t(out$S[i1, i1, drop = FALSE]), 
                     t(solve(D, C[, i1,drop = FALSE]))))
  
  Ah = rbind( cbind(A[i2, i2, drop = FALSE], A[i2, i1, drop = FALSE]),
              cbind(A[i1, i2, drop = FALSE], A[i1, i1, drop = FALSE]) )
  Bh = rbind( B[i2, , drop = FALSE] %*% U$D, B[i1, , drop = FALSE] %*% U$D + U$B )
  Ch = cbind( C[, i2, drop = FALSE], C[, i1, drop = FALSE] )
  Dh = D %*% U$D
  
  return(stsp(Ah, Bh, Ch, Dh))
}



#' Reflect Poles of Rational Matrices
#' 
#' Given a square, rational matrix \eqn{x(z)} and a set of poles of \eqn{x(z)}, 
#' the function \code{reflect_poles} constructs an allpass matrix \eqn{U(z)} 
#' such that \eqn{x(z)U(z)} has poles which are the reciprocals of the 
#' selected poles and the other poles are not changed. 
#' Also the zeroes are not changed, i.e. \eqn{x} and \eqn{xU} have the same zeroes. 
#' 
#' @note 
#' The procedures are only implemented for rational matrices with real coefficients. 
#' Therefore complex poles occur in complex conjugated pairs and such pairs are
#' \bold{jointly} reflected to ensure that the result again has real coefficients. 
#' Note however, that the argument \code{poles} must only contain \bold{one} element 
#' of the pairs to be reflected.
#'  
#' In some degnerate cases the procedure(s) may run into numerical problems: 
#' For matrices in \code{rmfd} form, 
#' the procedure assumes that the matrix has no poles at \eqn{z=0}. 
#' For the matrices in statespace form, poles 
#' on (or close to) the unit circle are not allowed. In addition multiple 
#' poles may lead to unreliable results. 
#'
#' @seealso The procedure uses the helper functions \code{\link{blaschke}} 
#'          and \code{\link{make_allpass}}. For reflecting zeroes, 
#'          see \code{\link{reflect_zeroes}}.
#'          
#' @param x Square, rational matrix object (with real coefficients).
#' @param poles Complex or real vector of poles to be reflected. 
#'        It is assumed that for complex conjugated pairs of poles, 
#'        only \bold{one} is contained in this vector.
#' @param ... Not used.
#' @param check_poles If \code{TRUE} then the procedure checks that the 
#'        given poles are poles of the rational matrix x(z).
#' @param tol (double) Tolerance parameter used for the checks.
#'
#' @return Rational matrix object which represents the rational matrix 
#'         \eqn{x(z)U(z)} with the reflected poles. The object is of the 
#'         same class as the input object \code{x}.
#'         
#' @rdname reflect_poles
#' @export
reflect_poles = function(x, poles, ...) {
  UseMethod("reflect_poles")
}

#' @rdname reflect_poles
#' @export
#' 
#' @examples
#' # ###################################################################
#' # rational matrix in statespace form ('stsp' object)
#'
#' set.seed(12345)
#' # create random (2,2) rational matrix in state space form 
#' # with state dimension s=5
#' ( x = test_stsp(dim = c(2,2), s = 5) )
#' # poles of x(z)
#' ( x_poles = poles(x) )
#' 
#' # reflect all unstable poles (inside the unit circle) ###########
#' # note: for complex zeroes, select only one of the complex conjugated pair!
#' x1 = reflect_poles(x, poles = x_poles[(abs(x_poles) < 1) & (Im(x_poles) >= 0)])
#' 
#' r_poles = x_poles
#' r_poles[abs(r_poles) < 1] = 1 / r_poles[abs(r_poles) < 1]
#' (x1_poles = poles(x1))
#' j = match_vectors(r_poles, x1_poles)
#' all.equal(r_poles, x1_poles[j]) 
#' 
#' # Check that the transformation matrix U (x1 = x %r% U) is all-pass
#' all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#' 
#' set.seed(NULL)
reflect_poles.stsp = function(x, poles, tol = sqrt(.Machine$double.eps), ...) {
  poles = validate_stsp_reflect_input(x, poles, tol, "poles")
  if (is.null(poles)) return(x)

  A = x$A
  B = x$B
  C = x$C
  D = x$D
  s = nrow(B)
  m = ncol(B)
  
  # transform A matrix to lower block diagonal matrix, 
  # where the top diagonal block corresponds to the selected poles!
  # ordered schur decomposition of A'
  out = try(schur(t(A), 1/poles))
  if (inherits(out, 'try-error')) stop('ordered schur decomposition failed.')
  
  # A' = U S U' => A = U S' U' 
  A = t(out$S)
  B = t(out$U) %*% B
  C = C %*% out$U
  k = out$k
  
  i1 = (1:k)
  i2 = iseq((k+1), s)
  
  # create all-pass function
  U = make_allpass(A[i1, i1, drop = FALSE], B[i1,,drop = FALSE])
  U = U^{-1}
  
  Ah = rbind( cbind(A[i2, i2, drop = FALSE], 
                    A[i2, i1, drop = FALSE] + B[i2, ,drop = FALSE] %*% U$C),
              cbind(matrix(0, nrow = k, ncol = s-k), U$A) )
  Bh = rbind( B[i2, , drop = FALSE] %*% U$D, U$B )
  Ch = cbind( C[, i2, drop = FALSE], C[ , i1, drop = FALSE] + D %*% U$C )
  Dh = D %*% U$D
  
  return(stsp(Ah, Bh, Ch, Dh))
}

#' @rdname reflect_poles
#' @export
#' 
#' @examples
#' 
#' 
#' # ###################################################################
#' # rational matrix in RMFD form
#'
#' set.seed(12345)
#' # create random (2 x 2) rational matrix in RMFD form with degrees (2,1)
#' ( x = test_rmfd(dim = c(2,2), degree = c(2,1)) )
#' # poles of x(z)
#' ( x_poles = poles(x) )
#' 
#' # reflect all unstable poles (inside the unit circle) ###########
#' # note: for complex zeroes, select only one of the complex conjugated pair!
#' x1 = reflect_poles(x, poles = x_poles[(abs(x_poles) < 1) & (Im(x_poles) >= 0)])
#' 
#' r_poles = x_poles
#' r_poles[abs(r_poles) < 1] = 1 / r_poles[abs(r_poles) < 1]
#' (x1_poles = poles(x1))
#' j = match_vectors(r_poles, x1_poles)
#' all.equal(r_poles, x1_poles[j]) 
#' 
#' # Check that the transformation matrix U (x1 = x %r% U) is all-pass
#' all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#' 
#' set.seed(NULL)
reflect_poles.rmfd = function(x, poles, tol = sqrt(.Machine$double.eps), 
                              check_poles = TRUE, ...) {
  # check inputs ....
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('argument "x" must be a square rational matrix (in stsp form).')
  }
  
  poles = as.vector(poles)
  k = length(poles)
  if (k == 0) {
    # nothing to do
    return(x)
  }
  
  c = x$c
  d = x$d
  
  c = t(reflect_zeroes(t(c), poles, check_zeroes = check_poles, 
                              tol = tol))
  return(rmfd(c = c, d= d))
}

#' Polynomial Roots as List
#'
#' This helper function coerces a vector of (complex or real) roots into a list. 
#' For a complex conjugate pair of roots, only the one with a positive 
#' imaginary part is retained. 
#' 
#' The routine assumes that complex roots appear, up to some small numerical
#' errors, in complex conjugate pairs. Therefore the procedure tries to match
#' the roots with their complex conjugates. If this is not possible, then an
#' error is thrown.
#' 
#' Roots, which are classified as real, are replaced by their real part.
#' Roots, which are classified as complex, are replaced by the mean of the 
#' root and the best matching conjugate root.
#' 
#' @seealso The matching of conjugate roots is done by the internal 
#'          helper function \code{\link{match_vectors}}.
#' 
#' @param roots Vector of cplx or doubles. Obtained e.g. from a call 
#'        to \code{\link{zeroes}}.
#' @param tol Double. Tolerance parameter used to decide whether roots 
#'        and conjugate roots "match". 
#'
#' @return List of roots. For complex conjugated pairs of roots, only the ones 
#'         with positive part are returned. The procedure orders the roots 
#'         according to their imaginary part, and thus real roots come first.  
#' @export
#' 
#' @examples
#' set.seed(12345)
#' p = 5
#' a = rnorm(p+1)   # coefficients of a random polynomial a(z) of degree p = 5
#' z = polyroot(a)  # compute the roots of a(z)
#' z
#' 
#' # try to match roots and conjugate roots
#' j = match_vectors(z, Conj(z))
#' # z is approximately equal to Conj( z[j] )
#' print(data.frame(z = z, j = j, `Conj(z[j])` = Conj(z[j]), d = z-Conj(z[j])))
#' # z[1] and z[5] are complex conjugates (up to numerical errors)
#' # z[2] and z[3] are complex conjugates (up to numerical errors)
#' # z[4] is real (up to numerical errors)
#' 
#' # coerce the vector "z" to a list
#' (z_list = roots_as_list(z))
#' # the first slot contains the real root (and thus is of class "numeric")
#' # z_list[[1]] = Re(z[4])
#' # the slots 2 and 3 contain the complex roots and thus are of class "complex"
#' # z_list[[2]] = (z[2] + Conj(z[3])/2
#' # z_list[[3]] = (z[1] + Conj(z[5])/2
#' 
#' # The routine zeroes() uses the function eigen() 
#' # (to compute the eigenvalues of the companion matrix) 
#' # and thus returns exact conjugate pairs:
#' (z = zeroes(polm(a)))
#' 
#' # match roots and conjugate roots
#' j = match_vectors(z, Conj(z))
#' # z is equal to Conj(z[j])
#' print(j)
#' all.equal(z, Conj(z[j]))
#' 
#' # coerce the vector "z" to a list
#' (z_list = roots_as_list(z))
#' 
#' set.seed(NULL)
roots_as_list = function(roots, tol = sqrt(.Machine$double.eps)) {
  roots = as.vector(roots)
  if (is.numeric(roots)) {
    return(as.list(roots))
  }
  if (is.complex(roots)) {
    p = length(roots)
    if (p == 0) return(as.list(roots))
    
    j = match_vectors(roots, Conj(roots))
    i = (1:p)
    if (max(abs(roots[i] - Conj(roots[j]))) > tol) {
      stop('could not match pairs of complex conjugate roots')
    }
    
    # if roots[k] is a "real" root, 
    #   then j[k] = k should hold 
    # if roots[k1], roots[k2] is a pair of complex conjugate roots, 
    #   then j[k1] = k2 and j[k2] = k1 should hold.  
    # Together this means jj = j[j] = (1:p) must hold!
    jj = j[j]
    
    if (any(jj != i)) {
      stop('could not match pairs of complex conjugate roots')
    }
    
    # make sure that "real" roots have imaginary part = 0 
    # and that pairs of roots are complex conjugates of each other 
    roots = (roots + Conj(roots[j]))/2 
    
    # skip roots with negative imaginary part
    roots = roots[Im(roots) >= 0]
    # order by imaginary part, => real roots come first
    roots = roots[order(Im(roots), Re(roots))]
    
    # convert to list
    roots = as.list(roots)
    
    # convert roots with zero imaginary part to 'numeric'
    roots = lapply(roots, function(x) ifelse(Im(x) > 0, x, Re(x)))
    
    return(roots)
  }
  stop('"roots" must be a vector of class "numeric" or "complex"!')
}
