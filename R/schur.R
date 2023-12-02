# schur.R
# 
# Schur decomposition and related tools and methods

#' Schur Decomposition
#' 
#' The Schur decomposition of a (real, square) matrix \eqn{A} is 
#' \deqn{A = U S U'}
#' where \eqn{U} is an orthogonal matrix and \code{S} is 
#' upper quasi-triangular with 1-by-1 and 2-by-2 blocks on the diagonal. 
#' The 2-by-2 blocks correspond to the complex eigenvalues of \eqn{A}. 
#' 
#' The optional parameter \code{select} determines which eigenvalues, respectively 
#' 1-by-1 and 2-by-2 blocks, should be put to the top of the matrix \eqn{S}. 
#' 
#' \itemize{
#' \item{\code{select='iuc'}: Select the eigenvalues with moduli less than one 
#'       (inside the unit circle).}
#' \item{\code{select='ouc'}: Select the eigenvalues with moduli greater than one 
#'       (outside the unit circle).}
#' \item{\code{select='lhf'}: Select the eigenvalues with a negative imaginary part 
#'       (left half plane).}
#' \item{\code{select='rhf'}: Select the eigenvalues with a positive imaginary part 
#'       (right half plane).}
#' \item{\code{select='real'}: Select the real eigenvalues.}
#' \item{\code{select='cplx'}: Select the complex eigenvalues.}
#' \item{\code{select} is a (complex) vector of eigenvalues. The function checks whether
#'       the "target eigenvalues", i.e. the entries of \code{select},  match the eigenvalues 
#'       of the matrix \eqn{A}. In addition the procedure also makes sure 
#'       that complex eigenvalues are selected in complex conjugated pairs.} 
#' }
#' 
#' The function \code{schur} is simply a wrapper for \code{\link[QZ]{qz.dgees}} 
#' and \code{\link[QZ]{qz.dtrsen}}. 
#' 
#' @param A square, non-empty, real matrix.
#' @param select If non NULL, then the Schur decomposition is reordered such that 
#'        the selected cluster of eigenvalues appear in the leading 
#'        diagonal blocks of the quasi-triangular matrix \eqn{S}. See the details below.
#' @param tol tolerance used to decide whether the target eigenvalues match the 
#'        eigenvalues of \eqn{A}. See the details below.
#'
#' @return List with components \code{U}, \code{S}, \code{lambda} and \code{k}, 
#'         where \code{lambda} contains the (computed) eigenvalues of \eqn{A} and 
#'         \code{k} indicates how many eigenvalues have been selected (put to the top).
#' @export
#'
#' @examples
#' # generate a "random" 6-by-6 matrix A
#' m = 6
#' set.seed(1532)
#' A = matrix(stats::rnorm(m*m, sd = 0.5), nrow = m, ncol = m)
#' set.seed(NULL)
#' 
#' # compute the Schur decomposition of A (and its eigenvalues)
#' out = schur(A)
#' lambda = out$lambda # eigenvalues of A
#' 
#' # check A = U S U' 
#' all.equal(A, out$U %*% out$S %*% t(out$U))
#' print(out$S)
#' print(lambda)
#' 
#' # compute an "ordered" Schur decomposition where the eigenvalues 
#' # inside the unit circle are put to the top of S:
#' out = schur(A, 'iuc')
#' print(out$S)
#' print(out$k) # three eigenvalues are inside the unit circle.
#' print(out$lambda) 
#' 
#' # compute an "ordered" Schur decomposition where the eigenvalues 
#' # lambda[5] and lambda[6] apear in the top. Note that 
#' # lambda[5] is complex and hence the procedure also selects 
#' # the conjugate of lambda[5]:  
#' out = schur(A, lambda[c(6,5)])
#' print(out$S)
#' print(out$k) # three eigenvalues have been selected
#' print(out$lambda) 
#' 
#' \dontrun{
#' # If the "target" eigenvalues do not match the eigenvalues of A 
#' # then "schur" throws an error:
#' out = schur(A, select = lambda[1:3]+ 1)
#' }
schur = function(A, select = NULL, tol = sqrt(.Machine$double.eps)) {
  if ( !is.numeric(A) || !is.matrix(A) || (nrow(A) != ncol(A)) ) {
    stop('input "A" must be a square, complex or real valued matrix.')
  }
  if ( any(!is.finite(A)) ) stop('input "A" contains missing or infinite elements.')
  
  m = nrow(A)
  if (m == 0) stop('input "A" has zero rows/columns.')
  
  out = QZ::qz.dgees(A)
  if (out$INFO != 0) stop('Schur decomposition failed. Error code of "dgees.f": ', out$INFO)
  lambda = complex(real = out$WR, imaginary = out$WI) # eigenvalues 

  selected = logical(m)
  if (!is.null(select)) {
    if (is.character(select)) {
      selected = switch(select,
                        iuc = (abs(lambda) < 1), 
                        ouc = (abs(lambda) > 1), 
                        lhp = (out$WR < 0), 
                        rhp = (out$WR > 0),
                        real = (out$WI == 0),
                        cplx = (out$WI != 0))
    } else {
      select = as.vector(select)
      if ( length(select) > m ) stop('The input vector "select" has more than ', m, ' entries.')
      if (length(select) > 0) {
        C = abs(matrix(select, nrow = length(select), ncol = m) - 
                  matrix(lambda, nrow = length(select), ncol = m, byrow = TRUE))
        match = munkres(C)
        if ( match$c > length(select)*tol ) stop('could not match the target eigenvalues to the eigenvalues of "A"!')
        selected[match$a[,2]] = TRUE 
        
        # make sure that complex conjugated pairs are both selected
        i = which( (out$WI > 0) & (selected) ) # positive imaginary parts come first
        selected[i+1] = TRUE
        i = which( (out$WI < 0) & (selected) ) # positive imaginary parts come first
        selected[i-1] = TRUE
      }
    }
  }
  k = sum(selected)
  if ((k > 0) && (k < m)) {
    # reorder diagonal blocks 
    out = QZ::qz.dtrsen(out$T, out$Q, selected, job = "N", want.Q = TRUE)
    if (out$INFO != 0) stop('Reordering of Schur decomposition failed. Error code of "dtrsen.f": ', out$INFO)
    lambda = complex(real = out$WR, imaginary = out$WI) # eigenvalues 
  }
  return(list(S = out$T, U = out$Q, lambda = lambda, k = k))
}
