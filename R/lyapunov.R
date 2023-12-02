# lyapunov.R 
# 
# These are essentially R wrapper functions for the Rcpp routines!

#' Lyapunov Equation
#'
#' This function solves the Lyapunov equation
#' \deqn{P = A P A' + Q}
#' where \eqn{A,Q} are real valued, square matrices and \eqn{Q} is symmetric.
#' The Lyapunov equation has a unique solution if 
#' \eqn{\lambda_i(A)\lambda_j(A) \neq 1}{\lambda_i(A)\lambda_j(A) != 1} holds 
#' for all eigenvalues of \eqn{A}.  
#' If \eqn{A} is stable (i.e. the spectral radius of \eqn{A} is less than one) and
#' \eqn{Q} is positive semidefinite then the solution \eqn{P} is also positive semidefinite.
#' \cr
#' The procedure uses the Schur decomposition(s) of \eqn{A} and computes 
#' the solution "column by column",
#' see \insertCite{Kitagawa77,Hammarling1982}{rationalmatrices}.
#'
#' @param A,Q  \eqn{(m,m)} matrices. Note that the routine silently assumes that  
#'        \eqn{Q} is symmetric (and hence the solution \eqn{P} is also symmetric). 
#' @param non_stable (character string) indicates what to do, when \eqn{A} is not stable.
#' @param attach_lambda (boolean) if yes, then the eigenvalues of \eqn{A} are
#'        attached to the solution \eqn{P} as an attribute.
#'
#' @return P (\eqn{(m,m)} matrix) the  solution of the Lyapunov equation. 
#'
#' @references
#' \insertRef{Kitagawa77}{rationalmatrices}
#'
#' \insertRef{Hammarling1982}{rationalmatrices}
#'
#' @export
#'
#' @examples
#' # A is stable and Q is positve definite
#' m = 4
#' A = diag(runif(m, min = -0.9, max = 0.9))
#' V = matrix(rnorm(m*m), nrow = m, ncol = m)
#' A = V %*% A %*% solve(V)
#' B = matrix(rnorm(m*m), nrow = m, ncol = m)
#' Q = B %*% t(B)
#' P = lyapunov(A, Q)
#' all.equal(P, A %*% P %*% t(A) + Q)
#'
#' # unstable matrix A
#' A = diag(runif(m, min = -0.9, max = 0.9))
#' A[1,1] = 2
#' V = matrix(rnorm(m*m), nrow = m, ncol = m)
#' A = V %*% A %*% solve(V)
#' P = lyapunov(A, Q)
#' all.equal(P, A %*% P %*% t(A) + Q)
#' # note that the solution P (in general) is not positive semidefinite
#' eigen(P, only.values= TRUE, symmetric = TRUE)$values
#'
#' # attach the eigenvalues of A to the solution P
#' P = lyapunov(A, Q, attach = TRUE)
#' print(P)
#'
#' # issue a warning message
#' P = lyapunov(A, Q, non_stable = 'warn')
#'
#' \dontrun{
#' # throw an error
#' P = lyapunov(A, Q, non_stable = 'stop')
#' }
lyapunov = function(A, Q, 
                    non_stable = c("ignore", "warn", "stop"),
                    attach_lambda = FALSE) {
  # check inputs 
  if ( (!is.numeric(A)) || (!is.matrix(A)) || (nrow(A) != ncol(A)) ) {
    stop('"A" must be a square, numeric  matrix ')
  }
  
  if ( (!is.numeric(Q)) || (!is.matrix(Q)) || any(dim(Q) != dim(A)) ) {
    stop('"Q" must be a numeric matrix with the same dimension as "A"')
  }
  
  non_stable = match.arg(non_stable)
  attach_lambda = as.logical(attach_lambda)[1]
  
  m = ncol(A)
  if (m == 0) {
    stop('A,Q are "empty"')
  }

  P = matrix(0, nrow = m, ncol = m)
  lambda_r = numeric(m)
  lambda_i = numeric(m)
  stop_if_non_stable = (non_stable == 'stop')
  
  # call RcppArmadillo routine
  is_stable = .Call(`_rationalmatrices_lyapunov_cpp`, A, Q, P, lambda_r, lambda_i, stop_if_non_stable)
  
  if ((stop_if_non_stable) && (!is_stable)) stop('"A" matrix is not stable')
  if ((non_stable == 'warn') && (!is_stable)) warning('"A" matrix is not stable')
  if (attach_lambda) {
      attr(P,'lambda') = complex(real = lambda_r, imaginary = lambda_i)
    }

  return(P)
}


#' Jacobian of the Solution of the Lyapunov Equation
#'
#' This (internal helper) function considers the solution of a Lyapunov equation
#' \deqn{P = A P A' + Q}
#' where \eqn{A,Q} are real valued, square matrices and \eqn{Q} is symmetric. 
#' The \emph{directional derivative} of \eqn{P} along (the matrices) \eqn{dA}, \eqn{dQ} 
#' is given by the solution of the Lyapunov Equation: 
#' \deqn{dP = A dP A' + dA P A' + A P dA' + dQ}
#' 
#' @param A,Q  \eqn{(m,m)} matrices. Note that the routine silently assumes that  
#'        \eqn{Q} is symmetric (and hence the solution \eqn{P} is also symmetric). 
#' @param dA,dQ  \eqn{(m^2,n)} matrices. Each column of \eqn{dA}, \eqn{dQ} determines a 
#'        direction along which the derivative is computed. Note that the routine 
#'        silently assumes that each column of \eqn{dQ} represents a symmetric matrix
#'        (and hence \eqn{dP} is also symmetric). 
#' @param non_stable (character string) indicates what to do, when \eqn{A} is not stable.
#'
#' @return List with the slots
#' \item{P}{(\eqn{(m^2,n)}-dimensional matrix) Solution of the lyapunov equation}
#' \item{J}{(\eqn{(m^2,n)}-dimensional matrix) Jacobian of the vectorised solution
#'          of the Lyapunov equation. Each column of \eqn{J} is the 
#'          directional derivative of \eqn{vec(P)} along the respective columns 
#'          of \eqn{dA} and \eqn{dQ}.} 
#' \item{lambda}{Eigenvalues of \eqn{A}.}
#' \item{is_stable}{(boolean) Is \eqn{A} stable or not?}
#'
#' @seealso \code{\link{lyapunov}}. 
#' 
#' @export
#' @keywords internal 
#'
#' @examples
#' m = 5
#' 
#' A = matrix(rnorm(m^2), nrow = m, ncol = m)
#' Q = crossprod(matrix(rnorm(m^2), nrow = m, ncol = m))
#' P = lyapunov(A, Q)
#' all.equal(Q, P - A %*% P %*% t(A))
#' 
#' n = 6
#' dA = matrix(rnorm((m^2)*n), nrow = m^2, ncol = n)
#' dQ = matrix(rnorm((m^2)*n), nrow = m^2, ncol = n)
#' for (i in (1:n)) {
#'   # make sure that each column of dQ corresponds to a symmetric matrix!
#'   junk = matrix(dQ[,i], nrow = m, ncol = m)
#'   junk = junk + t(junk)
#'   dQ[,i] = junk
#' }
#' out = lyapunov_Jacobian(A, Q, dA, dQ)
#' all.equal(out$P, P)
#'  
#' eps = 1e-8
#' theta = rnorm(n)
#' matrix(out$J %*% theta, nrow = m, ncol = m)
#' 
#' # compute the derivative via "finite differences"
#' dP = lyapunov(A + matrix(dA %*% theta, nrow = m, ncol = m)*eps, 
#'               Q + matrix(dQ %*% theta, nrow = m, ncol = m)*eps)
#' all.equal(matrix(out$J %*% theta, nrow = m, ncol = m), 
#'           (dP - P)/eps, scale = mean(abs(out$J)), tol = 1e-6)
lyapunov_Jacobian = function(A, Q, dA, dQ, 
                             non_stable = c("ignore", "warn", "stop")) {
  # check inputs 
  if ( (!is.numeric(A)) || (!is.matrix(A)) || (nrow(A) != ncol(A)) ) {
    stop('"A" must be a square, numeric  matrix ')
  }
  
  if ( (!is.numeric(Q)) || (!is.matrix(Q)) || any(dim(Q) != dim(A)) ) {
    stop('"Q" must be a numeric matrix with the same dimension as "A"')
  }
  
  m = nrow(A)
  
  if ( (!is.numeric(dA)) || (!is.matrix(dA)) || (nrow(dA) != (m^2) ) ) {
    stop('"dA" must be a (m^2, k) dimensional numeric  matrix ')
  }
  
  if ( (!is.numeric(dQ)) || (!is.matrix(dQ)) || any(dim(dA) != dim(dQ)) ) {
    stop('"dQ" is not compatible with "A" or "dA"')
  }
  
  n = ncol(dA)
  
  non_stable = match.arg(non_stable)

  if (m*n == 0) {
    stop('A, Q, dA or dQ are "empty"')
  }
  
  P = matrix(0, nrow = m, ncol = m)
  J = matrix(0, nrow = m^2, ncol = n)
  lambda_r = numeric(m)
  lambda_i = numeric(m)
  stop_if_non_stable = (non_stable == 'stop')
  
  # call RcppArmadillo routine
  is_stable = .Call(`_rationalmatrices_lyapunov_Jacobian_cpp`, A, Q, P, 
                    dA, dQ, J, lambda_r, lambda_i, stop_if_non_stable)
  
  if ((stop_if_non_stable) && (!is_stable)) stop('"A" matrix is not stable')
  if ((non_stable == 'warn') && (!is_stable)) warning('"A" matrix is not stable')

  return(list(P = P, J = J, lambda = complex(real = lambda_r, imaginary = lambda_i), is_stable = is_stable))
}