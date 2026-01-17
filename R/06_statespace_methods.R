# ===================================================================
# Functional Area #6: State-Space Tools (CLAUDE.md)
#
# Purpose: Implement state-space specific operations and analysis
#   - Grammians: Controllability and observability Grammians
#   - Balancing: Transform to balanced realization (diagonalize Grammians)
#   - Balanced truncation: Approximate high-order systems
#   - Controllability/Observability matrices: System property analysis
#   - System identification and realization algorithms
#
# Related Files: realization_10_utilities_tools.R (Hankel matrices), 07_numerical_lyapunov.R (Lyapunov equation)
# ===================================================================

# stsp_methods. R #########################################
# special methods/operations on statespace realizations

#' Controllability and Observability Matrix
#' 
#' The controllability matrix of a statespace realisation \eqn{k(z)=C(I-Az^{-1})^{-1}B + D} is 
#' the matrix
#' \deqn{
#' [B,AB,\dots,A^{o-1}B]
#' }
#' and the observability matrix is 
#' \deqn{
#' [C',A'C',\dots,(A')^{o-1}C']'
#' }
#' 
#'
#' @param A either a \code{\link{stsp}} object or a square \eqn{(s,s)} dimensional matrix. 
#' @param B \eqn{(s,n)} dimensional matrix. This argument is ignored if \code{A} is a 
#'          \code{stsp} object.
#' @param C \eqn{(m,s)} dimensional matrix. This argument is ignored if \code{A} is a 
#'          \code{stsp} object.
#' @param o (non negative) integer. The default value is \eqn{o=s}.  
#'
#' @return Controllability or observability matrix.
#' @export
ctr_matrix = function(A, B, o = NULL) {
  if (missing(A)) stop('parameter A is missing')
  if (inherits(A, 'stsp')) {
    B = A$B
    A = A$A
    s = nrow(B)
    n = ncol(B)
  } else {
    if ( !( is.numeric(A) || is.complex(A) ) ) stop('parameter A is not numeric or complex')
    if ( (!is.matrix(A)) || (nrow(A) != ncol(A)) ) stop('parameter A is not a square matrix')
    s = nrow(A)
    if (missing(B)) stop('parameter B is missing')
    if ( !( is.numeric(B) || is.complex(B) ) ) stop('parameter B is not numeric or complex')
    if ( (!is.matrix(B)) || (nrow(B) != s) ) stop('parameters A,B are not compatible')
    n = ncol(B)
  }
  if (is.null(o)) o = s
  o = as.integer(o)[1]
  if (o < 0) stop('o must be a non negative integer')
  
  Cm = array(0, dim = c(s,n,o))
  Cm[,,1] = B
  for (i in iseq(2,o)) {
    B = A %*% B
    Cm[,,i] = B
  }
  dim(Cm) = c(s,o*n)
  return(Cm)
}


#' @export
#' @rdname ctr_matrix
obs_matrix = function(A, C, o = NULL) {
  if (missing(A)) stop('parameter A is missing')
  if (inherits(A, 'stsp')) {
    C = A$C
    A = A$A
    s = ncol(C)
    m = nrow(C)
  } else {
    if ( !( is.numeric(A) || is.complex(A) ) ) stop('parameter A is not numeric or complex')
    if ( (!is.matrix(A)) || (nrow(A) != ncol(A)) ) stop('parameter A is not a square matrix')
    s = nrow(A)
    if (missing(C)) stop('parameter C is missing')
    if ( !( is.numeric(C) || is.complex(C) ) ) stop('parameter C is not numeric or complex')
    if ( (!is.matrix(C)) || (ncol(C) != s) ) stop('parameters A,C are not compatible')
    m = nrow(C)
  }
  if (is.null(o)) o = s
  o = as.integer(o)[1]
  if (o < 0) stop('o must be a non negative integer')
  
  Om = array(0, dim = c(s,m,o))
  A = t(A)
  C = t(C)
  Om[,,1] = C
  for (i in iseq(2,o)) {
    C = A %*% C
    Om[,,i] = C
  }
  dim(Om) = c(s,o*m)
  return(t(Om))
}


#' Check the Minimality of a Statespace Realization
#' 
#' Check whether a statespace realization of a rational matrix is minimal. 
#' The procedure constructs the Hankel matrix of the impulse response coefficients with \eqn{s} block rows and \eqn{s} block columns, 
#' where \eqn{s} is the statespace dimension of the given statespace realization. 
#' If this statespace realization is minimal then the Hankel matrix has rank \eqn{s}. 
#' Therefore the procedure returns \code{TRUE} if the \eqn{s}-th singular values of the Hankel matrix is larger than \code{tol}. 
#' 
#' The procedure does not check whether the statespace realization is \emph{observable} and/or 
#' \emph{controllable}. To this end one may compute the observability/controllability matrices (\code{\link{obs_matrix}}, \code{\link{ctr_matrix}}) or the corresponding Grammians (\code{\link{grammians}}). 
#'
#' @param x (\code{\link{stsp}} object) rational matrix in statespace form.
#' @param tol	a tolerance parameter, which is used to decide the rank of 
#'              the Hankel matrix of the impulse response coefficients.
#' @param only.answer	if TRUE, just return a logical (TRUE or FALSE). 
#'                    Otherwise a list with additional information is returned.
#' @param ... not used.
#' 
#' @seealso \code{\link{pseries}} and \code{\link{bhankel}}.
#' 
#' @note
#' This procedure returns different objects, depending on the parameter \code{only.answer}.
#'
#' @return If \code{only.answer} is true then a logical (\code{TRUE} or \code{FALSE}) is returned.
#'   Otherwise, a list with the following slots is returned. 
#'   \item{answer}{A boolean as above.}
#'   \item{H}{Hankel matrix of impulse response coefficients (with s block rows and s block columns).}
#'   \item{sv}{The singular values of H.}
#'   \item{s0}{(integer) (estimate of the) rank of H, i.e. the minimal statespace dimension.}
#'   
#' @export
#'
#' @examples
#' x = test_stsp(dim = c(2,2), s = 2)
#' is.minimal(x)
#' # note that operations on "stsp" objects may return non minimal realizations
#' # 
#' is.minimal(rbind(x, x), only.answer = FALSE)[c('answer','sv','s0')]
#' is.minimal(x %r% (x^(-1)), only.answer = FALSE)[c('answer','sv','s0')]
#' 
#' is.minimal(test_stsp(dim = c(2,0), s = 2))
#' is.minimal(test_stsp(dim = c(0,2), s = 0))
is.minimal = function(x, ...) {
  UseMethod("is.minimal", x)
}

#' @rdname is.minimal
#' @export
is.minimal.stsp = function(x, tol = sqrt(.Machine$double.eps), only.answer = TRUE, ...) {
  d = dim(x)
  s = unname(d[3])
  if (prod(d) == 0) {
    H = matrix(0, nrow = d[1]*s, ncol = d[2]*s)
    svH = numeric(0)
    s0 = 0
  } else {
    ir = unclass(pseries(x, lag.max = 2*s-1))[,,-1]
    H = bhankel(ir)
    svH = svd(H, nu = 0, nv = 0)$d
    s0 = sum(svH > tol)
  }
  is_minimal = (s == s0)
  
  if (only.answer) return(is_minimal)
  
  return(list(answer = is_minimal, H = H, sv = svH, s0 = s0))
}


#' State Transformation
#'
#' Applies a "state transformation" to a given rational matrix \eqn{K(z) = C(Iz^{-1} - A)^{-1}B +D} in statespace form. 
#' The parameter matrices are transformed as
#' \eqn{A \rightarrow T A T^{-1}}{A --> T A T^{-1}},
#' \eqn{B \rightarrow T B}{B --> T B} and
#' \eqn{C \rightarrow C T^{-1}}{C --> C T^{-1}}, where \eqn{T} is a non-singular (square) matrix. 
#'
#' @param obj \code{\link{stsp}} object, represents a rational matrix \eqn{K(z)}.
#' @param T state transformation matrix. The parameter \code{T} must be a non-singular (s-by-s) matrix, or a vector
#'        of length s^2 (in this case T is coerced to an s-by-s matrix) or a vector of length s. In the latter case
#'        T is coerced to a diagonal s-by-s matrix.
#' @param inverse if TRUE, the transformation is reversed, i.e. 
#' \eqn{A \rightarrow T^{-1} A T}{A --> T^{-1} A T},
#' \eqn{B \rightarrow T^{-1} B}{B --> T^{-1} B} and
#' \eqn{C \rightarrow C T}{C --> C T}.
#'
#' @return \code{\link{stsp}} object which represents the transformed state space realization.
#' @export
#'
#' @examples
#' obj = stsp(A = c(0,0.2,1,-0.5),
#'            B = c(1,1), C = c(1,0))
#'            
#' # random state transformation
#' T = stats::rnorm(4)
#' obj1 = state_trafo(obj, T)
#' 
#' # obj and obj1 are equivalent, they produce the same IRF
#' testthat::expect_equivalent(pseries(obj), pseries(obj1))
#' 
#' # diagonal state transformation matrix
#' T = stats::rnorm(2)
#' obj1 = state_trafo(obj, T)
#' 
#' # revert the transformation
#' testthat::expect_equivalent(obj, state_trafo(obj1, T, inverse = TRUE))
#' \dontrun{
#' state_trafo(obj, stats::rnorm(9)) # dimension of T does not fit
#' state_trafo(obj, c(1,0))          # T = diag(c(1,0)) is singular!
#' }
state_trafo = function(obj, T, inverse = FALSE) {
  if (!inherits(obj, 'stsp')) stop('argument "obj" must be "stsp" objcet')
  
  d = unname(dim(obj))
  m = d[1]
  n = d[2]
  s = d[3]
  
  if (s == 0) {
    if (length(T) == 0 ) return(obj)
    stop('The argument "T" is not compatible with "obj"')
  }
  
  if ( !(is.numeric(T) || is.complex(T)) ) stop('T must be a numeric (or complex) vector or a matrix!')
  if (is.vector(T)) {
    if (length(T) == s) {
      T = diag(x = T, nrow = s)
    } else {
      if (length(T) == (s^2)) {
        T = matrix(T, nrow = s, ncol = s)
      } else stop('T is not a compatible vector')
    }
  }
  if ( (!is.matrix(T)) || (ncol(T) != s) || (nrow(T) != s) ) stop('T must be a square, non-singular and compatible matrix!')
  
  obj = unclass(obj)
  
  if (inverse) {
    junk = try(solve(T, obj[1:s,,drop = FALSE]))
    if (inherits(junk, 'try-error')) stop('T is singular')
    obj[1:s,] = junk
    obj[,1:s] = obj[,1:s,drop = FALSE] %*% T
  } else {
    obj[1:s,] = T %*% obj[1:s,,drop = FALSE]
    junk = try(t(solve(t(T), t(obj[,1:s, drop = FALSE]))))
    if (inherits(junk, 'try-error')) stop('T is singular')
    obj[,1:s] = junk
  }
  obj = structure(obj, order = as.integer(c(m,n,s)),  class = c('stsp','ratm'))

  return(obj)
}

#' Grammians
#'
#' The procedure computes "grammians" of a statespace realization, which may e.g. be used for 
#' balancing the statespace realization. 
#' 
#' The \emph{controllability Grammian} \eqn{P} of a (stable) statespace realization 
#' \deqn{K(z) = C(Iz^{-1} - A)^{-1}B + D} 
#' is the solution of 
#' the Lyapunov equation \eqn{P = APA' + BB'}. The \emph{observability Grammian} is the solution 
#' of the Lyapunov equation \eqn{Q = A'QA + C'C}. If the statespace realization is \emph{stable} 
#' (the moduli of the eigenvalues of \code{A} are less than one) then \eqn{P,Q} are positive 
#' semidefinite 
#' and \eqn{P} is non singular if and only if the statespace realization is \emph{controllable} 
#' and \eqn{Q} is non singular if and only if the statespace realization is \emph{observable}. 
#' Hence the grammians may also be used to check whether the statespace realization is minimal 
#' (controllable \emph{and} observable). 
#' 
#' If the rational matrix is (strictly) minimum phase (i.e. \eqn{K(z)} is a square, invertible 
#' matrix and the eigenvalues of the matrix \eqn{(A - BD^{-1}C)} have moduli less than one) 
#' then we may also compute the controllability and the observability Grammian 
#' of the statespace realization 
#' \deqn{K^{-1}(z) = -D^{-1}C (Iz^{-1} - (A - BD^{-1}C))^{-1}BD^{-1} + D^{-1}.} 
#' of the inverse matrix \eqn{K^{-1}(z)}. These grammians have a similar interpretation. 
#' 
#' The above described grammians may be selected by setting the parameter \code{which} to 
#' \code{'ctr'}, \code{'obs'}, \code{'ctr_in'} or \code{'obs_inv'} respectively. 
#' 
#' For \emph{balancing} a statespace realization one needs a suitable pair of grammians. 
#' Two popular choices have been implemented: For \code{which = 'lyapunov'} the procedure 
#' returns the controllability and the observability Grammian and for \code{which = 'minimum phase'} 
#' the controllability matrix of the system and the observability Grammian of the inverse system 
#' are returned.
#' 
#' The procedure throws an error if the state space realization is not stable, respectively not 
#' minimum phase. 
#' 
#' @param obj (\code{\link{stsp}} object) rational matrix in statespace form.
#' @param which (character string) specifies the type of Grammian(s) to be computed.See 
#'        below for more details.
#'
#' @seealso \code{\link{ctr_matrix}}, \code{\link{obs_matrix}}, \code{\link{lyapunov}} 
#'          and \code{\link{balance}}. 
#' 
#' @return Either the selected Grammian (if \code{which} is one of \code{'ctr', 'obs', 'ctr_inv', 'obs_inv'}) 
#'         or a list with two components \code{P} and \code{Q} (for the case \code{which = 'lyapunov'} 
#'         or \code{which = 'miniumum phase'}).
#' @export
#'
#' @examples
#' # create a random, (3 by 2) rational matrix, 
#' # with a stable and minimum phase statespüace realization
#' obj = test_stsp(dim = c(3,2), s = 5, bpoles = 1, bzeroes = 1)
#' gr = grammians(obj, which = 'lyapunov')
#' gr
#' 
#' # we could also compute these grammians seperately 
#' all.equal(gr$P, grammians(obj,'ctr'))
#' all.equal(gr$Q, grammians(obj,'obs'))
#' 
#' # create a random (3 by 3) rational matrix, 
#' # with a stable and minimum phase statespüace realization
#' # Note: for the choice "minimum phase" the rational matrix 
#' # must be square and invertible.
#' obj = test_stsp(dim = c(3,3), s = 5, bpoles = 1, bzeroes = 1)
#' gr = grammians(obj, which = 'minimum phase')
#' gr
#' 
#' # we could also compute these grammians seperately 
#' all.equal(gr$P, grammians(obj,'ctr'))
#' all.equal(gr$Q, grammians(obj,'obs_inv'))
grammians = function(obj, which = c('lyapunov','minimum phase',
                                    'ctr','obs','obs_inv','ctr_inv')) {
  if (!inherits(obj, 'stsp')) stop('argument "obj" must a "stsp" object')
  
  which = match.arg(which)
  if (which == 'ctr') {
    P = try(lyapunov(obj$A, obj$B %*% t(obj$B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not stable')
    }
    return(P)
  }
  if (which == 'obs') {
    Q = try(lyapunov(t(obj$A), t(obj$C)%*% obj$C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      stop('statespace realization is not stable')
    }
    return(Q)
  }
  if (which == 'ctr_inv') {
    d = dim(obj)
    if (d[1] != d[2]) stop('the rational matrix must be square')
    
    # B matrix of inverse 
    B = t(solve(t(obj$D), t(obj$B)))
    if (inherits(B, 'try-error')) stop('obj is not invertible')
    # A matrix of inverse 
    A = obj$A - B %*% obj$C
    
    P = try(lyapunov(A, B %*% t(B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not minimum phase')
    }
    return(P)
  }
  if (which == 'obs_inv') {
    d = dim(obj)
    if (d[1] != d[2]) stop('the rational matrix must be square')
    
    # C matrix of inverse 
    C = -solve(obj$D, obj$C)
    if (inherits(C, 'try-error')) stop('obj is not invertible')
    # A matrix of inverse 
    A = obj$A + obj$B %*% C
    
    Q = try(lyapunov(t(A), t(C) %*% C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      stop('statespace realization is not minimum phase')
    }
    return(Q)
  }
  if (which == 'lyapunov') {
    P = try(lyapunov(obj$A, obj$B %*% t(obj$B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not stable')
    }
    # this is not efficient, since we compute the schur decomposition of A (above)
    # and then the schur decomposition of A' (below)
    Q = try(lyapunov(t(obj$A), t(obj$C) %*% obj$C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      # this should not happen
      stop('statespace realization is not stable')
    }
    return(list(P = P, Q = Q))
  }
  if (which == 'minimum phase') {
    d = dim(obj)
    if (d[1] != d[2]) stop('the rational matrix must be square')
    
    P = try(lyapunov(obj$A, obj$B %*% t(obj$B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not stable')
    }
    # C matrix of inverse 
    C = -solve(obj$D, obj$C)
    if (inherits(C, 'try-error')) stop('obj is not invertible')
    # A matrix of inverse 
    A = obj$A + obj$B %*% C
    
    Q = try(lyapunov(t(A), t(C) %*% C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      stop('statespace realization is not minimum phase')
    }
    return(list(P = P, Q = Q))
  }
  
  stop('this should not happen!')  
}


#' Balanced Realization and Balanced Truncation
#'
#' Compute a balanced realization or a balanced truncated statespace realization
#' of a rational matrix in statespace form.
#' 
#' Let \eqn{P,Q} denote the controllability and the observability Grammian of a statespace 
#' realization 
#' \deqn{K(z) = C(Iz^{-1} - A)^{-1}B +D} 
#' of a rational matrix \eqn{K(z)} with statespace dimension \eqn{s}. 
#' The matrix \eqn{PQ} is diagonalizable and the eigenvalues are real and non negative. 
#' Let \eqn{\sigma_1^2 \geq \cdots \geq \sigma_s^2 \geq 0}{\sigma[1]^2 \ge \dots \ge \sigma[s]^2 \ge 0}
#' denote the eigenvalues of \eqn{PQ} and suppose that \eqn{PQ} has rank \eqn{k\leq s}{k \le s}, 
#' i.e. \eqn{\sigma_i=0}{\sigma[i]=0} for \eqn{i>k}. 
#' Then there exists a statespace transformation \eqn{T} which renders both Grammians into diagonal form 
#' and furthermore we have \eqn{p_{ii}=q_{ii}=\sigma_i}{p[i,i]=q[i,i]=\sigma[i]} for 
#' the first \eqn{k} diagonal entries of the (transformed) Grammians.
#' 
#' The (non zero) \eqn{\sigma_i}{\sigma[i]}'s are equal to the (non zero) singular values of the 
#' Hankel matrix of the impulse response coefficients of \eqn{K} and hence are called 
#' \emph{Hankel singular values}. These singular values are also the singular values of the 
#' product \eqn{P^{1/2}Q^{1/2}} where \eqn{P^{1/2}} and \eqn{Q^{1/2}} denote "square roots" 
#' of the Grammians \eqn{P} and \eqn{Q} respectively.
#'  
#' The procedure \code{balance(obj,...)} computes a somewhat simplified "balanced" form where 
#' the transformed Grammians \eqn{P,Q} are block diagonal with two diagonal blocks of dimension 
#' \eqn{s_0}{s[0]} and \eqn{s_1 = s - s_0}{s[1] = s - s[0]} respectively. The two 
#' upper, left blocks of \eqn{P} and \eqn{Q} are diagonal and the diagonal entries are  
#' \eqn{p_{ii}=q_{ii}=\sigma_i}{p[i,i]=q[i,i]=\sigma[i]} for  \eqn{i=1,\ldots,s_0}{i=1,...,s[0]}.
#' If \eqn{s_0=k}{s[0]=k} is equal to the rank of \eqn{PQ} 
#' then the product of the two lower, right blocks of \eqn{P} and \eqn{Q} is zero 
#' (up to numerical errors). 
#' 
#' Note that \eqn{k=\mathrm{rk}(PQ)}{k = rk(PQ)} is equal to the minimal statespace dimension, i.e. 
#' there exists a statespace realization with statespace dimension \eqn{k} and 
#' any statespace realization of \eqn{K} has a statespace dimension \eqn{geq k}{\ge k}. 
#' Such a \emph{minimal} statespace realization now be constructed by "truncating" the balanced 
#' realization:  
#' \deqn{K(z) = C_{1}(Iz^{-1} - A_{11})^{-1} B_{1} + D}{K(z) = C[1](Iz^{-1} - A[1,1])^{-1} B[1] + D}
#' where \eqn{A_{11}}{A[1,1]} is the left, upper, \eqn{(k,k)} dimensional block of \eqn{A}, 
#' \eqn{B_1}{B[1]} is the upper, \eqn{(k,n)} dimensional block of \eqn{B} and 
#' \eqn{C_1}{C[1]} is the left, \eqn{(m,k)} dimensional block of \eqn{C}. 
#' 
#' This balanced truncated statespace realization is returned if the 
#' optional parameter \code{truncate} is \code{TRUE}. Note that in the case where 
#' \code{s0} is less than the rank of \eqn{PQ} the truncated realization 
#' is just an approximate realization of \eqn{K(z)}. The approximation error 
#' depends on the size of "neglected" Hankel singular values. Note also that the 
#' truncated realization is not balanced (if \eqn{s_0<k}{s[0]<k}). 
#' 
#' If the "target" statespace dimension \eqn{s_0}{s[0]} is not given then the 
#' procedure sets \eqn{s_0}{s[0]} equal to an estimate of the rank of 
#' \eqn{PQ}. This estimate is computed by inspecting the singular values 
#' \eqn{\sigma_i}{\sigma_[i]} of the product of the square roots 
#' \eqn{P^{1/2}} and \eqn{Q^{1/2}}. 
#'
#' The above discussion deals with the controllability and the observability 
#' Grammian of the statespace realization. However, one may also use 
#' other pairs of Grammians, e.g. for invertible matrices \eqn{K} 
#' one may use the controllabaility Grammian of the statespace realizatton 
#' and the observability Grammian of the statespace realization of the 
#' inverse \code{K^{-1}(z)}. 
#'
#' @param obj (\code{\link{stsp}} object) rational matrix in statespace form.
#' @param gr a list with two components \code{P} and \code{Q} which contain two Grammians. 
#'           These Grammians may e.g. be computed with \code{\link{grammians}}.
#' @param tol a tolerance parameter used for the determination of the rank \eqn{PQ}. 
#'            If \code{s0 = NULL} then the procedure estimates the rank of \eqn{PQ} 
#'            and sets \code{s0} equal to this estimate. If \code{s0} is given, 
#'            then \code{tol} is ignored. 
#' @param s0  determines the  size of the two diagonal blocks of transformed 
#'            Grammians \eqn{P}, \eqn{Q}, respectively (for \code{truncate = TRUE} 
#'            the statespace dimension of the balanced truncated statespace realization.
#' @param truncate (boolean) If true then the balanced truncated model is returned.
#'
#' @return
#' A list with components
#' \item{obj}{(\code{\link{stsp}} object) represents the balanced (truncated) statespace 
#'            realization.}
#' \item{T,Tinv}{the state transformation and the inverse state transformation matrix. 
#'               Note that for the case \code{truncate=TRUE} the matrix \code{T} is 
#'               an \code{s0 x s} matrix and \code{Tinv} is of dimension \code{s x s0}.}
#' \item{P,Q}{the transformed (and truncated) Grammians.}
#' \item{sigma}{the vector of singular values of the matrix \eqn{(P^{1/2}Q^{1/2})}.}
#' @export
#'
#' @examples
#' # example A ############################################################
#' 
#' # "obj" is a (1 by 1) rational matrix in statespace form, 
#' # with stespace dimension s = 2. 
#' 
#' obj = stsp(A = c(0,0.2,1,-0.5),
#'            B = c(1,1), C = c(1,0))
#' gr = grammians(obj, 'lyapunov')
#' bal = balance(obj, gr)
#'
#' print(cbind(bal$P, bal$Q, diag(bal$sigma, nrow = 2, ncol = 2)))
#' all.equal(grammians(bal$obj), bal[c('P','Q')])
#' 
#' # example B (non minimal statespace realization #########################
#' 
#' # The "rbind" operation below returns a statespace realization with 
#' # statespace dimension s = 4. However the minimal statespace dimensions 
#' # is s0 = 2. 
#' obj = rbind(obj, obj)
#' gr = grammians(obj, 'lyapunov')
#' bal = balance(obj, gr, s0 = 2, truncate = FALSE)
#'
#' # the upper (2 by 2) block of the (transformed) controllability 
#' # Grammian is diagonal, the lower (2 by 2) block is "zero". 
#' # This shows that the (balanced) realization is not controllable. 
#' print(bal$P)  
#'               
#' # the upper (2 by 2) block of the (transformed) observability 
#' # Grammian is diagonal and equal to the upper block of bal$P.
#' print(bal$Q) 
#' 
#' # the product of the (transformed) controllability and observability 
#' # Grammians is (approximately) diagonal and the diagonal entries are  
#' # the squares of the Hankel singular values.
#' print(bal$P %*% bal$Q)
#' print(bal$sigma^2)    
#' 
#' all.equal(grammians(bal$obj), bal[c('P','Q')])
#' 
#' # we may construct a minimal realization by 'balanced truncation'.
#' # note that we let the procedure determine the minimal statespace dimension 
#' trunc = balance(obj, gr)  
#' print(trunc$obj)
#' # compare with the above balanced realization 
#' print(bal$obj)
#' # check 
#' all.equal(pseries(obj), pseries(trunc$obj))
#' 
#' # example C (balanced truncation) ##########################
#' 
#' # construct a random rational matrix with statespace dimension s=10
#' obj = test_stsp(dim = c(2,2), s = 10, bpoles = 1, bzeroes = 1)
#' # compute an approximate realization with s0 = 8
#' gr = grammians(obj, 'minimum phase')
#' trunc = balance(obj, gr, s0 = 5)
#' print(trunc$sigma)
#' 
#' max(abs(unclass(pseries(obj, lag.max = 25)) - 
#'         unclass(pseries(trunc$ob, lag.max = 25))))
#' plot(pseries(obj, lag.max = 25), x_list= list(pseries(trunc$obj, lag.max = 25)), 
#'      type = c('l','p'), legend = c('s=10', 's=5')) 
balance = function(obj, gr, tol = 10*sqrt(.Machine$double.eps), s0 = NULL, truncate = TRUE) {
  # check inputs
  if (!inherits(obj,'stsp')) stop('obj must be an "stsp" object!')

  d = unname(dim(obj))
  m = d[1]
  n = d[2]
  s = d[3] # statespace dimension of the given statespace realization "obj"
  
  P = gr$P
  if ( (!is.numeric(P)) || (!is.matrix(P)) || (any(dim(P) != s)) ) {
    stop('argument "P" is not a compatible matrix')
  }
  Q = gr$Q
  if ( (!is.numeric(Q)) || (!is.matrix(Q)) || (any(dim(Q) != s)) ) {
    stop('argument "Q" is not a compatible matrix')
  }
  
  if (s == 0) {
    return(list(obj = obj, T = matrix(0, nrow = 0, ncol = 0), 
                Tinv = matrix(0, nrow = 0, ncol = 0), 
                P = P, Q = Q, sigma = numeric(0)))
  }

  # construct square roots of P and Q  
  # is it better to use eigen()?
  # with eigen we could check that P,Q are positive semidefinite

  out = svd(P, nv = 0)
  P2 = out$u %*% diag(x = sqrt(out$d), nrow = s, ncol = s) %*% t(out$u)
  # print(all.equal(P, P2 %*% P2))
  
  out = svd(Q, nv = 0)
  Q2 = out$u %*% diag(x = sqrt(out$d), nrow = s, ncol = s) %*% t(out$u)
  # print(all.equal(Q, Q2 %*% Q2))

  # svd of the product P2*Q2  
  out = svd(P2 %*% Q2)
  # (Hankel) singular values   
  sigma = out$d

  if (!is.null(s0)) {
    s0 = as.integer(s0)[1]
    if ((s0 < 0) || (s0 > s)) stop('illegal (target) statespace dimension "s0"')
  } else {
    s0 = sum(sigma > (tol * sigma[1]))
  }

  junk = matrix(1/sqrt(sigma[1:s0]), nrow = s0, ncol = s)  
  T = (t(out$v[, 1:s0, drop = FALSE]) %*% Q2) * junk
  S = (P2 %*% out$u[, 1:s0, drop = FALSE]) * t(junk)
  # print(all.equal(diag(s0), T %*% S))
  
  if ((truncate) || (s0 == s)) {
    # balance and truncate 
    if (s0 == 0) {
      return(list(obj = stsp(D = obj$D), T = matrix(0, nrow = 0, ncol = s0), 
                  Tinv = matrix(0, nrow = s0, ncol = 0), 
                  P = matrix(0, nrow = 0, ncol = 0), Q = matrix(0, nrow = 0, ncol = 0), 
                  sigma = sigma))        
    }

    obj = stsp(T %*% obj$A %*% S, T %*% obj$B, obj$C %*% S, obj$D) 
    return(list(obj = obj, T = T, Tinv = S, P = diag(x = sigma[1:s0], nrow = s0, ncol = s0), 
                Q = diag(sigma[1:s0], nrow = s0, ncol = s0), sigma = sigma))
  }
  
  # just balance 

  # extend T and S to square matrices
  out = svd(S %*% T)
  # print(out)
  T2 = t(out$u[, (s0+1):s, drop = FALSE])
  S2 = out$v[, (s0+1):s, drop = FALSE]
  # print(T %*% S2)
  # print(T2 %*% S)
  
  out = svd(T2 %*% S2)
  junk = matrix(1/sqrt(out$d), nrow = (s-s0), ncol = s)
  T2 = (t(out$u) %*% T2) * junk
  S2 = (S2 %*% out$v) * t(junk)
  T = rbind(T, T2)
  S = cbind(S, S2)
  # print( T %*% S)
  
  obj = stsp(T %*% obj$A %*% S, T %*% obj$B, obj$C %*% S, obj$D) 
  Pb = diag(x = sigma, nrow = s, ncol = s)
  Qb = Pb
  
  
  # print(P)
  # print(Pb)
  # print(Pb[(s0+1):s, (s0+1:s)])
  # print(T2)
  Pb[(s0+1):s, (s0+1):s] = T2 %*% P %*% t(T2)
  Qb[(s0+1):s, (s0+1):s] = t(S2) %*% Q %*% S2

  return(list(obj = obj, T = T, Tinv = S, P = Pb, Q = Qb, sigma = sigma))
}
