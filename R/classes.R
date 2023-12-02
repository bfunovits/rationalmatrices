#' Constructor for Laurent Polynomial Matrices
#'
#' \code{lpolm} objects represent Laurent polynomial matrices of the form
#' \deqn{a(z) = a_{-q} z^{-q} + \cdots + a_{-1} z^{-1} + a_0 + a_1 z + \cdots + a_p z^p}{
#'       a(z) = a[-q] (1/z)^q + \cdots + a[-1] 1/z + a[0] + a[1]z + \dots + a[p]z^p.}
#' If \eqn{a(z)} is an \eqn{(m,n)} dimensional Laurent polynomial matrix
#' (i.e. the coefficients \eqn{a_i}{a[i]} are \eqn{(m,n)} dimensional real or complex valued matrices),
#' then the \code{lpolm} object stores the coefficients in 
#' an \code{(m,n,q+p+1)}-dimensional (numeric or complex) array \eqn{a_{-q}, \ldots, a_{-1}, a_{0}, a_{1}, \ldots, a_{p}},
#' together with a class attribute \code{c("lpolm","ratm")}. 
#' \cr
#' The constructor function \code{lpolm(a, min_deg)} takes an integer \code{min_deg} (default is zero) and 
#' a (numeric or complex) vector, matrix or 3-dimensional array and returns an \code{lpolm} object.
#' 
#' Any of the dimensions of the 3-dimensional array may also be zero. 
#' In particular, if the third dimension is zero, then the \code{lpolm} object is interpreted as the zero Laurent polynomial.
#' 
#' This class is special in the sense that it is not possible to \emph{upgrade} it to an \code{\link{lmfd}}, \code{\link{rmfd}}, \code{\link{stsp}}, or \code{\link{pseries}} object.
#' 
#' For important methods and functions for this class have a look at the "see also" section.
#' Note that some functions are only written for \code{\link{polm}} objects: \code{\link{degree}}, \code{\link{col_end_matrix}}, normal forms like \code{\link{snf}}, \code{\link{hnf}}, \code{\link{whf}}.
#' 
#' @param a either a (numeric or complex) vector, matrix or 3-D array. A vector is coerced to a scalar 
#'          (i.e. \eqn{(1,1)}-dimensional) polynomial and a matrix gives a Laurent polynomial matrix whose only coefficient matrix is of degree \code{min_deg}. 
#' @param min_deg Integer. Default set to zero. Smallest degree in the Laurent polynomial.
#'                Negative for Laurent polynomials which cannot be coerced to \code{\link{polm}} objects (\code{\link{as.polm}} throws an error).
#'
#' @return An object of class \code{c("lpolm", "ratm")}. 
#'         When non-negative, the object is still of this class. 
#'         However, it can then be coerced to an object of class \code{c("polm", "ratm")} using \code{\link{as.polm}}
#' @export
#' 
#' @seealso 
#' \itemize{
#'    \item \code{\link{polm}}.
#'    \item \code{\link{test_lpolm}} generates random polynomials. 
#'    \item checks: \code{\link{is.lpolm}}
#'    \item generic S3 methods: \code{\link[rationalmatrices]{dim}}, 
#'          \code{\link[rationalmatrices]{str}} and \code{\link[rationalmatrices]{print}}.
#'    \item arithmetics: \code{\link[rationalmatrices]{Ops.ratm}}, matrix multiplication \code{\link{\%r\%}}
#'    \item matrix operations: \code{\link[rationalmatrices]{t.lpolm}}, 
#'          \code{\link[rationalmatrices]{bind}}, \code{\link{[.lpolm}}, \code{\link{[<-.lpolm}}, ...
#'    \item \code{\link{prune}} discards leading and trailing zero matrices of Laurent polynomials. 
#'    \item \code{\link{get_bwd}} discards all coefficient matrices pertaining to \strong{negative} powers and returns a polm object.
#'    \item \code{\link{get_fwd}} discards all coefficient matrices pertaining to \strong{positive} powers. It returns thus an lpolm object.
#'    \item \code{\link{polm2fwd}} performs the transformation \eqn{p(z) \rightarrow p(z^{-1})} on a polm object and returns an lpolm object.
#' }
#' 
#' @examples 
#' # (1 x 1) Laurent polynomial matrix a(z) =  3z^{-2} + 2z^{-1} + 1
#' lpolm(3:1, min_deg = -2)
#' 
#' # Non-negative minimal degrees are allowed too (no implicit coercion to polm object)
#' lpolm(3:1, min_deg = 2)
#' 
#' lpolm(matrix(1:4,2,2), min_deg = -2)
#' 
lpolm = function(a, min_deg = 0) {
  
  stopifnot("Only array, matrix, or vector objects can be supplied to this constructor." = any(class(a) %in% c("matrix", "array")) || is.vector(a))
  
  stopifnot('Input "a" must be a numeric or complex vector/matrix/array!' = (is.numeric(a) || is.complex(a)))
  
  if (is.vector(a)) {
    dim(a) = c(1,1,length(a))
  }
  if (is.matrix(a)) {
    dim(a) = c(dim(a),1)
  }
  
  stopifnot('could not coerce input parameter "a" to a valid 3-D array!' = is.array(a) && (length(dim(a)) == 3),
            "Minimal degree 'min_deg' must be set, numeric, and of length 1!" = !is.null(min_deg) && is.numeric(min_deg) && length(min_deg) == 1)
  
  # achtung
  min_deg = as.integer(min_deg)[1]

  # Initially, lpolm object were implicitly coerced to polm object when deg_min >= 0. 
  # This was not optimal in terms of type consistency. 
  # It is now necessary to coerce such lpolm objects to polm objects explicitly with as.polm.lpolm()
  # The code below is used in as.polm.lpolm()
  # if (min_deg >= 0) {
  #   polm_offset = array(0, dim = c(dim(a)[1:2], min_deg))
  #   return(polm(dbind(d = 3, polm_offset, a)))
  # }
  
  return(structure(a, class = c("lpolm", 'ratm'), min_deg = min_deg)) 
}


#' Constructor for Polynomial Matrices
#'
#' \code{polm} objects represent polynomial matrices
#' \deqn{a(z) = a_0 + a_1 z + \cdots + a_p z^p}{a(z) = a[0] + a[1]z + \dots + a[p]z^p}
#' If the matrix \eqn{a(z)} is an \eqn{(m,n)} dimensional polynomial matrix (i.e. the coefficients 
#' \eqn{a_i}{a[i]} are \eqn{(m,n)} dimensional real or complex valued matrices) then the 
#' \code{polm} object stores the coefficients in an \code{(m,n,p+1)} dimensional (numeric or complex) 
#' array together with a class attribute \code{c("polm","ratm")}. 
#' \cr
#' The constructor function \code{polm(a)} takes a (numeric or complex) vector, matrix or 
#' 3-dimensional array and returns a \code{polm} object.
#' 
#' Any of the dimensions of the 3-dimensional array may also be zero. In particular, 
#' if the third dimension is zero, then the \code{polm} object is interpreted as the zero polynomial.
#' 
#' For important methods and functions for this class have a look at the "see also" section.
#' 
#' @param a either a (numeric or complex) vector, matrix or 3-D array. A vector is coerced to a scalar 
#'          (i.e. \eqn{(1,1)}-dimensional) polynomial and a matrix gives a polynomial matrix of zero degree. 
#'
#' @return An object of class \code{polm}. 
#' @export
#' 
#' @seealso 
#' \itemize{
#'    \item \code{\link{lpolm}} objects allow for coefficient matrices pertaining to negative powers of \eqn{z}.
#'    \item \code{\link{test_polm}} generates random polynomials. 
#'    \item checks: \code{\link{is.polm}}, \code{\link{is.miniphase}}, 
#'          and \code{\link{is.coprime}}. As a byproduct \code{is.coprime(a)} 
#'          computes the zeroes of a square polynomial matrix \eqn{a(z)}. 
#'    \item generic S3 methods: \code{\link[rationalmatrices]{dim}}, 
#'          \code{\link[rationalmatrices]{str}} and \code{\link[rationalmatrices]{print}}.
#'    \item arithmetics: \code{\link[rationalmatrices]{Ops.ratm}}, matrix multiplication \code{\link{\%r\%}}, 
#'          polynomial division \code{\%/\%}, 
#'          polynomial remainder \code{\%\%}, ...
#'    \item matrix operations: \code{\link[rationalmatrices]{t.polm}}, 
#'          \code{\link[rationalmatrices]{bind}}, \code{\link{[.polm}}, \code{\link{[<-.polm}}, ...
#'    \item \code{\link{degree}} returns the degree, \code{\link{col_end_matrix}} computes the 
#'          \emph{column end matrix} and \code{\link{prune}} "simplifies" a polynomial.
#'          Note that the degree of the zero polynomial is implemented as being equal to \eqn{-1}, see \code{\link{degree}}!
#'    \item \code{\link{reflect_zeroes}} may be used 
#'          to reflect zeroes of a polynomial matrix by multiplication
#'          with allpass rational matrices.
#'    \item normal forms: Hermite normal form \code{\link{hnf}}, Smith normal form \code{\link{snf}},  
#'                        column reduced form \code{\link{col_reduce}} and Wiener Hopf factorization 
#'                        \code{\link{whf}}.
#'    \item \code{\link{companion_matrix}}, \code{\link{zeroes}}, \code{\link{pseries}}, \code{\link{zvalues}}, ...  
#' }
#' 
#' @examples 
#' # (1 x 1) polynomial matrix a(z) =  0 + 1z + 2z^2
#' polm(0:2)
#' 
#' # (2 x 3) polynomial matrix a(z) = a0 (degree is zero)
#' polm(diag(1, nrow = 2, ncol = 3))
#' 
#' # random (2 x 3) polynomial matrix a(z) = a0 + a1 z + a2 z^2 + a3 z^3 (degree = 3)
#' polm(array(stats::rnorm(2*3*4), dim = c(2,3,4)))
#' 
#' # random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
#' a = polm(array(complex(real = stats::rnorm(2*1*3), 
#'                imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
#' is.polm(a)
#' dim(a)
#' str(a)
#' print(a, digits = 3)
polm = function(a) {
  
  stopifnot("Only array, matrix, or vector objects can be supplied to this constructor." = any(class(a) %in% c("matrix", "array")) || is.vector(a))
  
  if (!(is.numeric(a) || is.complex(a))) {
    stop('input "a" must be a numeric or complex vector/matrix/array!')
  }
  if (is.vector(a)) {
    dim(a) = c(1,1,length(a))
  }
  if (is.matrix(a)) {
    dim(a) = c(dim(a),1)
  }
  
  if ( (!is.array(a)) || (length(dim(a)) !=3) ) {
    stop('could not coerce input parameter "a" to a valid 3-D array!')
  }
  class(a) = c('polm','ratm')
  return(a)
}

#' Constructor for Left Matrix Fraction Descriptions (LMFDs)
#'
#' A Left Matrix Fraction Description (LMFD) of a rational matrix, \eqn{x(z)} say, is a
#' pair \eqn{(a(z),b(z))} of polynomial matrices, such
#' that \eqn{x(z) = a^{-1}(z) b(z)}. 
#' The polynomial matrix \eqn{a(z)} must be square and invertible.
#' 
#' Suppose that \eqn{x(z)=a^{-1}(z) b(z)} is an \eqn{(m,n)}-dimensional matrix and that 
#' \eqn{a(z)} and \eqn{b(z)} have degrees \eqn{p} and \eqn{q} respectively. 
#' The corresponding \code{lmfd} object stores the coefficients of the polynomials \eqn{a(z), b(z)} in 
#' an \eqn{(m,m(p+1)+n(q+1))} dimensional (real or complex valued) matrix together with an 
#' attribute \code{order = c(m,n,p,q)} and a class attribute \code{c("lmfd", "ratm")}. 
#' 
#' For a valid LMFD we require \eqn{m>0} and \eqn{p\geq 0}{p \ge 0}. 
#'
#' @param a,b \code{\link{polm}} objects, or objects which may be coerced to a \code{\link{polm}} object,
#'            via \code{x = polm(x)}. Either of the two arguments may be omitted.
#'
#'
#' @seealso Useful methods and functions for the \code{lmfd} class are: 
#' \itemize{
#'    \item \code{\link{test_lmfd}} generates random rational matrices in LMFD form.
#'    \item checks: \code{\link{is.lmfd}}, \code{\link{is.miniphase}}, \code{\link{is.stable}} 
#'          and \code{\link{is.coprime}}.
#'    \item generic S3 methods: \code{\link[rationalmatrices]{dim}}, 
#'          \code{\link[rationalmatrices]{str}} and \code{\link[rationalmatrices]{print}}.
#'    \item arithmetics: \code{\link[rationalmatrices]{Ops.ratm}}.
#'    \item matrix operations: \code{\link[rationalmatrices]{bind}}.
#'    \item extract the factors \eqn{a(z)} and \eqn{b(z)} with \code{\link{$.lmfd}}.      
#'    \item \code{\link{zeroes}}, \code{\link{pseries}}, \code{\link{zvalues}}, ...  
#' }
#' 
#' @return An object of class \code{lmfd}.
#' @export
#' 
#' @examples 
#' ### (1 x 1) rational matrix x(z) = (1+z+z^2)^(-1) (3+2z+z^2)
#' lmfd(c(1,1,1), c(3,2,1)) %>% print(format = 'c')
#' 
#' ### (1 x 1) rational matrix x(z) = (3+2z+z^2)
#' lmfd(b = c(3,2,1)) %>% print(format = 'c')
#' 
#' ### (1 x 1) rational matrix x(z) = (1+z+z^2)^(-1)
#' lmfd(c(1,1,1)) %>% print(format = 'c')
#' 
#' ### (2 x 3) rational matrix with degrees p=1, q=1
#' x = lmfd(array(rnorm(2*2*2), dim = c(2,2,2)), 
#'          array(rnorm(2*3*2), dim = c(2,3,2)))
#' is.lmfd(x)
#' dim(x)
#' str(x)
#' print(x, digits = 2)
#' 
#' \dontrun{
#' ### the following calls to lmfd() throw an error 
#' lmfd() # no arguments!
#' lmfd(a = test_polm(dim = c(2,3), degree = 1))  # a(z) must be square 
#' lmfd(a = test_polm(dim = c(2,2), degree = -1)) # a(z) must have degree >= 0
#' }
lmfd = function(a, b) {
  if (missing(a) && missing(b)) {
    stop('no arguments have been provided')
  }
  if (!missing(a)) {
    if (!inherits(a,'polm')) {
      a = try(polm(a))
      if (inherits(a, 'try-error')) {
        stop('could not coerce "a" to a "polm" object!')
      }
    }
    dim_a = dim(unclass(a))
    if ((dim_a[1] == 0) || (dim_a[1] != dim_a[2]) || (dim_a[3] == 0)) {
      stop('"a" must represent a square polynomial matrix with degree >= 0')
    } 
  }
  if (!missing(b)) {
    if (!inherits(b,'polm')) {
      b = try(polm(b))
      if (inherits(b, 'try-error')) {
        stop('could not coerce "b" to a "polm" object!')
      }
    }
    dim_b = dim(unclass(b))
  }
  if (missing(b)) {
    b = polm(diag(dim_a[2]))
    dim_b = dim(unclass(b))
  }
  if (missing(a)) {
    a = polm(diag(dim_b[1]))
    dim_a = dim(unclass(a))
  }
  if (dim_a[2] != dim_b[1]) {
    stop('the arguments "a", "b" are not compatible')
  }
  
  c = matrix(0, nrow = dim_b[1], ncol = dim_a[2]*dim_a[3] + dim_b[2]*dim_b[3])
  if ((dim_a[2]*dim_a[3]) > 0) c[, 1:(dim_a[2]*dim_a[3])] = a
  if ((dim_b[2]*dim_b[3]) > 0) c[, (dim_a[2]*dim_a[3]+1):(dim_a[2]*dim_a[3] + dim_b[2]*dim_b[3])] = b
  
  c = structure(c, order = as.integer(c(dim_b[1], dim_b[2], dim_a[3]-1, dim_b[3] - 1)),  
                class = c('lmfd','ratm'))
  return(c)
}


#' Constructor for Right Matrix Fraction Descriptions (RMFDs)
#'
#' A Right Matrix Fraction Description (RMFD) of a rational matrix, \eqn{k(z)}
#' say, is  pair \eqn{(c(z),d(z))} of polynomial matrices, such that
#' \eqn{k(z) = d(z) c^{-1}(z)}. The polynomial matrix \eqn{c(z)} must be square
#' and invertible.
#'
#' Suppose that \eqn{k(z) = d(z) c^{-1}(z)} is an \eqn{(m,n)}-dimensional matrix
#' and that \eqn{c(z)} and \eqn{d(z)} have degrees \eqn{p} and \eqn{q}
#' respectively. The corresponding \code{rmfd} object stores the coefficients of
#' the polynomials \eqn{c(z), d(z)} in an \eqn{(n(p+1)+m(q+1), n)} dimensional
#' (real or complex valued) matrix. The \code{rmfd} object also stores the
#' attribute \code{order = c(m,n,p,q)} and a class attribute \code{c("lmfd",
#' "ratm")}.
#' 
#' For a valid RMFD we require \eqn{m>0} and \eqn{p\geq 0}{p \ge 0}. 
#'
#' @param c,d \code{\link{polm}} objects, or objects which may be coerced to a \code{\link{polm}} object, via \code{x = polm(x)}. 
#'   Either of the two arguments may be omitted.
#'
#' @seealso Useful methods and functions for the \code{rmfd} class are: 
#' \itemize{
#'    \item \code{\link{test_rmfd}} generates random rational matrices in RMFD form.
#'    \item checks: \code{\link{is.rmfd}}, \code{\link{is.miniphase}}, \code{\link{is.stable}} 
#'          and \code{\link{is.coprime}}.
#'    \item generic S3 methods: \code{\link[rationalmatrices]{dim}}, 
#'          \code{\link[rationalmatrices]{str}} and \code{\link[rationalmatrices]{print}}.
#'    \item arithmetics: \code{\link[rationalmatrices]{Ops.ratm}}.
#'    \item matrix operations: \code{\link[rationalmatrices]{bind}}.
#'    \item extract the factors \eqn{c(z)} and \eqn{d(z)} with \code{\link{$.rmfd}}.      
#'    \item \code{\link{zeroes}}, \code{\link{pseries}}, \code{\link{zvalues}}, ...  
#' }
#' 
#' @return An object of class \code{rmfd}.
#' @export
#' 
#' @examples 
#' ### (1 x 1) rational matrix k(z) =  (3+2z+z^2) (1+z+z^2)^(-1)
#' rmfd(c(1,1,1), c(3,2,1)) %>% print(format = 'c')
#' 
#' ### (1 x 1) rational matrix k(z) = (3+2z+z^2)
#' rmfd(c = c(3,2,1)) %>% print(format = 'c')
#' 
#' ### (1 x 1) rational matrix k(z) = (1+z+z^2)^(-1)
#' rmfd(c(1,1,1)) %>% print(format = 'c')
#' 
#' ### (3 x 2) rational matrix with degrees p=1, q=1
#' x = rmfd(array(rnorm(2*2*2), dim = c(2,2,2)), 
#'          array(rnorm(3*2*2), dim = c(3,2,2)))
#' is.rmfd(x)
#' dim(x)
#' str(x)
#' print(x, digits = 2)
#' 
#' xr = test_rmfd(dim = c(3,2), degree = c(2,2))
#' 
#' \dontrun{
#' ### the following calls to rmfd() throw an error 
#' rmfd() # no arguments!
#' rmfd(c = test_polm(dim = c(2,3), degree = 1))  # c(z) must be square 
#' rmfd(c = test_polm(dim = c(2,2), degree = -1)) # c(z) must have degree >= 0
#' }
rmfd = function(c = NULL, d = NULL) {
  
  # Check inputs: At least one polm object must be non-empty ####
  
  if (is.null(c) && is.null(d)) {
    stop('At least one of c(z) or d(z) needs to be provided.')
  }
  
  # Convert array to polm objects ####
  if (!is.null(c)) {
    if (!inherits(c,'polm')) {
      c = try(polm(c))
      if (inherits(c, 'try-error')) {
        stop('Could not coerce "c" to a "polm" object!')
      }
    }
    dim_c = dim(unclass(c))
    if ((dim_c[1] == 0) || (dim_c[1] != dim_c[2]) || (dim_c[3] == 0)) {
      stop('"c" must represent a square polynomial matrix with degree >= 0')
    } 
  }
  
  if (!is.null(d)) {
    if (!inherits(d,'polm')) {
      d = try(polm(d))
      if (inherits(d, 'try-error')) {
        stop('Could not coerce "d" to a "polm" object!')
      }
    }
    dim_d = dim(unclass(d))
  }
  
  # If one polm object is NULL, set to identity matrix ####
  
  if (is.null(d)) {
    d = polm(diag(dim_c[2]))
    dim_d = dim(unclass(d))
  }
  
  if (is.null(c)) {
    c = polm(diag(dim_d[2]))
    dim_c = dim(unclass(c))
  }
  
  # Check input: Dimensions ####
  
  if (dim_c[1] != dim_d[2]) {
    stop('The dimensions of "c" and "d" are not compatible.')
  }
  
  h = matrix(0, 
             nrow = dim_c[2], # number of cols! 
             ncol = dim_c[1]*dim_c[3] + dim_d[1]*dim_d[3]) # (cols of c(z))* (lag.max + 1) + (cols of d(z))* (lag.max + 1) 
  if ((dim_c[1]*dim_c[3]) > 0) h[, 1:(dim_c[1]*dim_c[3])] = aperm(c, c(2,1,3))
  if ((dim_d[1]*dim_d[3]) > 0) h[, (dim_c[1]*dim_c[3]+1):(dim_c[1]*dim_c[3] + dim_d[1]*dim_d[3])] = aperm(d, c(2,1,3))
  h = t(h)
  
  h = structure(h, 
                order = as.integer(c(dim_d[1], dim_d[2], dim_c[3]-1, dim_d[3] - 1)),  
                class = c('rmfd','ratm'))
  return(h)
}

#' Constructor for Statespace Realizations
#' 
#' Any rational (\eqn{(m,n)}-dimensional) matrix \eqn{x(z)} which has no pole at 
#' \eqn{z=0} may be represented as 
#' \deqn{x(z) = C(z^{-1} I_s - A)^{-1} B + D}{x(z) = C(z^{-1} I - A)^{-1} B + D}
#' Here \eqn{I_s}{I} denotes the \eqn{(s,s)}-dimensional identity matrix and \eqn{A,B,C,D} 
#' are (real or complex valued) matrices of size \eqn{(s,s)}, \eqn{(s,n)}, \eqn{(m,s)} and \eqn{(m,n)}
#' respectively. The integer \eqn{s \geq 0}{s \ge 0} is called the \emph{state dimension} 
#' of the above \emph{statespace realization} of \eqn{x(z)}. 
#' 
#' Internally statespace realizations are stored as an \eqn{(s+m, s+n)}-dimensional matrix 
#' with an attribute \code{order = c(m,n,s)} and a class attribute 
#' \code{c("stsp","ratm")}.
#' 
#' Any of the integers \eqn{m,n,s} may be zero. 
#'
#' If the arguments \code{A,B,C} are missing, then a statespace realization with 
#' statespace dimension \eqn{s=0} is constructed. In this case \code{D} must a be matrix or 
#' a scalar. 
#' If \code{A,B,C} are given (and compatible) and \code{D} is missing, then 
#' \code{D = diag(x = 1, nrow = m, ncol = n)} is used.
#' 
#' @param A \eqn{(s,s)} matrix (or a vector of length \eqn{(s^2)})
#' @param B \eqn{(s,n)} matrix (or a vector of length \eqn{(sn)})
#' @param C \eqn{(m,s)} matrix (or a vector of length \eqn{(ms)})
#' @param D \eqn{(m,n)} matrix (or a vector of length \eqn{(mn)})
#'
#' @return An object of class \code{stsp}. 
#' @export
#' 
#' @seealso 
#' \itemize{
#'    \item \code{\link{test_stsp}} generates random polynomials. 
#'    \item checks: \code{\link{is.stsp}}, \code{\link{is.miniphase}}, 
#'          \code{\link{is.stable}} and \code{\link{is.minimal}}.
#'    \item generic S3 methods: \code{\link[rationalmatrices]{dim}}, 
#'          \code{\link[rationalmatrices]{str}} and \code{\link[rationalmatrices]{print}}.
#'    \item arithmetics: \code{\link[rationalmatrices]{Ops.ratm}}, matrix multiplication \code{\link{\%r\%}}, ... 
#'    \item matrix operations: transpose \code{\link[rationalmatrices]{t.stsp}}, Hermitean transpose 
#'          \code{\link[rationalmatrices]{Ht}}, \code{\link[rationalmatrices]{bind}}, 
#'          \code{\link{[.stsp}}, ...
#'    \item extract the system matrices \eqn{A,B,C,D} with \code{\link{$.stsp}}.      
#'    \item statespace realizations related tools: state transfromation \code{\link{state_trafo}}, 
#'          observability and controllability matrices \code{\link{obs_matrix}} and 
#'          \code{\link{ctr_matrix}},  Grammian matrices \code{\link{grammians}} and 
#'          balanced realizations \code{\link{balance}}.. 
#'    \item \code{\link{reflect_poles}} and \code{\link{reflect_zeroes}} may be used 
#'          to reflect poles and zeroes of a rational matrix in statespace form by multiplication
#'          with allpass rational matrices.
#'    \item \code{\link{zeroes}}, \code{\link{pseries}}, \code{\link{zvalues}}, ...  
#' }
#' 
#' @examples 
#' ### x(z) =  I
#' stsp(D = diag(3))
#' 
#' ### random (2 x 3) rational matrix in statespace form with state dimension s = 4
#' x = stsp(A = stats::rnorm(4*4), B = stats::rnorm(4*3), C = stats::rnorm(2*4))
#' 
#' is.stsp(x)
#' dim(x)
#' str(x)
#' print(x, digits = 3)
stsp = function(A, B, C, D) {
  # only D has been given => statespace dimension s = 0
  if (missing(A) && missing(B) && missing(C)) {
    if (missing(D)) stop('no parameters')
    if (!( is.numeric(D) || is.complex(D) )) stop('parameter D is not numeric or complex')
    if ( (is.vector(D)) && (length(D) == 1) ) {
      D = matrix(D, nrow = 1, ncol = 1)
    }
    if (!is.matrix(D)) {
      stop('parameter D must be a numeric or complex matrix')
    }
    m = nrow(D)
    n = ncol(D)
    s = 0
    ABCD = structure(D, order = as.integer(c(m,n,s)),  class = c('stsp','ratm'))
    return(ABCD)
  }
  
  if (missing(A)) stop('parameter A is missing')
  if ( !( is.numeric(A) || is.complex(A) ) ) stop('parameter A is not numeric or complex')
  if ( is.vector(A) && (length(A) == (as.integer(sqrt(length(A))))^2) )  {
    A = matrix(A, nrow = sqrt(length(A)), ncol = sqrt(length(A)))
  }
  if ( (!is.matrix(A)) || (nrow(A) != ncol(A)) ) stop('parameter A is not a square matrix')
  s = nrow(A)
  
  if (missing(B)) stop('parameter B is missing')
  if (!( is.numeric(B) || is.complex(B) )) stop('parameter B is not numeric or complex')
  if ( is.vector(B) &&  (s > 0) && ((length(B) %% s) == 0) ) {
    B = matrix(B, nrow = s)
  }
  if ( (!is.matrix(B)) || (nrow(B) != s) ) stop('parameter B is not compatible to A')
  n = ncol(B)
  
  if (missing(C)) stop('parameter C is missing')
  if (!( is.numeric(C) || is.complex(C) )) stop('parameter C is not numeric or complex')
  if ( is.vector(C) &&  (s > 0) && ((length(C) %% s) == 0) ) {
    C = matrix(C, ncol = s)
  }
  if ( (!is.matrix(C)) || (ncol(C) != s)) stop('parameter C is not compatible to A')
  m = nrow(C)
  
  if (missing(D)) D = diag(x = 1, nrow = m, ncol = n)
  if (!( is.numeric(D) || is.complex(D) )) stop('parameter D is not numeric or complex')
  if ( is.vector(D) && (length(D) == (m*n)) ) {
    D = matrix(D, nrow = m, ncol = n)
  }
  if ( (!is.matrix(D)) || (nrow(D) != m)  || (ncol(D) != n) ) {
    stop('parameter D is not compatible to B,C')
  }
  
  ABCD = rbind( cbind(A,B), cbind(C, D) )
  ABCD = structure(ABCD, order = as.integer(c(m,n,s)),  class = c('stsp','ratm'))
  return(ABCD)
}
