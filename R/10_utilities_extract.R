# ===================================================================
# Functional Area #10: Utilities & Helpers (CLAUDE.md)
#
# Purpose: Implement subsetting and extraction methods
#   - []: Subsetting and element extraction
#   - $: Component extraction
#   - Support for all rational matrix representations
#
# Related Files: 10_utilities_dim.R (dimensions), 10_utilities_print.R (display)
# ===================================================================

# subsetting for rational matrices ###############################################

# helper function:
# check the result of subsetting x[]
# based on the result of subsetting an "ordinary" matrix 
#' Helper for Extracting indices 
#' 
#' Linear indices for extracting elements of a matrix are obtained, 
#' depending on the dimension of the matrix, number of arguments \code{n_args}, the missingness \code{is_missing} of these arguments, and their contents \code{i} and \code{j}.
#' 
#' No input checks are performed because this is an internal function. 
#' Consequently, some input combination are contradictory without alerting the user.
#'
#' @param m Integer. Output dimension.
#' @param n Integer. Input dimension.
#' @param n_args Integer. One, or two. Corresponds to the case where one set of indices or two sets are given.
#'               Corresponds to x[i] and x[i,j]
#' @param is_missing Boolean vector of dimension \code{n_args}, i.e. one or two. 
#'                   TRUE when this function is used fo x[] and x[,]
#' @param i First set of indices. Can be integers or booleans.
#' @param j Second set of indices. Like above.
#'
#' @return Indices (integers) which will be further used to extract the relevant elements
#' @keywords internal
#' @export
#'
#' @examples
#' # x[,]
#' extract_matrix_(6, 5, 
#'                 n_args = 2,
#'                 is_missing = c(TRUE, TRUE),
#'                 i = 1:3,
#'                 j = c(2,5))
#'                 
#' # x[i,j]
#' extract_matrix_(6, 5, 
#'                 n_args = 2,
#'                 is_missing = c(FALSE, FALSE),
#'                 i = 1:3,
#'                 j = c(2,5))
#'                 
#' # x[,j]
#' extract_matrix_(6, 5, 
#'                 n_args = 2,
#'                 is_missing = c(TRUE, FALSE),
#'                 i = 1:3,
#'                 j = c(2,5))
#'                 
#' # x[i,]
#' extract_matrix_(6, 5, 
#'                 n_args = 2,
#'                 is_missing = c(FALSE, TRUE),
#'                 i = 1:3,
#'                 j = c(2,5))
#'                 
#' # x[i], j is ignored if available (doesn't happen because it will not be called in this way)
#' extract_matrix_(6, 5, 
#'                 n_args = 1,
#'                 is_missing = c(FALSE, TRUE),
#'                 i = 1:10,
#'                 j = c(2,5))
#'                 
#' # x[i], j is ignored if available (doesn't happen because it will not be called in this way)
#' extract_matrix_(6, 5, 
#'                 n_args = 1,
#'                 is_missing = c(FALSE, TRUE),
#'                 i = c(3,6,7,1),
#'                 j = c(2,5))
extract_matrix_ = function(m, n, n_args, is_missing, i, j) {
  # Write linear indices
  idx = matrix(iseq(1, m*n), nrow = m, ncol = n)
  
  # Case 1: Everything missing: x[] and x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(idx)
  
  # Case 2: One argument: x[i] (input j is ignored if available (which will not be the case))
  if (n_args == 1) {
    # x[i]
    idx = idx[i]
    if (any(is.na(idx))) stop('index out of bounds')
    idx = matrix(idx, nrow = length(idx), ncol = 1)
    return(idx)
  }
  
  # Case 3a: Two arguments, first missing: x[,j]
  if (is_missing[1]) {
    # x[,j]
    return(idx[, j, drop = FALSE])
  }
  
  # Case 3b: Two arguments, second missing: x[i,]
  if (is_missing[2]) {
    # x[i,]
    return(idx[i, , drop = FALSE])
  }
  
  # Case 3c: Two arguments, none missing: x[i,j]
  return(idx[i, j, drop = FALSE])
}


#' Replace Parts of a (Laurent) Polynomial Matrix
#' 
#' The assignment operation \code{x[,] <- value} for (Laurent) polynomial matrices works quite analogously 
#' to the assignment operation of "ordinary" matrices. 
#' \cr
#' In the case of Laurent polynomial objects \code{\link{lpolm}} and if \code{value} is not an \code{\link{lpolm}} object itself 
#' (i.e. \code{value} is a vector, matrix, or array), \code{value} will first be coerced to an \code{\link{lpolm}} object with \code{min_deg = 0}.
#' 
#' Note that "named" arguments are not supported (in order to simplify the coding). 
#'
#' @param x \code{\link{polm}} or \code{\link{lpolm}} object
#' @param i,j indices
#' @param value  Either a \code{\link{polm}} or \code{\link{lpolm}} object, or a vector/matrix/array which may be coerced to 
#'               a \code{polm} or \code{\link{lpolm}} object by \code{polm(value)} or \code{lpolm(value)}.
#'
#' @export
#' @rdname replace
#' @name replace
#'
#' @examples
#' a = test_polm(dim = c(3,2), degree = 1)
#' print(a, format = 'c')
#' 
#' a[FALSE] = 0           # no items to replace,
#' print(a, format = 'c') # a is not changed
#' 
#' a[lower.tri(matrix(0, nrow = 3, ncol = 2))] = 0 # set elements below the 
#' print(a, format = 'c')                          # diagonal equal to zero
#'
#' a[3,1] = c(1,-1)       # set (3,1) element
#' print(a, format = 'c') # equal to (1 - z)
#' 
#' a[1:2, 2:1] = c(0,1)   # set the elements in rows 1,2 and coluimns 1,2
#' print(a, format = 'c') # equal to z
#' 
#' a[2, ] = test_polm(dim = c(1,2), degree = 4)
#' print(a, format = 'c')
#' 
#' a[, 1] = test_polm(dim = c(2,1), degree = 4) # this gives a warning
#' print(a, format = 'c')
#' 
#' \dontrun{
#' a[i=1] = 0   # named arguments are not supported
#' }  
"[<-.polm" = function(x,i,j,value) {
  names_args = names(sys.call())
  # print(sys.call())
  # print(names_args)
  if (!all(names_args[-length(names_args)] == '')) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 2
  # print(n_args)
  is_missing = c(missing(i), missing(j))
  # print(is_missing)
  
  x = unclass(x)
  dx = dim(x)
  m = dx[1]
  n = dx[2]
  p = dx[3] - 1 
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  idx = as.vector(idx)
  # print(length(idx))
  
  if (!inherits(value,'polm')) {
    # coerce right hand side 'value' to polm object  
    value = try( polm(value) )
    if (inherits(value, 'try-error')) {
      stop('Could not coerce the right hand side to a polm object!')
    }
  }
  value = unclass(value)
  dv = dim(value)
  
  # no items to replace: return original object
  if (length(idx) == 0) return(polm(x))
  
  if ((dv[1]*dv[2]) == 0) stop('replacement has length 0')
  
  # bring degrees of 'x' and of 'value' in line
  if (dv[3] > dx[3]) {
    x = dbind(d = 3, x, array(0, dim = c(dx[1], dx[2], dv[3] - dx[3])) )
    p = dv[3] - 1
  }
  if (dv[3] < dx[3]) {
    value = dbind(d = 3, value, array(0, dim = c(dv[1], dv[2], dx[3] - dv[3])) )
  }
  
  # coerce 'x' and 'value' to "vectors"
  dim(x) = c(m*n, p+1)
  dim(value) = c(dv[1]*dv[2], p+1)
  # print(value)
  
  # extend 'value' if needed
  if ( (length(idx) %% nrow(value)) != 0 ) {
    warning('number of items to replace is not a multiple of replacement length')
  }
  value = value[rep_len(1:nrow(value), length(idx)),,drop = FALSE]
  # print(value)
  
  # plug in new values
  x[idx,] = value
  # reshape
  dim(x) = c(m, n, p+1)
  # re-coerce to polm
  x = polm(x)
  
  return(x) 
}

#' @rdname replace
#' @export
#' 
#' @examples 
#' (lp = lpolm(1:3, min_deg = -1))
#' lp[1] = 0
#' lp
#' 
#' (lp = lpolm(1:3, min_deg = -1))
#' lp[1] = polm(1:4)
#' lp
#' 
#' (lp = lpolm(1:3, min_deg = -1))
#' lp[1] = lpolm(1, min_deg = 0)
#' lp
'[<-.lpolm' = function(x, i, j, value){
  
  # Check inputs
  names_args = names(sys.call())
  if (!all(names_args[-length(names_args)] == '')) {
    stop('named dimensions are not supported')
  }

  # Info about inputs
  n_args = nargs() - 2
  is_missing = c(missing(i), missing(j))
  
  # Extract attributes and transform to array
  attr_x = attributes(x)
  attributes(x) = NULL
  dim_x = attr_x$dim
  dim(x) = dim_x
  min_deg_x = attr_x$min_deg
  
  # x = unclass(x)
  # dx = dim(x)
  # m = dx[1]
  # n = dx[2]
  # p = dx[3] - 1 
  
  # Extract linear indices
  idx = try(extract_matrix_(dim_x[1], dim_x[2], n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  idx = as.vector(idx)

  if (inherits(value, "polm")){
    value = as.lpolm(value)
  }
  
  if (!inherits(value,'lpolm')) {
    # coerce right hand side 'value' to lpolm object  
    value = try(lpolm(value, min_deg = 0) )
    if (inherits(value, 'try-error')) {
      stop('Could not coerce the right hand side to a lpolm object with min_deg = 0!')
    }
  }
  attr_v = attributes(value)
  attributes(value) = NULL
  dim_v = attr_v$dim
  dim(value) = dim_v
  min_deg_v = attr_v$min_deg
  
  # value = unclass(value)
  # dim_v = dim(value)
   
  # no items to replace: return original object
  if (length(idx) == 0) return(lpolm(x, min_deg = min_deg_x))
  
  if ((dim_v[1]*dim_v[2]) == 0) stop('replacement has length 0')

  p1 = dim_x[3]-1+min_deg_x
  p2 = dim_v[3]-1+min_deg_v 
  
  min_q = min(min_deg_x, min_deg_v)
  max_p = max(p1, p2)
  
  x = dbind(d = 3, 
             array(0, dim = c(dim_x[1], dim_x[2], -min_q+min_deg_x)),
             x,
             array(0, dim = c(dim_x[1], dim_x[2], max_p-p1)))
  value = dbind(d = 3, 
             array(0, dim = c(dim_v[1], dim_v[2], -min_q+min_deg_v)),
             value,
             array(0, dim = c(dim_v[1], dim_v[2], max_p-p2)))
  dim_x = dim(x)
  dim_v = dim(value)
# 
#   
#   # bring degrees of 'x' and of 'value' in line
#   if (dim_v[3] > dim_x[3]) {
#     x = dbind(d = 3, 
#               x, 
#               array(0, dim = c(dim_x[1], dim_x[2], dim_v[3] - dim_x[3])))
#     dim_x[3] = dim_v[3]
#   }
#   if (dim_v[3] < dim_x[3]) {
#     value = dbind(d = 3, 
#                   value, 
#                   array(0, dim = c(dim_v[1], dim_v[2], dim_x[3] - dim_v[3])))
#     dim_v[3] = dim_x[3]
#   }
  
  # coerce 'x' and 'value' to "vectors"
  dim(x) = c(dim_x[1]*dim_x[2], max(dim_x[3], dim_v[3]))
  #cat(dim_v)
  #cat(dim(value))
  dim(value) = c(dim_v[1]*dim_v[2], dim_v[3])

  # extend 'value' if needed
  if ( (length(idx) %% nrow(value)) != 0 ) {
    warning('number of items to replace is not a multiple of replacement length')
  }
  value = value[rep_len(1:nrow(value), length(idx)),,drop = FALSE]
  # print(value)
  
  # plug in new values
  x[idx,] = value
  # reshape
  dim(x) = c(dim_x[1], dim_x[2], dim_x[3])
  # re-coerces to polm
  x = lpolm(x, min_deg = min_q)
  
  return(x) 
  
}



#' Extract Parts of a Rational Matrix
#' 
#' The subsetting operation \code{x[,]} for rational matrices works analogously to the subsetting of 
#' ordinary matrices. However, this operator is only partly implemented 
#' for \code{\link{lmfd}} and \code{\link{lmfd}} objects. See the details below.
#' \cr
#' The \code{$} operator may be used to extract the polynomial factors of a left/right 
#' matrix fraction description (\code{\link{lmfd}}, \code{\link{rmfd}} object). Furthermore  
#' one may retrieve the parameter matrices \eqn{A,B,C,D} of a state space representation 
#' (\code{\link{stsp}} object). For an \code{\link{zvalues}} object, we may 
#' access the complex numbers at which the rational matrix has been evaluated.
#' 
#' \itemize{
#' \item \code{x[]} or \code{x[,]} simply return the original object.
#' \item \code{x[i]} returns a "vector", i.e. an \eqn{(s,1)} dimensional matrix. 
#' \item \code{x[i,]}, \code{x[,j]} or \code{x[i,j]} return a rational matrix 
#'       with rows selected by \code{i} and columns selected by \code{j}. 
#' \item \strong{Note:} for \code{\link{lmfd}} objects (\code{\link{rmfd}} objects) 
#'       only extraction of columns (respectively rows) is implemented. Therefore, 
#'       e.g. \code{x[i,]} throws an error if \code{x} is an \code{lmfd} object.
#' \item Note that "named" arguments are not supported (in order to simplify the coding). 
#' \item In order to have a finer control, one may e.g. use \code{unclass(x)[,,]}.
#' \item \code{x$a}, \code{x$b} returns the left, respectively right factor of an 
#'       LMFD \eqn{x(z) = a^{-1}(z)b(z)}. (\code{x} is an \code{\link{lmfd}} object).
#' \item \code{x$c}, \code{x$d} returns the right, respectively left factor of an 
#'       RMFD  \eqn{x(z) = d(z)c^{-1}(z)}. (\code{x} is an \code{rmfd} object.)
#' \item \code{x$A}, \code{x$B}, \code{x$C}, \code{x$D} return the 
#'       parameter matrices of the statespace realization  
#'       \eqn{x(z) = D + z (I- Az)^{-1}B}. (\code{x} is an \code{stsp} object.)
#' \item If \code{x} is an \code{\link{zvalues}} object, then 
#'       \code{x$z} returns the vector of complex points at which the 
#'       rational matrix \eqn{x} has been evaluated. Furthermore \code{x$f} gives 
#'       the corresponding "frequencies", i.e. \code{x$f = -Arg(x$z)/(2*pi)}.
#' }       
#'       
#'  
#' 
#' 
#'
#' @param x a rational matrix, i.e. a \code{\link{polm}}, \code{\link{lpolm}}, 
#'          \code{\link{lmfd}}, \code{\link{rmfd}}, \code{\link{stsp}}, 
#'          \code{\link{pseries}} or \code{\link{zvalues}} object.
#' @param i,j indices (integer or boolean vector) 
#' @param name character: A,B,C,D for \code{stsp} objects, 
#'             a,b for \code{lmfd} objects, 
#'             c,d for \code{rmfd} objects and 
#'             z,f for \code{zvalues} objects.
#'
#' @return 
#' The subsetting operation \code{x[i,j]} returns a rational matrix of the same class 
#' as the input \code{x}. The mode of the output of the \code{$} operator depends on 
#' which "component" of the rational matrix is extracted. See the details above.
#'  
#' @rdname extract
#' @name extract
#' @export
#' 
#' @examples 
#' # polynomial matrices 
#' a = test_polm(dim = c(3,2), degree = 1)
#' a[]          # returns a
#' a[,]         # returns a
#' a[c(1,3,6)]  # returns a "vector" with the (1,1), (3,1) and (3,2) element of a
#' a[1,]        # returns the first row of a 
#' a[,2]        # returns the second column of a 
#' a[c(TRUE,FALSE,TRUE),c(FALSE, TRUE)] # returns a 2 by 1 matrix 
#' a[c(1,1),c(2,1)] # returns a 2 by 2 matrix 
#' # check with pseries/zvalues
#' all.equal(pseries(a[c(1,1),c(2,1)]), pseries(a)[c(1,1),c(2,1)])
#' all.equal(zvalues(a[c(1,1),c(2,1)]), zvalues(a)[c(1,1),c(2,1)])
#' 
#' \dontrun{
#' a[i=1, j=2] # throws an error, since "named" arguments are not allowed.
#' }
'[.polm' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)
  
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  
  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(x)
  }
  
  if (n_args == 1) { 
    dim(x) = c(m*n,1,p+1)
    x = polm(x[i,1,,drop = FALSE])
    return(x)
  }
  
  if (is_missing[1]) {
    # x[,j]
    x = polm(x[,j,,drop = FALSE])
    return(x)
  }
  
  if (is_missing[2]) {
    # x[i,]
    x = polm(x[i,,,drop=FALSE])
    return(x)
  }
  
  # x[i,j]
  x = polm(x[i,j,,drop=FALSE])
  return(x)
}

#' @rdname extract
#' @export
'[.lpolm' = function(x,i,j){

  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)
  
  attr_x = attributes(x)
  min_deg = attr_x$min_deg
  attributes(x) = NULL
  d = attr_x$dim
  dim(x) = d
  
  # x = unclass(x)
  # d = dim(x)
  # m = d[1]
  # n = d[2]
  # p = d[3] - 1
  
  idx = try(extract_matrix_(d[1], d[2], n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  
  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(x)
  }
  
  if (n_args == 1) { 
    dim(x) = c(d[1]*d[2],1,d[3])
    x = lpolm(x[i,1,,drop = FALSE], min_deg = attr_x$min_deg)
    return(x)
  }
  
  if (is_missing[1]) {
    # x[,j]
    x = lpolm(x[,j,,drop = FALSE], min_deg = attr_x$min_deg)
    return(x)
  }
  
  if (is_missing[2]) {
    # x[i,]
    x = lpolm(x[i,,,drop=FALSE], min_deg = attr_x$min_deg)
    return(x)
  }
  
  # x[i,j]
  x = lpolm(x[i,j,,drop=FALSE], min_deg = attr_x$min_deg)
  return(x)


}

#' @rdname extract
#' @export
#' @examples
#' 
#' # the subsetting operator [,] is only implemented for "lmfd" columns
#' (l = test_lmfd(dim = c(2,2), degrees = c(1,1)))
#' l[,1]
'[.lmfd' = function(x,i,j) {
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  if(!all(is_missing == c(TRUE, FALSE)) && n_args != 2){
    stop("Only columns of lmfd()s can be subset.")
  }
  
  lmfd_obj = x
  x = lmfd_obj$b
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  
  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(lmfd(a = lmfd_obj$a, b = x))
  }
  
  if (is_missing[1]) {
    # x[,j]
    x = polm(x[,j,,drop = FALSE])
    return(lmfd(a = lmfd_obj$a, b = x))
  }
  
  stop("This should not happen.")  
}

#' @rdname extract
#' @export
#' @examples
#' 
#' # the subsetting operator [,] is only implemented for "rmfd" rows
#' (r = test_rmfd(dim = c(2,2), degrees = c(1,1)))
#' r[1,]
'[.rmfd' = function(x,i,j) {
  
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  if(!all(is_missing == c(FALSE, TRUE)) && n_args != 2){
    stop("Only rows of rmfd()s can be subset.")
  }
  
  rmfd_obj = x
  x = rmfd_obj$d
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  
  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(rmfd(c = rmfd_obj$c, d = x))
  }
  
  if (is_missing[2]) {
    # x[i,]
    x = polm(x[i,,,drop = FALSE])
    return(rmfd(c = rmfd_obj$c, d = x))
  }
  
  stop("This should not happen.")  
}

#' @rdname extract
#' @export
'[.stsp' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # print(is_missing)
  
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  
  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = ncol(idx)),
             C = matrix(0, nrow = nrow(idx), ncol = 0), D = matrix(0, nrow = nrow(idx), ncol = ncol(idx)))
    return(x)
  }
  
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)
  
  x = unclass(x)
  A = x[iseq(1,s), iseq(1,s), drop = FALSE]
  B = x[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = x[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = x[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]
  
  if (n_args == 1) {
    transposed = FALSE
    if (m < n) {
      # in order to get a statespace model with a "small" statespace dimension
      # transpose matrix
      transposed = TRUE
      A = t(A)
      D = t(D)
      junk = B
      B = t(C)
      C = t(junk)
      i = matrix(1:(m*n), nrow = m, ncol = n)
      i = t(i)[idx]
      m = nrow(D)
      n = ncol(D)
    }
    A = do.call('bdiag', rep(list(A), n))
    C = do.call('bdiag', rep(list(C), n))
    dim(B) = c(s*n,1)
    dim(D) = c(m*n,1)
    if (transposed) {
      A = t(A)
      D = t(D)
      junk = B
      B = t(C)
      C = t(junk)
    }
    C = C[i,,drop = FALSE]
    D = D[i,,drop = FALSE]
    x = stsp(A = A, B = B, C = C, D = D)
    return(x)
  }
  
  if (is_missing[1]) {
    # x[,j]
    B = B[,j,drop = FALSE]
    D = D[,j,drop = FALSE]
    x = stsp(A = A, B = B, C = C, D = D)
    return(x)
  }
  if (is_missing[2]) {
    # x[i,]
    C = C[i,,drop = FALSE]
    D = D[i,,drop = FALSE]
    x = stsp(A = A, B = B, C = C, D = D)
    return(x)
  }
  # x[i,j]
  B = B[,j,drop = FALSE]
  C = C[i,,drop = FALSE]
  D = D[i,j,drop = FALSE]
  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}

#' @rdname extract
#' @export
'[.pseries' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)
  
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  lag.max = d[3] - 1
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  
  if (length(idx) == 0) {
    # result is an "empty" impulse response
    x = array(0, dim = c(nrow(idx), ncol(idx), lag.max+1))
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }
  
  if (n_args == 1) { 
    dim(x) = c(m*n, 1, lag.max + 1)
    x = x[i,1,,drop = FALSE]
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }
  
  if (is_missing[1]) {
    # x[,j]
    x = x[,j,,drop = FALSE]
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }
  
  if (is_missing[2]) {
    # x[i,]
    x = x[i,,,drop=FALSE]
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }
  
  # x[i,j]
  x = x[i,j,,drop=FALSE]
  x = structure(x, class = c('pseries', 'ratm'))
  return(x)
}

#' @rdname extract
#' @export
'[.zvalues' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)
  
  x = unclass(x)
  z = attr(x,'z')
  d = dim(x)
  m = d[1]
  n = d[2]
  n.z = d[3]
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  
  if (length(idx) == 0) {
    # result is an "empty" frequency response
    x = array(0, dim = c(nrow(idx), ncol(idx), n.z))
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }
  
  if (n_args == 1) { 
    dim(x) = c(m*n, 1, n.z)
    x = x[i,1,,drop = FALSE]
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }
  
  if (is_missing[1]) {
    # x[,j]
    x = x[,j,,drop = FALSE]
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }
  
  if (is_missing[2]) {
    # x[i,]
    x = x[i,,,drop=FALSE]
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }
  
  # x[i,j]
  x = x[i,j,,drop=FALSE]
  x = structure(x, z = z, class = c('zvalues','ratm'))
  return(x)
}


#' @rdname extract
#' @export
'$.lmfd' = function(x, name) {
  i = match(name, c('a','b'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  x = unclass(x)
  if (i == 1) {
    return(polm(array(x[,iseq(1,m*(p+1))], dim = c(m,m,p+1))))
  }
  if (i == 2) {
    return(polm(array(x[,iseq(m*(p+1)+1,m*(p+1)+n*(q+1))], dim = c(m,n,q+1))))
  }
  # this should not happen
  stop('unknown reference')
}


#' @rdname extract
#' @export
'$.rmfd' = function(x, name) {
  i = match(name, c('c','d'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  x = t(unclass(x)) # transposing is necessary because RMFDs are stacked one above the other
  if (i == 1) {
    w = array(x[,iseq(1,n*(p+1))], dim = c(n,n,p+1))
    return(polm(aperm(w, c(2,1,3))))
  }
  if (i == 2) {
    w = array(x[,iseq(n*(p+1)+1,n*(p+1)+m*(q+1))], dim = c(n,m,q+1))
    return(polm(aperm(w, c(2,1,3))))
  }
  # this should not happen
  stop('unknown reference')
}

#' @rdname extract
#' @export
'$.stsp' = function(x, name) {
  i = match(name, c('A','B','C','D'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  s = d[3]
  x = unclass(x)
  if (i == 1) {
    return(x[iseq(1,s),iseq(1,s),drop = FALSE])
  }
  if (i == 2) {
    return(x[iseq(1,s),iseq(s+1,s+n),drop = FALSE])
  }
  if (i == 3) {
    return(x[iseq(s+1,s+m),iseq(1,s),drop = FALSE])
  }
  if (i == 4) {
    return(x[iseq(s+1,s+m),iseq(s+1,s+n),drop = FALSE])
  }
  # this should not happen
  stop('unknown reference')
}

#' @rdname extract
#' @export
'$.zvalues' = function(x, name) {
    i = match(name, c('z','f'))
    if (is.na(i)) stop('reference to "',name, '" is not defined')
    z = attr(x, 'z')
    if (i == 1) {
      return(z)
    }
    if (i == 2) {
      return( -Arg(z)/(2*pi) )
    }
    # this should not happen
    stop('unknown reference')
  }
  