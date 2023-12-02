# dim.____ methods ##############################################################

#' Dimensions of Objects
#' 
#' Retrieve the dimension and degrees of rational matrix objects. 
#' E.g. for a polynomial matrix \eqn{x(z)} (i.e. for a \code{\link{polm}} object),
#' \code{dim(x)} returns the vector  \code{c(m,n,p)}, 
#' where \code{m,n} are the respective number of rows and columns of the matrix and \code{p} is the (polynomial) degree of the matrix. 
#' For a Laurent polynomial, there is an additional element in the vector \code{c(m,n,p,min_deg)} returned by \code{dim(x)} pertaining to the (possibly negative) minimal degree. 
#' For a rational matrix in RMFD form \eqn{x(z)=d(z)c^{-1}(z)} (i.e. for an \code{\link{rmfd}} object),
#' \code{dim(x)} returns the vector \code{c(m,n,p,q)},
#' where \code{m,n} are the respective number of rows and columns of the matrix and 
#' \code{p,q} are the (polynomial) degrees of the polynomial matrices \eqn{c(z)} and \eqn{d(z)} respectively. 
#'
#' @param x Object of type \code{\link{polm}}, \code{\link{lmfd}}, 
#'          \code{\link{rmfd}}, \code{\link{stsp}}, \code{\link{pseries}}, 
#'          or \code{\link{zvalues}}.
#'
#' @return Returns a named vector of integers (\code{m,n} always refer to the number of 
#'         rows and columns of the rational matrix)
#'         
#' \item{\code{c(m,n,p,min_deg}}{for Laurent polynomial matrices 
#'   (i.e. for \code{\link{lpolm}} objects).}
#' \item{\code{c(m,n,p)}}{for polynomial matrices 
#'   (i.e. for \code{\link{polm}} objects).}
#' \item{\code{c(m,n,p,q)}}{for left matrix fraction descriptions 
#'   (i.e. for \code{\link{lmfd}} objects).}
#' \item{\code{c(m,n,p,q)}}{for right matrix fraction descriptions 
#'   (i.e. for \code{\link{rmfd}} objects).}
#' \item{\code{c(m,n,s)}}{for state space representations 
#'   (i.e. for \code{\link{stsp}} objects).}
#' \item{\code{c(m,n,lag.max)}}{for power series expansions  
#'   (i.e. for \code{\link{pseries}} objects).}
#' \item{\code{c(m,n,n.f)}}{for frequency response functions 
#'   (i.e. for \code{\link{zvalues}} objects).}
#' 
#' @rdname dim
#' @name dim methods
#' @export
dim.polm = function(x) {
  d = dim(unclass(x))
  d[3] = d[3] - 1
  names(d) = c('m','n','p')
  return(d)
}

#' @rdname dim
#' @export
dim.lpolm = function(x) {
  d = dim(unclass(x))
  min_deg = attr(x, which = "min_deg")
  d[3] = d[3] - 1 + min_deg
  d = c(d, min_deg)
  names(d) = c('m','n','p', 'min_deg')
  return(d)
}

#' @rdname dim
#' @export
dim.lmfd = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','p','q')
  return(d)
}


#' @rdname dim
#' @export
dim.rmfd = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','p','q')
  return(d)
}



#' @rdname dim
#' @export
dim.stsp = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','s')
  return(d)
}

#' @rdname dim
#' @export
dim.pseries = function(x) {
  d = dim(unclass(x))
  d[3] = d[3] - 1
  names(d) = c('m','n','lag.max')
  return(d)
}

#' @rdname dim
#' @export
dim.zvalues = function(x) {
  d = dim(unclass(x))
  names(d) = c('m','n','n.f')
  return(d)
}
