# is.____ methods ##############################################################

#' Check Objects
#' 
#' Check if the argument is a valid object. 
#'
#' @param x object
#'
#' @return \code{TRUE} if the argument is a valid object.
#' @rdname is
#' @name is.class methods
#' @export
is.polm = function(x) {
  if (!inherits(x,'polm')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }
  
  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }
  
  return(TRUE)
}

#' @rdname is
#' @export
is.lpolm = function(x) {
  if (!inherits(x,'lpolm')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }
  
  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }
  
  min_deg = attr(x,'min_deg')
  if ( (is.null(min_deg)) || (length(min_deg) != 1) ) {
    return(FALSE)
  }

  return(TRUE)
}


#' @rdname is
#' @export
is.lmfd = function(x) {
  if (!inherits(x,'lmfd')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }
  
  if (!is.matrix(x)) {
    return(FALSE)
  }
  
  order = attr(x,'order')
  if ( (is.null(order)) || (!is.integer(order)) || (length(order) != 4) || 
       (nrow(x) == 0) || (order[1] != nrow(x)) || (order[3] < 0) || 
       (order[2] < 0) || (order[4] < -1) || 
       ((order[1]*(order[3]+1) + order[2]*(order[4]+1)) != ncol(x)) ) {
    return(FALSE)
  }
  
  return(TRUE)
}


#' @rdname is
#' @export
is.rmfd = function(x) {
  if (!inherits(x,'rmfd')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }
  
  if (!is.matrix(x)) {
    return(FALSE)
  }
  
  order = attr(x,'order')
  if ( (is.null(order)) || (!is.integer(order)) || (length(order) != 4) || 
       (ncol(x) == 0) || (order[2] != ncol(x)) || (order[3] < 0) || 
       (order[1] < 0) || (order[4] < -1) || 
       ((order[1]*(order[4]+1) + order[2]*(order[3]+1)) != nrow(x)) ) {
    return(FALSE)
  }
  
  return(TRUE)
}


#' @rdname is
#' @export
is.stsp = function(x) {
  if (!inherits(x,'stsp')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }
  
  if (!is.matrix(x)) {
    return(FALSE)
  }
  
  order = attr(x, 'order')
  if ( (is.null(order)) || (!is.integer(order)) || (length(order) != 3) || 
       (min(order) < 0) || ((order[1]+order[3]) != nrow(x)) || 
       ((order[2]+order[3]) != ncol(x)) ) {
    return(FALSE)
  }
  
  return(TRUE)
}

#' @rdname is
#' @export
is.pseries = function(x) {
  if (!inherits(x,'pseries')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if ( !( is.numeric(x) || is.complex(x) ) ) {
    return(FALSE)
  }
  
  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }
  
  return(TRUE)
}


#' @rdname is
#' @export
is.zvalues = function(x) {
  if (!inherits(x,'zvalues')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)
  
  x = unclass(x)
  if ( !( is.numeric(x) || is.complex(x) ) ) {
    return(FALSE)
  }
  
  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }
  
  z = attr(x, 'z')
  if ( (is.null(z)) || ( !( is.numeric(z) || is.complex(z) ) ) || 
       (!is.vector(z)) || (length(z) != dim(x)[3]) ) {
    return(FALSE)
  } 

  return(TRUE)
}


#' Check Properties of Rational Matrices
#'
#' @param x a rational matrix object, i.e. a \code{\link{polm}}, \code{\link{lmfd}}, 
#'          \code{\link{stsp}}, \code{\link{pseries}}, or \code{\link{zvalues}} object.
#' @param ... not used. 
#'
#' @return Boolean. Note \code{is.stable}, \code{is.miniphase} return \code{NA} 
#'         for objects of type \code{pseries} and \code{zvalues}.
#' @rdname check
#' @name check objects
#' @export
#'
#' @examples
#' a = test_polm(dim = c(2,2), deg = 1, random = TRUE)
#' b = test_polm(dim = c(2,3), deg = 1, random = TRUE)
#' c = lmfd(a,b)
#' 
#' is.stable(a)
#' is.stable(c)
#' is.stable(as.stsp(c))
#' is.stable(pseries(c))
#'
#' is.miniphase(b[,1:2])
#' is.miniphase(as.stsp(c)[,1:2])
#' 
#' \dontrun{
#' is.miniphase(b)
#' is.miniphase(c)
#' }
is.stable = function(x, ...) {
    UseMethod("is.stable", x)
}

#' @rdname check
#' @export
is.stable.ratm = function(x, ...) {
  if (class(x)[1] %in% c('polm','lmfd','rmfd','stsp')) {
    z = poles(x)
    return(all(abs(z) > 1))
  }
  return(NA)
}

#' @rdname check
#' @export
is.miniphase = function(x, ...) {
  UseMethod("is.miniphase", x)
}

#' @rdname check
#' @export
is.miniphase.ratm = function(x, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1]==0)) stop('the miniphase property is only defined for square non-empty matrices')
  
  if (class(x)[1] %in% c('polm','lmfd','rmfd','stsp')) {
    z = zeroes(x)
    return(all(abs(z) > 1))
  }
  return(NA)
}

