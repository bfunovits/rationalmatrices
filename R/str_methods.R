# str.____ methods ##############################################################

#' Display the Structure of Objects
#'
#' @param object an object
#' @param ... not used
#'
#' @return invisible(NULL)
#' 
#' @rdname str
#' @name str methods
#' @export
str.lpolm = function(object, ...) {
  d = dim(unclass(object))
  min_deg = attr(object, "min_deg")
  cat('( ',d[1],' x ',d[2],' ) Laurent polynomial matrix with degree <= ', 
      d[3]-1+min_deg, ', and minimum degree >= ', min_deg, '\n', sep = '')
  return(invisible(NULL))  
}

#' @rdname str
#' @export
str.polm = function(object, ...) {
  d = dim(unclass(object))
  cat('( ',d[1],' x ',d[2],' ) matrix polynomial with degree <= ', d[3]-1,'\n', sep = '')
  return(invisible(NULL))  
}

#' @rdname str
#' @export
str.lmfd = function(object, ...) {
  d = attr(object, 'order')
  cat('( ',d[1],' x ',d[2],' ) left matrix fraction description with degrees (p = ', 
      d[3], ', q = ', d[4],')\n', sep = '')
  return(invisible(NULL))  
}

#' @rdname str
#' @export
str.rmfd = function(object, ...) {
  d = attr(object, 'order')
  cat('( ',d[1],' x ',d[2],' ) right matrix fraction description with degrees (deg(c(z)) = p = ', 
      d[3], ', deg(d(z)) = q = ', d[4],')\n', sep = '')
  return(invisible(NULL))  
}


#' @rdname str
#' @export
str.stsp = function(object, ...) {
  d = attr(object, 'order')
  cat('( ',d[1],' x ',d[2],' ) statespace realization with s = ', d[3], ' states\n', sep = '')
  return(invisible(NULL))
}

#' @rdname str
#' @export
str.pseries = function(object, ...) {
  d = dim(unclass(object))
  cat('( ',d[1],' x ',d[2],' ) power series parameters with maximum lag = ', d[3]-1, '\n', sep = '')
  return(invisible(NULL))
}

#' @rdname str
#' @export
str.zvalues = function(object, ...) {
  d = dim(unclass(object))
  cat('( ',d[1],' x ',d[2],' ) functional values with ', d[3], ' frequencies/points\n', sep = '')
  return(invisible(NULL))
}
