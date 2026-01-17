# ===================================================================
# Functional Area #10: Utilities & Helpers (CLAUDE.md)
#
# Purpose: Implement display and printing methods
#   - print(): Format and display rational matrices
#   - show(): S4 method display
#   - Character string conversion
#
# Related Files: 10_utilities_str.R (structure display), 10_utilities_extract.R (subsetting)
# ===================================================================

# Helpers ##############################################################

#' Coerce Scalar Polynomials to Character Strings
#'
#' This utility coerces a scalar polynomial (given by the vector of real coefficients) or a scalar Laurent polynomial to a character string. 
#' If the vector should be interpreted as a Laurent polynomial, the minimal degree should be given in the argument \code{laurent} described below.
#' The following "formats" are implemented.
#' \itemize{
#'   \item \code{syntax = "txt"} returns a simple text representation,
#'   \item \code{syntax = "TeX"} renders the coefficients to string in "TeX" syntax and
#'   \item \code{syntax = "expression"} gives a string which may be rendered to
#'         an \code{R} expression with \code{\link[base]{parse}}. This expression
#'         may be used to evaluate the polynomial and for annotating plots,
#'         see \code{\link[grDevices]{plotmath}} and the examples below.
#' }
#'
#' @param coefs Vector of doubles (complex elements are not allowed) representing a univariate polynomial (first element corresponds to power zero).
#'          When it represents a Laurent polynomial, the first element corresponds to the minimal degree.
#' @param syntax (character string) determines the format of the output string.
#' @param x (character string) polynomial variable.
#' @param laurent Boolean or integer. Default set to FALSE.
#'                If one deals with Laurent polynomials, an integer corresponding to the minimal degree of the Laurent polynomial should be supplied.
#'                
#' @return Character string used for printing univariate (Laurent) polynomials.
#' @export
#'
#' @examples
#' coefs = c(1, 2.3, 0, -1, 0)
#'
#' as_txt_scalarpoly(coefs, syntax = 'txt', x = 'x')
#' as_txt_scalarpoly(coefs, syntax = 'TeX', x = '\\alpha')
#' as_txt_scalarpoly(coefs, syntax = 'expression', x = 'z')
#' as_txt_scalarpoly(coefs = sample((-10):10, 7, replace = TRUE), 
#'                   syntax = 'txt', x = 'x', 
#'                   laurent = -3)
#'
#' \dontrun{
#' # the case syntax = "expression" may be used e.g. as follows
#'
#' # make_polyfun creates a "closure" which evaluates the polynomial at given points
#' # note that this simple version does not work for zero polynomials!
#' make_polyfun = function(coefs) {
#'   expr = parse(text = as_txt_scalarpoly(coefs, 'expression', 'x'))
#'   fun = function(x) {
#'     return(eval(expr))
#'   }
#'   return(fun)
#' }
#'
#' a = make_polyfun(coefs)
#' a(1)   # return the value  of the polynomial at x = 1
#' a(1:5) # return the values of the polynomial at x = 1,2,3,4,5
#'
#' # create a plot
#' x_grid = seq(from = -1, to = 1, length.out = 101)
#' plot(x_grid, a(x_grid), type = 'l', xlab = 'x', ylab = 'a(x)',
#'      main = parse(text = paste('a(x) == ',
#'           as_txt_scalarpoly(coefs, syntax = 'expression', x = 'x'))))
#' }
as_txt_scalarpoly = function(coefs, 
                             syntax = c('txt', 'TeX', 'expression'), x = 'z', 
                             laurent = FALSE) {
  coefs = as.vector(coefs)
  if (!is.numeric(coefs)) {
    stop('"coefs" must be a numeric vector')
  }
  syntax = match.arg(syntax)
  
  if ((length(coefs) == 0) || all(coefs == 0)) {
    return('0')
  }
  
  p = length(coefs)-1
  if (laurent){
    powers = laurent:(laurent+p)
  } else {
    powers = (0:p)
  }
  
  # skip zero coefficients
  non_zero = (coefs != 0)
  coefs = coefs[non_zero]
  powers = powers[non_zero]
  
  # convert powers to character strings
  if (syntax == 'txt') {
    # x^k
    powers_txt = paste(x, '^', powers, sep = '')
  } else {
    # x^{k}
    powers_txt = paste(x, '^{', powers, '}', sep = '')
    # fmt = 'x^{k}' # never used according to RStudio
  }
  powers_txt[powers == 0] = ''
  powers_txt[powers == 1] = x
  powers = powers_txt
  
  signs = ifelse(coefs < 0, '- ', '+ ')
  signs[1] = ifelse(coefs[1] < 0, '-', '')
  
  # convert coefficients to character strings
  coefs = paste(abs(coefs))
  coefs[ (coefs == '1') & (powers != '') ] = ''
  
  if (syntax == 'expression') {
    mults = rep('*', length(coefs))
    mults[ (coefs == '') | (powers == '') ] = ''
  } else {
    mults = rep('', length(coefs))
  }
  
  txt = paste(signs, coefs, mults, powers, sep = '', collapse = ' ')
  return(txt)
  
}

#' Coerce Scalar Polynomial Filters to Character Strings
#'
#' This utility coerces a scalar polynomial filter (given by the vector
#' of coefficients) to a character string. The following "formats" are implemented.
#' \code{syntax = "txt"} returns a simple text representation,
#' \code{syntax = "TeX"} renders the coefficients to string in "TeX" syntax and
#' \code{syntax = "expression"} gives a string which may be rendered to
#' an \code{R} expression with \code{\link[base]{parse}}. This expression
#' may be used to evaluate the filter and for annotating plots,
#' see \code{\link[grDevices]{plotmath}} and the examples below.
#'
#' @param coefs (numeric) vector of coefficients.
#' @param syntax (character string) determines the format of the output string.
#' @param x (character string) names the "input" series.
#' @param t (character string) names the "time-index".
#'
#' @return character string.
#' @export
#'
#' @seealso \code{\link{as_txt_scalarpoly}} and \code{\link{as_tex_matrixfilter}}.
#'
#' @examples
#' coefs = c(1, 2.3, 0, -1, 0)
#'
#' as_txt_scalarfilter(coefs, syntax = 'txt', x = 'x', t = 't')
#' as_txt_scalarfilter(coefs, syntax = 'TeX', x = 'x', t = 's')
#' as_txt_scalarfilter(coefs, syntax = 'expression', x = 'x', t = 'k')
#'
#' \dontrun{
#' # the case syntax = "expression" may be used e.g. as follows
#'
#' # make_filterfun creates a "closure" which computes the filter-output
#' # note that this simple version does not work for zero filters!
#' make_filterfun = function(coefs) {
#'   p = length(coefs) - 1
#'   expr = parse(text = as_txt_scalarfilter(coefs, 'expression', 'x', 't'))
#'   fun = function(x, t) {
#'     # x, t must be vectors
#'     y = rep(NA_real_, length(t))
#'     t0 = t
#'     y = rep(NA_real_, length(t))
#'
#'     i = ((t0 > p) & (t0 <= length(x)))
#'     t = t0[i]
#'     if (any(i)) y[i] = eval(expr)
#'
#'     return(y)
#'   }
#'   return(fun)
#' }
#'
#' coefs = rep(1, 4) / 4  # represents a moving average of length 4.
#' a = make_filterfun(coefs)
#' u = rnorm(100)       # input series
#' a(u, 1)    # return the value of the output series at t = 1
#'            # this value is not defined due to missing initial values
#' a(u, 1:10) # return the values of the output series at t = 1,..,10
#'
#' # create a plot
#' plot(1:length(u), u, type = 'n', xlab = 'time', ylab = '')
#' grid()
#' lines(1:length(u), u, col = 'black', lwd = 1)
#' lines(1:length(u), a(u, 1:length(u)), col = 'red', lwd = 2)
#' legend('topright', bty = 'n',
#'        fill = c('black', 'red'),
#'        legend = c(expression(u[t]),
#'                   parse(text = paste('x[t] == ',
#'                      as_txt_scalarfilter(coefs, 'expression', 'u','t')))) )
#' }
as_txt_scalarfilter = function(coefs, syntax = c('txt', 'TeX', 'expression'),
                               x = 'z', t = 't') {
  coefs = as.vector(coefs)
  if (!is.numeric(coefs)) {
    stop('"coefs" must be a numeric vector')
  }
  syntax = match.arg(syntax)
  
  if ((length(coefs) == 0) || all(coefs == 0)) {
    return('0')
  }
  lags = (0:(length(coefs)-1))
  
  # skip zero coefficients
  non_zero = (coefs != 0)
  coefs = coefs[non_zero]
  lags = lags[non_zero]
  
  # convert lags to character strings
  if (syntax == 'TeX') {
    # x_{t-k}
    lags_txt = paste(x, '_{', t, '-', lags, '}', sep = '')
    lags_txt[lags == 0] = paste(x, '_{', t, '}', sep = '')
  } else {
    # x[t-k]
    lags_txt = paste(x, '[', t, '-', lags, ']', sep = '')
    lags_txt[lags == 0] = paste(x, '[', t, ']', sep = '')
  }
  lags = lags_txt
  
  signs = ifelse(coefs < 0, '- ', '+ ')
  signs[1] = ifelse(coefs[1] < 0, '-', '')
  
  # convert coefficients to character strings
  coefs = paste(abs(coefs))
  coefs[ (coefs == '1') & (lags != '') ] = ''
  
  if (syntax == 'expression') {
    mults = rep('*', length(coefs))
    mults[ (coefs == '') | (lags == '') ] = ''
  } else {
    mults = rep('', length(coefs))
  }
  txt = paste(signs, coefs, mults, lags, sep = '', collapse = ' ')
  return(txt)
  
}

#' TeX Matrix
#'
#' @param x matrix, where \code{paste(x[i,j])} returns a valid "TeX" string.
#'
#' @return character string.
#' @export
#'
#' @examples
#' as_tex_matrix(diag(1:2, nrow = 2, ncol = 3))
as_tex_matrix = function(x) {
  if ( !is.matrix(x) ) stop('"x" must be a matrix')
  
  m = nrow(x)
  n = ncol(x)
  
  if (length(x) == 0) return('\\begin{pmatrix}\n\\end{pmatrix}')
  
  tex = '\\begin{pmatrix}\n'
  for (i in (1:m)) {
    tex = paste(tex, paste(x[i,], collapse = ' & '), '\\\\\n', sep = '  ')
  }
  tex = paste(tex, '\\end{pmatrix}', sep = '')
  return(tex)
}


#' TeX Matrix Polynomials
#'
#' @param coefs 3-dimensional array with the coefficients of the matrix polynomial.
#' @param x (character string) polynomial variable or process variable.
#' @param as_matrix_of_polynomials boolean.
#'
#' @return character string.
#' @export
#'
#' @examples
#' coefs = array(round(rnorm(2*3*1), 1), dim = c(2,3,2))
#'
#' as_tex_matrixpoly(coefs)
#' as_tex_matrixpoly(coefs, x = 'x', as_matrix_of_polynomials = FALSE)
as_tex_matrixpoly = function(coefs, x = 'z', as_matrix_of_polynomials = TRUE) {
  # only some basic checks
  if ( (!is.array(coefs)) || (length(dim(coefs)) != 3) || (!is.numeric(coefs)) ) {
    stop('"coefs" must be 3-dimensional numeric array')
  }
  
  d = dim(coefs)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  
  if ((m*n) == 0) {
    return('\\begin{pmatrix}\n\\end{pmatrix}')
  }
  
  if ((m*n) == 1) {
    return(as_txt_scalarpoly(coefs, syntax = 'TeX', x = x))
  }
  
  if ((p < 0) || all(coefs == 0)) {
    return(as_tex_matrix(matrix(0, nrow = m, ncol = n)))
  }
  
  if (as_matrix_of_polynomials) {
    tex = apply(coefs, MARGIN = c(1,2), FUN = as_txt_scalarpoly,
                syntax = 'TeX', x = x)
    tex = as_tex_matrix(tex)
    return(tex)
  }
  
  # print as polynomial with matrix coefficients
  
  powers = (0:p)
  # coerce powers to character strings of the form x^{k}
  powers_txt = paste(x, '^{', powers, '}', sep = '')
  powers_txt[powers == 0] = ''
  powers_txt[powers == 1] = x
  powers = powers_txt
  
  tex = ''
  for (k in (0:p)) {
    a = matrix(coefs[,,k+1], nrow = m, ncol = n)
    if ( !all(a == matrix(0, nrow = m, ncol = n)) ) {
      # non-zero coefficient matrix
      if (tex != '' ) tex = paste(tex, '+\n')
      if ( (m == n) && all(a == diag(m)) ) {
        # coefficient matrix is identity matrix
        tex = paste(tex, ' I_{', m, '} ', powers[k+1], sep = '')
      } else {
        tex = paste(tex, as_tex_matrix(a), powers[k+1])
      }
    }
  }
  
  return(tex)
}


#' TeX Matrix Polynomial Filters
#'
#' @param coefs 3-dimensional array with the coefficients of the filter.
#' @param x (character string) polynomial variable or process variable.
#' @param t (character string) time/index variable.
#'
#' @return character string.
#' @export
#'
#' @examples
#' coefs = array(round(rnorm(2*3*1), 1), dim = c(2,3,2))
#'
#' as_tex_matrixfilter(coefs, x = '\\epsilon', t = 's')
as_tex_matrixfilter = function(coefs, x = 'z', t = 't') {
  # only some basic checks
  if ( (!is.array(coefs)) || (length(dim(coefs)) != 3) || (!is.numeric(coefs)) ) {
    stop('"coefs" must be 3-dimensional numeric array')
  }
  
  d = dim(coefs)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  
  if ((m*n) == 0) {
    return('\\begin{pmatrix}\n\\end{pmatrix}')
  }
  
  if ((m*n) == 1) {
    return(as_txt_scalarfilter(coefs, syntax = 'TeX', x = x, t = t))
  }
  
  if ((p < 0) || all(coefs == 0)) {
    tex = as_tex_matrix(matrix(0, nrow = m, ncol = n))
    return( paste(tex, ' ', x, '_{', t, '}', sep = '') )
  }
  
  lags = (0:p)
  # coerce lags to character strings of the form x_{t-k}
  lags_txt = paste(x, '_{', t, '-', lags, '}', sep = '')
  lags_txt[lags == 0] = paste(x, '_{', t, '}', sep = '')
  lags = lags_txt
  
  tex = ''
  for (k in (0:p)) {
    a = matrix(coefs[,,k+1], nrow = m, ncol = n)
    if ( !all(a == matrix(0, nrow = m, ncol = n)) ) {
      # non-zero coefficient matrix
      if (tex != '' ) tex = paste(tex, '+\n')
      if ( (m==n) && all(a == diag(m)) ) {
        # coefficient matrix is identity matrix
        tex = paste(tex, ' I_{', m, '} ', lags[k+1], sep = '')
      } else {
        tex = paste(tex, as_tex_matrix(a), lags[k+1])
      }
    }
  }
  
  return(tex)
}


#' Helper for Printing Matrix (Laurent) Polynomials
#' 
#' Called by the \code{rationalmatrices} methods for the \link{print} S3 generic.
#'
#' @param a Array with 3 dimensions, representing a polynomial or Laurent polynomial
#' @param digits Default NULL. Otherwise an integer specifying the number of digits of the coefficients that should be printed, 
#'               see \code{\link{round}}.
#' @param format Default "i|jz", i.e. the matrix polynomial is printed that the "z part is the slowest moving index.
#'               The vertical bar designates separation across array dimensions (one bar = two array dimensions, two bars = array of dimension 3).
#'               In addition to multiple ordering options regarding row-index \code{i}, column-index \code{j}, and polynomial power \code{z}, 
#'               there is also the option \code{character} (partial matching is enabled such that \code{format = "c"} gives the desired result)
#'               which prints each univariate (Laurent) polynomial into a matrix with appropriate numbers of rows and columns.
#'               Note that \code{\link{stsp}} objects have no format option. 
#'               The option \code{'character'} is only implemented for polynomials, Laurent polynomials, LMFDs and RMFDs with real coefficients.
#' @inheritParams as_txt_scalarpoly  
#'
#' @return Printed (Laurent) polynomial matrix
#' @export
#' @keywords internal
print_3D = function(a, digits = NULL, 
                    format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), 
                    laurent = FALSE) {
  dim = dim(a)
  m = dim[1]
  n = dim[2]
  p = dim[3]
  # empty array -> do nothing
  if (min(dim) == 0) return(invisible(NULL))
  
  # a must have full 'dimnames'
  names = dimnames(a)
  inames = names[[1]]
  jnames = names[[2]]
  znames = names[[3]]
  
  # round 
  if (!is.null(digits)) a = round(a, digits)
  
  format = match.arg(format)
  
  if (format == 'character') {

    # convert vector of coefficients to character representation of a polynomial
    a = apply(a, MARGIN = c(1,2), FUN = as_txt_scalarpoly, syntax = "txt", x = "z", 
              laurent = laurent)
    
    # add column names (jnames)
    a = rbind(jnames, a)
    
    # add row names (inames)
    a = cbind( c('',inames), a)
    
    # right justify columns 
    w = nchar(a)
    w = apply(w, MARGIN = 2, FUN = max)
    for (j in (1:(n+1))) {
      fmt = paste('%', w[j], 's', sep='')
      pad = function(s) { sprintf(fmt, s) }
      a[,j] = apply(a[,j,drop = FALSE], MARGIN = 1, FUN = pad)
    }

    # convert matrix a to a string
    a = apply(a, MARGIN = 1, FUN = paste, collapse = '  ')
    a = paste(a, collapse = '\n')
    cat(a,'\n')
  }
  
  
  if (format == 'i|jz') {
      # create a vector of the form 
      # j[1],...,j[n],j[1],...,j[n],...
      jnames = rep(jnames, p)
      # create a vector of the form 
      # z[1],'',...,'',z[2],'',...,'',
      if (n > 1) {
        znames = as.vector(rbind(znames, 
                                 matrix('', nrow = n-1, ncol = p)))
      }
      
      dim(a) = c(m,n*p)
      rownames(a) = inames
      colnames(a) = paste(znames, jnames, sep = ' ')
      print(a)  
  }
  
  if (format == 'i|zj') {
    # create a vector of the form 
    # z[1],...,z[p],z[1],...,z[p],...
    znames = rep(znames, n)
    # create a vector of the form 
    # j[1],'',...,'',j[2],'',...,'',
    if (p > 1) {
      jnames = as.vector(rbind(jnames, 
                               matrix('', nrow = p-1, ncol = n)))
    }
    
    a = aperm(a, c(1,3,2))
    dim(a) = c(m,p*n)
    rownames(a) = inames
    colnames(a) = paste(jnames, znames, sep = ' ')
    print(a)  
  }
  
  if (format == 'iz|j')  {
    # create a vector of the form 
    # i[1],...,i[m],i[1],...,i[m],...
    inames = rep(inames, p)
    # create a vector of the form 
    # z[1],'  ',...,'  ',z[2],'  ',...,'  ',
    if (m > 1) {
      znames = as.vector(rbind( znames, 
                                matrix(' ', nrow = m-1, ncol = p)))
    }
    # right justify
    fmt = paste('%', max(nchar(znames)), 's', sep='')
    pad = function(s) { sprintf(fmt, s) }
    znames = as.vector(apply(matrix(znames, ncol = 1), MARGIN = 1, FUN = pad))
    
    a = aperm(a, c(1,3,2))
    dim(a) = c(m*p, n)
    rownames(a) = paste(znames, inames, sep = ' ')
    colnames(a) = jnames
    print(a)
  }
  
  if (format == 'zi|j')  {
    # create a vector of the form 
    # z[1],...,z[p],z[1],...,z[p],...
    znames = rep(znames, m)
    # create a vector of the form 
    # i[1],'  ',...,'  ',i[2],'  ',...,'  ',
    if (p > 1) {
      inames = as.vector(rbind( inames, 
                                matrix(' ',nrow = p-1, ncol = m)))
    }
    # right justify
    fmt = paste('%', max(nchar(inames)), 's', sep='')
    pad = function(s) { sprintf(fmt, s) }
    inames = as.vector(apply(matrix(inames, ncol = 1), MARGIN = 1, FUN = pad))

    a = aperm(a, c(3,1,2))
    dim(a) = c(p*m, n)
    rownames(a) = paste(inames, znames, sep = ' ')
    colnames(a) = jnames
    print(a)
  }
  
  if (format == 'i|j|z') {
    # the last case 'i|j|z' just uses the R default print of 3D array
    print(a)
  }
    
  return(invisible(NULL))
}

# print.___() for rationalmatrices objects ####

#' Print Methods
#' 
#' Printing rational matrix objects.
#'
#' @param x rational matrix object, i.e. a \code{\link{polm}}, \code{\link{lpolm}}, \code{\link{lmfd}}, \code{\link{rmfd}},
#'          \code{\link{stsp}}, \code{\link{pseries}} or \code{\link{zvalues}} object.
#' @param digits (integer) if non \code{NULL} then correspondingly rounded numbers are printed, 
#'        see \code{\link{round}}.
#' @param format (character string) selects specific output formats. Note that 
#'        \code{\link{stsp}} objects have no format option. The option \code{'character'} 
#'        is only implemented for polynomials, Laurent polynomials, LMFDs and RMFDs with real coefficients, se
#' @param ... Further parameters are ignored.
#'
#' @return \code{invisible(x)}
#'
#' @rdname print
#' @name print methods
#' 
#' @examples
#' # for polynomials six different print formats are implemented ###################
#' a = test_polm(dim = c(2,3), degree = 2, random = TRUE)
#' 
#' for (fmt in c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character")) {
#'    cat('\nformat =', fmt, '\n')
#'    print(a, digits = 2, format = fmt)
#' }
#' 
#' # "empty" (2 x 0) polynomial matrix (degree = 2)
#' a = test_polm(dim = c(2,0), degree = 0)
#' print(a)
#' 
#' # random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
#' a = polm(array(complex(real = stats::rnorm(2*1*3), 
#'                        imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
#' print(a, digits = 2)
#' \dontrun{
#' # the format option 'character' is only implemented for polynomials matrices 
#' # with real coefficients!
#' print(a, digits = 2, format = 'character')
#' }
#'
#' # print a rational matrix in statespace form
#' a = test_stsp(dim = c(3,3), s = 2)
#' print(a, digits = 2)
#' 
#' # print a rational matrix in 'lmfd' form 
#' a = test_lmfd(dim = c(2,3), degrees = c(2,1))
#' print(a, digits = 2, format = 'character')
#' 
#' # print impulse response 
#' print(pseries(a), format = 'i|zj', digits = 2)
#' 
#' # print frequency response 
#' print(zvalues(a), format = 'iz|j', digits = 2)
#' 
NULL


#' @rdname print
#' @export
print.lpolm = function(x, digits = NULL, 
        format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  a = unclass(x)
  min_deg = attr(x, which = "min_deg")
  attr(a, 'min_deg') = NULL # remove 'min_deg' attribute
  m = dim(a)[1]
  n = dim(a)[2]
  p = dim(a)[3]-1+min_deg
  
  cat('( ', m, ' x ', n,' ) Laurent polynomial matrix with degree <= ', p, 
      ', and minimal degree >= ', min_deg, '\n', sep = '')
  if ((m*n*(dim(a)[3])) == 0) {
    return(invisible(x))
  }
  
  if ((format == 'character') && (is.complex(a))) {
    stop(paste('the format option "character" is only implemented',
               'for Laurent polynomials with real coefficients'))
  }
  
  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('z^', min_deg:p, sep = ''))
  print_3D(a, digits, format, laurent = min_deg)
  
  invisible(x)
}

#' @rdname print
#' @export
print.polm = function(x, digits = NULL, 
        format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  p = dim(a)[3]-1
  
  cat('(',m,'x',n,') matrix polynomial with degree <=', p,'\n')
  if ((m*n*(p+1)) == 0) {
    return(invisible(x))
  }
  
  if ((format == 'character') && (is.complex(a))) {
    stop(paste('the format option "character" is only implemented',
               'for polynomials with real coefficients'))
  }
  
  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('z^',0:p, sep = ''))
  print_3D(a, digits, format)
  
  invisible(x)
}


#' @rdname print
#' @export
print.lmfd = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  cat('( ', m, ' x ', n,' ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = ', 
      p, ', q = ', q, ')\n', sep = '')
  
  if ((format == 'character') && (is.complex(unclass(x)))) {
    stop('the format option "character" is only implemented for LMFDs with real coefficients')
  }

  if ((m*m*(p+1)) > 0) {
    cat('left factor a(z):\n')
    
    a = unclass(x$a)
    
    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(a, digits, format)
  }
  
  if ((m*n*(q+1)) > 0) {
    cat('right factor b(z):\n')
    
    a = unclass(x$b)
    
    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }
  
  invisible(x)
}

#' @rdname print
#' @export
print.rmfd = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  cat('( ', m, ' x ', n,' ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = ', 
      p, ', deg(d(z)) = q = ', q, '\n', sep = '')
  
  if ((format == 'character') && (is.complex(unclass(x)))) {
    stop('the format option "character" is only implemented for RMFDs with real coefficients')
  }
  
  
  if ((m*n*(q+1)) > 0) {
    cat('left factor d(z):\n')
    
    d = unclass(x$d)
    
    # use the above defined internal function print_3D
    dimnames(d) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(d, digits, format)
  }
  
  if ((n*n*(p+1)) > 0) {
    cat('right factor c(z):\n')
    
    c = unclass(x$c)
    
    # use the above defined internal function print_3D
    dimnames(c) = list(paste('[', 1:n, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(c, digits, format)
  }
  
  invisible(x)
}


#' @rdname print
#' @export
print.stsp = function(x, digits = NULL, ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]
  
  cat('statespace realization [', m, ',', n, '] with s = ', s, ' states\n', sep = '')
  
  a = unclass(x)
  attr(a, 'order') = NULL
  if (length(a) == 0) {
    return(invisible(x))
  }
  
  # rounding digits
  if (!is.null(digits)) {
    a = round(a, digits)
  }
  
  snames = character(s)
  if (s > 0) snames = paste('s[',1:s,']',sep = '')
  xnames = character(m)
  if (m > 0) xnames = paste('x[',1:m,']',sep = '')
  unames = character(n)
  if (n > 0) unames = paste('u[',1:n,']',sep = '')
  
  rownames(a) = c(snames, xnames)
  colnames(a) = c(snames, unames)
  print(a)
  
  invisible(x)
}


#' @rdname print
#' @export
print.pseries = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  lag.max = dim(a)[3]-1
  
  cat('(',m,'x',n,') impulse response with maximum lag =', lag.max,'\n')
  if ((m*n*(lag.max+1)) == 0) {
    return(invisible(x))
  }
  
  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('lag=',0:lag.max, sep = ''))
  print_3D(a, digits, format)

  invisible(x)
}

#' @rdname print
#' @export
print.zvalues = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  z = attr(x, 'z')
  n.z = length(z)
  
  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  attr(a, 'z') = NULL # remove 'z' attribute
  
  cat('(',m,'x',n,') frequency response\n')
  if ((m*n*n.z) == 0) {
    return(invisible(x))
  }
  
  # use the above defined internal function print_3D
  if ((format == 'i|jz') || (format == 'i|zj')) {
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z[',1:n.z, ']', sep = ''))
  } else {
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z=', round(z,3), sep = ''))
  }
  print_3D(a, digits, format)

  invisible(x)
}
