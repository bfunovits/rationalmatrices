context("test-matrix-multiplication")

# set.seed(04272019)

r_matrix = function(dim = NULL, classes = c('matrix','polm','lmfd','stsp')) {
  cl = sample(classes,1)
  if (is.null(dim)) {
    m = sample(0:4,1)
    n = sample(0:4,1)
    dim = c(m,n)
  }
  m = dim[1]
  n = dim[2]
  if (cl == 'matrix') {
    return(matrix(rnorm(m*n), nrow = m, ncol = n))
  }
  if (cl == 'polm') {
    p = sample(0:4,1)
    return(test_polm(dim = dim, degree = p, random = TRUE))
  }
  if (cl == 'lmfd') {
    p = sample(0:4,1)
    q = sample(0:4,1)
    if (m > 0) {
      return(test_lmfd(dim = c(m,n), degrees = c(p,q), bpoles = 1))
    } else {
      return(test_polm(dim = c(m,n), degree = q, random = TRUE))
    }
  }
  if (cl == 'stsp') {
    s = sample(0:4,1)
    return(test_stsp(dim = c(m,n), s= s, bpoles = 1))
  }
}

test_that('x %r% y', {
  for (i in (1:100)) {
    x0 = r_matrix(classes = c('polm','lmfd','stsp'))
    x1 = x0
    if (!inherits(x1,'ratm')) x1 = polm(x1)
    
    dim1 = dim(x1)[1:2]
    dim2 = c(dim1[2], sample(0:4,1))
    
    
    y0 = r_matrix(dim2)
    y1 = y0
    if (!inherits(y1,'ratm')) y1 = polm(y1)
    
    # if ((prod(dim1) > 0) && (prod(dim2) > 0)) {
    #   x2 = x1
    #   if (class(x1)[1] != 'stsp') x2 = as.stsp(x1)
    #   cat(dim(x1), dim(x2), 
    #       max(abs(unclass(zvalues(x1))-unclass(zvalues(x2)))), '\n')
    #   cat(min(c(Inf,abs(poles(x1)))),'\n')
    #   
    #   y2 = y1
    #   if (class(y1)[1] != 'stsp') y2 = as.stsp(y1)
    #   cat(dim(x1), dim(x2), 
    #       max(abs(unclass(zvalues(x1))-unclass(zvalues(x2)))), '\n')
    #   cat(min(c(Inf,abs(poles(y1)))),'\n')
    # }
    
    expect_equal(pseries(x0 %r% y0),  pseries(x1) %r%  pseries(y1))
    expect_equal(x0 %r% pseries(y1),  pseries(x1) %r%  pseries(y1))
    # make sure that zvalues is complex
    expect_equal(unclass(zvalues(x0 %r% y0)) + complex(real = 0), 
                 unclass(zvalues(x1) %r% zvalues(y1)) +  complex(real = 0))
    expect_equal(unclass(x0 %r% zvalues(y1)) + complex(real = 0), 
                 unclass(zvalues(x1) %r% zvalues(y1)) + complex(real = 0))
  }
})


set.seed(NULL)
