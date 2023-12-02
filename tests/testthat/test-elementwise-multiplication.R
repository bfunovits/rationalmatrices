context("test-elementwise-multiplication")

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
    return(test_stsp(dim = c(m,n), s = s, bpoles = 1))
  }
}


test_that('x * y', {
  for (i in (1:100)) {
    err = TRUE
    while (err) {
      x = r_matrix()
      dim = dim(x)[1:2]
      y = r_matrix(dim)
      # at least one of x,y must be "rational"
      err = !( inherits(x, 'ratm') || inherits(y, 'ratm'))
    }
    
    expect_equal(pseries(x * y),  pseries(x) *  pseries(y))
    expect_equal(x * pseries(y),  pseries(x) *  pseries(y))
    if(!any(is.matrix(x), is.matrix(y))){
      expect_equal(zvalues(x * y), zvalues(x) * zvalues(y))
      expect_equal(x * zvalues(y), zvalues(x) * zvalues(y))
    }
  }
})



test_that('x * y (scalar)', {
  for (i in (1:100)) {
    err = TRUE
    while (err) {
      x = r_matrix()
      y = r_matrix(dim = c(1,1))
      # at least one of x,y must be "rational"
      err = !( inherits(x, 'ratm') || inherits(y, 'ratm'))
    }
    
    expect_equal(pseries(x * y),  pseries(x) *  pseries(y))
    expect_equal(pseries(y * x),  pseries(y) *  pseries(x))
    if(!any(is.matrix(x), is.matrix(y))){
      expect_equal(zvalues(x * y), zvalues(x) * zvalues(y))
      expect_equal(zvalues(y * x), zvalues(y) * zvalues(x))
    }
  }
})


set.seed(NULL)
