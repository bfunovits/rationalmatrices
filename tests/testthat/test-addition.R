context("test-addition")

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

test_that('x + y', {
  for (i in (1:100)) {
    x0 = r_matrix(classes = c('polm','lmfd','stsp'))
    x1 = x0
    if (!inherits(x1,'ratm')) x1 = polm(x1)
    
    dim = dim(x1)[1:2]
    
    y0 = r_matrix(dim)
    y1 = y0
    if (!inherits(y1,'ratm')) y1 = polm(y1)
    
    z0 = r_matrix(dim = c(1,1))
    z1 = z0
    if (!inherits(z1,'ratm')) z1 = polm(z1)
    
    expect_equal(pseries(x0 + y0),  pseries(x1) +  pseries(y1))
    expect_equal(x0 + pseries(y1),  pseries(x1) +  pseries(y1))
    expect_equal(zvalues(x0 - y0), zvalues(x1) - zvalues(y1))
    expect_equal(x0 - zvalues(y1), zvalues(x1) - zvalues(y1))
    expect_equal(pseries(x0 + z0),  pseries(x1) +  pseries(z1))
    expect_equal(zvalues(z0 + x0), zvalues(z1) + zvalues(x1))
  }
})


set.seed(NULL)
