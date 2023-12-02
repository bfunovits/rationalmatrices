context("test-unary-operations")

# set.seed(04272019)

r_matrix = function(dim = NULL, classes = c('polm','lmfd','stsp')) {
  cl = sample(classes,1)
  if (is.null(dim)) {
    m = sample(0:4,1)
    n = sample(0:4,1)
    dim = c(m,n)
  }
  m = dim[1]
  n = dim[2]
  if (cl == 'polm') {
    p = sample(0:4,1)
    return(test_polm(dim = dim, degree = p, random = TRUE))
  }
  if (cl == 'lmfd') {
    p = sample(0:4,1)
    q = sample(0:4,1)
    if (m > 0) {
      return(test_lmfd(dim = c(m,n), degrees = c(p,q)))
    } else {
      return(test_polm(dim = c(m,n), degree = q, random = TRUE))
    }
  }
  if (cl == 'stsp') {
    s = sample(0:4,1)
    return(test_stsp(dim = c(m,n), s= s))
  }
}

test_that('+-x', {
  for (i in (1:100)) {
    x = r_matrix()
    
    expect_equal(x, +x)
    expect_equal(pseries(-x), -pseries(x))
    expect_equal(zvalues(-x), -zvalues(x))
  }
})

set.seed(NULL)

