context("test-power-operation")

# set.seed(04272019)

r_matrix = function(classes = c('polm','lmfd','stsp')) {
  cl = sample(classes,1)
  m = sample(0:4,1)
  n = m
  dim = c(m,n)
  
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
    return(test_stsp(dim = c(m,n), s = s))
  }
}

test_that('x^a', {
  for (i in (1:100)) {
    x = r_matrix()
    m = unname(dim(x))[1]
    
    expect_equal(pseries(x^0), pseries(polm(diag(m))))
    expect_equal(x^1, x)
    expect_equal(pseries(x^2), pseries(x)^2)
    expect_equal(pseries(x^3), pseries(x)^3)
    expect_equal(zvalues(x^2), zvalues(x)^2)
    expect_equal(zvalues(x^3), zvalues(x)^3)
    if (m > 0) {
      expect_equal(pseries(x^(-1)), pseries(x)^(-1))
      expect_equal(pseries(x^(-2)), pseries(x)^(-2))
      expect_equal(zvalues(x^(-1)), zvalues(x)^(-1))
      expect_equal(zvalues(x^(-2)), zvalues(x)^(-2))
    }
  }
})


set.seed(NULL)
