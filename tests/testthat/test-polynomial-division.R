context("test-polynomial-division")

# set.seed(04272019)

r_matrix = function(dim = NULL, classes = c('matrix','polm')) {
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
}

test_that('x // y', {
  for ( i in (1:100) ) {
    x = r_matrix()
    dim = unname(dim(x))[1:2]
    m = dim[1]
    n = dim[2]
    y = r_matrix(dim)
    z = r_matrix(c(1,1))
    
    x1 = x
    if (!inherits(x1,'polm')) x1 = polm(x)
    y1 = y
    if (!inherits(y1,'polm')) y1 = polm(y)
    z1 = z
    if (!inherits(z1,'polm')) z1 = polm(z)
    
    if ( (class(x)[1] == 'polm') || (class(y)[1] == 'polm') ) {
      q = x %/% y
      r = x %% y
      expect_equal(c(m,n), unname(dim(q)[1:2]))
      expect_equal(c(m,n), unname(dim(r)[1:2]))
      if (m*n > 0) {
        expect_true(all(degree(r) < degree(y1)))
        expect_equal(prune(x1) + 0, prune(y1 * q + r))
      }
    }
    
    if ( (class(x)[1] == 'polm') || (class(z)[1] == 'polm') ) {
      q = x %/% z
      r = x %% z
      expect_equal(c(m,n), unname(dim(q)[1:2]))
      expect_equal(c(m,n), unname(dim(r)[1:2]))
      if (m*n > 0) {
        expect_true(all(degree(r) < degree(z1)[1,1]))
        expect_equal(prune(x1) + 0, prune(z1 * q + r))
      }
      
      q = z %/% x
      r = z %% x
      expect_equal(c(m,n), unname(dim(q)[1:2]))
      expect_equal(c(m,n), unname(dim(r)[1:2]))
      if (m*n > 0) {
        expect_true(all(degree(r) < degree(x1)))
        expect_equal(prune(z1 * matrix(1, nrow = m, ncol = n)) + 0, prune(x1 * q + r))
      }
    }
  }
})

set.seed(NULL)