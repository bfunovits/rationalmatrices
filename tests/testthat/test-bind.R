context("test-bind")

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

test_that('rbind', {
  for (i in (1:100)) {
    err = TRUE
    while (err) {
      x = r_matrix()
      dim_x = dim(x)[1:2]
      
      dim_y = c(sample(0:4,1), dim_x[2])
      y = r_matrix(dim_y)
      
      dim_z = c(sample(0:4,1), dim_x[2])
      z = r_matrix(dim_z)
      
      # at least one of x,y,z must be a rational matrix object!
      err = !( inherits(x,'ratm') || inherits(y,'ratm') || inherits(z, 'ratm') )
    }
    
    expect_equal(pseries(rbind(x)), rbind(pseries(x)))
    if (!is.matrix(x)){
      expect_equal(unclass(zvalues(rbind(x))) + complex(real = 0), 
                   unclass(rbind(zvalues(x))) + complex(real = 0))
    }
    
    
    expect_equal(pseries(rbind(x,y)), rbind(pseries(x),pseries(y)))
    if (!any(is.matrix(x), is.matrix(y))){
      expect_equal(unclass(zvalues(rbind(x,y))) + complex(real = 0), 
                   unclass(rbind(zvalues(x),zvalues(y))) + complex(real = 0))
    }
    
    expect_equal(pseries(rbind(x,y,z)), rbind(pseries(x),pseries(y),pseries(z)))
    if (!any(is.matrix(x), is.matrix(y), is.matrix(z))){
      expect_equal(unclass(zvalues(rbind(x,y,z))) + complex(real = 0), 
                 unclass(rbind(zvalues(x),zvalues(y),zvalues(z))) + complex(real = 0))
    }
  }
})


test_that('cbind', {
  for (i in (1:100)) {
    err = TRUE
    while (err) {
      x = r_matrix()
      dim_x = dim(x)[1:2]
      
      dim_y = c(dim_x[1], sample(0:4,1))
      y = r_matrix(dim_y)
      
      dim_z = c(dim_x[1], sample(0:4,1))
      z = r_matrix(dim_z)
      
      # at least one of x,y,z must be a rational matrix object!
      err = !( inherits(x,'ratm') || inherits(y,'ratm') || inherits(z, 'ratm') )
    }
    
    expect_equal(pseries(cbind(x)), cbind(pseries(x)))
    if (!is.matrix(x)){
      expect_equal(unclass(zvalues(cbind(x))) + complex(real = 0), 
                   unclass(cbind(zvalues(x))) + complex(real = 0))
    }
    
    expect_equal(pseries(cbind(x,y)), cbind(pseries(x),pseries(y)))
    if (!any(is.matrix(x), is.matrix(y))){
      expect_equal(unclass(zvalues(cbind(x,y))) + complex(real = 0), 
                   unclass(cbind(zvalues(x),zvalues(y))) + complex(real = 0))
    }
    
    expect_equal(pseries(cbind(x,y,z)), cbind(pseries(x),pseries(y),pseries(z)))
    if (!any(is.matrix(x), is.matrix(y), is.matrix(z))){
      expect_equal(unclass(zvalues(cbind(x,y,z))) + complex(real = 0), 
                   unclass(cbind(zvalues(x),zvalues(y),zvalues(z))) + complex(real = 0))
    }
  }
})

set.seed(NULL)
