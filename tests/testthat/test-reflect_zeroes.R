test_that("reflect zeroes works", {
  # 
  set.seed(12345)
  
  for (isim in (1:20)) {
    # complex conjugated pair of zeroes with a "real kernel"
    # ######################################################
    m = 3
    alpha = complex(real = rnorm(1), imaginary = rnorm(1))
    a = array(0, dim = c(m, m, 3))
    a[,,1] = diag(m)
    a[1,1,1] = abs(alpha)^2
    a[1,1,2] = -2*Re(alpha)
    a[1,1,3] = 1
    a = polm(a)
    a = test_polm(dim = c(m,m), degree = 1, random = TRUE) %r% a %r% 
      matrix(rnorm(m*m), nrow = m, ncol = m)
    
    a1 = reflect_zeroes(a, alpha, tol = 1e-7)
    
    # Check zeroes 
    # Note that a(z) has only 5 finite zeroes. Therefore we only  
    # keep the 5 smallest ones.
    a_zeroes = zeroes(a, print_message = FALSE)
    a_zeroes = a_zeroes[order(abs(a_zeroes))]
    a_zeroes = a_zeroes[1:5]
    a1_zeroes = zeroes(a1, print_message = FALSE)
    a1_zeroes = a1_zeroes[order(abs(a1_zeroes))]
    a1_zeroes = a1_zeroes[1:5]
    
    r_zeroes = a_zeroes
    r_zeroes[ which.min(abs(r_zeroes - alpha)) ] = 1/alpha
    r_zeroes[ which.min(abs(r_zeroes - Conj(alpha))) ] = 1/Conj(alpha)
    
    j = match_vectors(r_zeroes, a1_zeroes)
    testthat::expect_equal(r_zeroes, a1_zeroes[j], tol = 1e-4, scale = 1) 
    # print(cbind(a_zeroes, r_zeroes, a1_zeroes[j], r_zeroes - a1_zeroes[j]))
    
    # Check that the transformation matrix U (x1 = x %r% U) is all-pass
    testthat::expect_equal(zvalues(a) %r% Ht(zvalues(a)), 
                           zvalues(a1) %r% Ht(zvalues(a1)))
    
    # complex conjugated pair of zeroes with 2-dimensional kernel
    # ###########################################################
    m = 4
    alpha = complex(real = rnorm(1), imaginary = rnorm(1))
    w = matrix(complex(real = rnorm(m*2), imaginary = rnorm(m*2)), nrow = m, ncol = 2)
    
    a = Re( cbind(w, Conj(w)) %*% 
              diag(c(alpha, alpha, Conj(alpha), Conj(alpha))) %*% 
              solve(cbind(w, Conj(w))) )
    a = polm(array(cbind(-a, diag(m)), dim = c(m, m, 2)))
    
    a = test_polm(dim = c(m,m), degree = 1, random = TRUE) %r% a %r% 
      test_polm(dim = c(m,m), degree = 1, random = TRUE)
    
    a1 = reflect_zeroes(a, alpha, tol = 1e-7)
    
    # Check zeroes 
    a_zeroes = zeroes(a, print_message = FALSE)
    r_zeroes = a_zeroes
    r_zeroes[ which.min(abs(r_zeroes - alpha)) ] = 1/alpha
    r_zeroes[ which.min(abs(r_zeroes - Conj(alpha))) ] = 1/Conj(alpha)
    
    a1_zeroes = zeroes(a1, print_message = FALSE)
    
    j = match_vectors(r_zeroes, a1_zeroes)
    testthat::expect_equal(r_zeroes, a1_zeroes[j]) 
    # print(cbind(a_zeroes, r_zeroes, a1_zeroes[j]))
    
    # Check that the transformation matrix U (x1 = x %r% U) is all-pass
    testthat::expect_equal(zvalues(a) %r% Ht(zvalues(a)), 
                           zvalues(a1) %r% Ht(zvalues(a1)))
    
    
    # reflect alpha twice!
    # ###########################################################
    
    a1 = reflect_zeroes(a, c(alpha, alpha), tol = 1e-7)
    
    # Check zeroes 
    a_zeroes = zeroes(a, print_message = FALSE)
    r_zeroes = a_zeroes
    r_zeroes[ which.min(abs(r_zeroes - alpha)) ] = 1/alpha
    r_zeroes[ which.min(abs(r_zeroes - alpha)) ] = 1/alpha
    r_zeroes[ which.min(abs(r_zeroes - Conj(alpha))) ] = 1/Conj(alpha)
    r_zeroes[ which.min(abs(r_zeroes - Conj(alpha))) ] = 1/Conj(alpha)
    
    a1_zeroes = zeroes(a1, print_message = FALSE)
    
    j = match_vectors(r_zeroes, a1_zeroes)
    testthat::expect_equal(r_zeroes, a1_zeroes[j]) 
    # print(cbind(a_zeroes, r_zeroes, a1_zeroes[j]))
    
    # Check that the transformation matrix U (x1 = x %r% U) is all-pass
    testthat::expect_equal(zvalues(a) %r% Ht(zvalues(a)), 
                           zvalues(a1) %r% Ht(zvalues(a1)))
    
    # try stsp method 
    # ##########################################################
    
    a = as.stsp(a)
    a1 = reflect_zeroes(a, c(alpha, alpha))
    
    # Check zeroes 
    a_zeroes = zeroes(a, print_message = FALSE)
    r_zeroes = a_zeroes
    r_zeroes[ which.min(abs(r_zeroes - alpha)) ] = 1/alpha
    r_zeroes[ which.min(abs(r_zeroes - alpha)) ] = 1/alpha
    r_zeroes[ which.min(abs(r_zeroes - Conj(alpha))) ] = 1/Conj(alpha)
    r_zeroes[ which.min(abs(r_zeroes - Conj(alpha))) ] = 1/Conj(alpha)
    
    a1_zeroes = zeroes(a1, print_message = FALSE)
    
    j = match_vectors(r_zeroes, a1_zeroes)
    testthat::expect_equal(r_zeroes, a1_zeroes[j]) 
    # print(cbind(a_zeroes, r_zeroes, a1_zeroes[j]))
    
    # Check that the transformation matrix U (x1 = x %r% U) is all-pass
    testthat::expect_equal(zvalues(a) %r% Ht(zvalues(a)), 
                           zvalues(a1) %r% Ht(zvalues(a1)))
    
  }
  
  set.seed(NULL)
})
