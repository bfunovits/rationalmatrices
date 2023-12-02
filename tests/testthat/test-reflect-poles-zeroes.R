test_that("reflect_poles_zeroes works", {
  # only check reflect unstable poles/zeroes for "stsp" objects
  set.seed(12345)
  for (i in (1:100)) {
    m = sample(1:5,1)
    s = sample(1:5,1) 
    
    K = test_stsp(dim = c(m,m), s = s)
    
    # poles of K(z)
    pK = poles(K)
    # zeroes of K(z)
    zK = zeroes(K)
    # frequency responses 
    Kf = zvalues(K)
    
    # reflect all unstable poles (inside the unit circle) ###########
    reflect = (abs(pK) < 1) & (Im(pK) >= 0)
    KU = reflect_poles(K, poles = pK[reflect])
    
    # frequency responses 
    KUf = zvalues(KU)
    
    # check poles 
    testthat::expect_gt(min(abs(poles(KU))), 1)
    
    # check that U(z) is all-pass
    testthat::expect_equal(Kf %r% Ht(Kf), KUf %r% Ht(KUf))
    
    
    # reflect all unstable zeroes (inside the unit circle) ###########
    reflect = (abs(zK) < 1) & (Im(zK) >= 0)
    KU = reflect_zeroes(K, zeroes = zK[reflect])
    
    # frequency responses 
    KUf = zvalues(KU)
    
    # check zeroes
    testthat::expect_gt(min(abs(zeroes(KU))), 1)
    
    # check that U(z) is all-pass
    testthat::expect_equal(Kf %r% Ht(Kf), KUf %r% Ht(KUf))
  }
  
  set.seed(NULL)
})
