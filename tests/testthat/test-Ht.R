test_that("Ht works", {
  set.seed(12345)
  for (isim in (1:20)) {
    # check statespace realizations with complex coefficients
    s = sample(0:5, 1)
    m = sample(0:5, 1)
    n = sample(0:5, 1)
    x = stsp(A = matrix(complex(real = rnorm(s*s), imaginary = rnorm(s*s)), nrow = s, ncol = s), 
             B = matrix(complex(real = rnorm(s*n), imaginary = rnorm(s*n)), nrow = s, ncol = n), 
             C = matrix(complex(real = rnorm(m*s), imaginary = rnorm(m*s)), nrow = m, ncol = s), 
             D = matrix(complex(real = rnorm(m*n), imaginary = rnorm(m*n)), nrow = m, ncol = n))
    expect_equal(Ht(zvalues(x)), zvalues(Ht(x)))
    # check statespace realizations with real coefficients
    s = sample(0:5, 1)
    m = sample(0:5, 1)
    n = sample(0:5, 1)
    x = stsp(A = matrix(rnorm(s*s), nrow = s, ncol = s), 
             B = matrix(rnorm(s*n), nrow = s, ncol = n), 
             C = matrix(rnorm(m*s), nrow = m, ncol = s), 
             D = matrix(rnorm(m*n), nrow = m, ncol = n))
    expect_equal(Ht(zvalues(x)), zvalues(Ht(x)))
  }
  set.seed(NULL)
})
