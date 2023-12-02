context("test-coprime")

# test is.coprime

# test some error conditions
expect_error(is.coprime(stsp(D=1)))
expect_error(is.coprime(a = test_polm(dim = c(2,1)), test_polm(dim = c(1,2))))
expect_error(is.coprime(a = test_polm(dim = c(0,1), degree = -1)))

# zeroe polynomials
expect_false( is.coprime(a = test_polm(dim = c(2,1), degree = -1)) )
expect_false( is.coprime(a = polm(matrix(0, nrow = 2, ncol = 3))) )

# constant polynomials
expect_false( is.coprime(a = polm(matrix(rnorm(3*2), nrow = 3, ncol = 2))) )
expect_true( is.coprime(a = polm(matrix(rnorm(3*2), nrow = 2, ncol = 3))) )
expect_false( is.coprime(a = polm(matrix(rnorm(2), nrow = 2, ncol = 1) %*% 
                                  matrix(rnorm(3), nrow = 1, ncol = 3))) )

# c0 = 0
expect_length(is.coprime(c(0,1,-1), only.answer = FALSE)$zeroes, 2)

a = array(rnorm(2*3*3), dim =c(2,3,3))
a[,,1] = 0
a = polm(a)
out = is.coprime(a, only.answer = FALSE, debug = FALSE)             
expect_false( out$answer )
