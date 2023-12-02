test_that("lyapunov works", {
  
  set.seed(04272019)
  
  # symmetric case
  for (i in (1:100)) {
    m = sample(1:10, 1)
    n = sample(1:5, 1)
    A = matrix(stats::rnorm(m*m), nrow = m, ncol = m)
    mstable = sample(c(TRUE,FALSE), 1)
    if (mstable) {
      A = A / (1.1*max(abs(eigen(A, only.values = TRUE)$values)))
    }
    Q = crossprod(matrix(stats::rnorm(m*m), nrow = m, ncol = m))
    
    lambda_A = eigen(A, symmetric = FALSE, only.values = TRUE)$values
    is_stable = (max(abs(lambda_A)) < 1)
    P = lyapunov(A, Q)
    testthat::expect_equal(P, A %*% P %*% t(A) + Q)
    if (!is_stable) {
      testthat::expect_warning(lyapunov(A, Q, non_stable = 'warn'))
    }
    
    dA = matrix(rnorm((m^2)*n), nrow = m^2, ncol = n)
    dQ = matrix(rnorm((m^2)*n), nrow = m^2, ncol = n)
    for (i in (1:n)) {
      # make sure that each column of dQ corresponds to a symmetric matrix!
      junk = matrix(dQ[,i], nrow = m, ncol = m)
      junk = junk + t(junk)
      dQ[,i] = junk
    }
    out = lyapunov_Jacobian(A, Q, dA, dQ)
    all.equal(out$P, P)
    
    eps = 1e-8
    theta = rnorm(n)
    dP = matrix(out$J %*% theta, nrow = m, ncol = m)
    
    # compute the derivative via "finite differences"
    P1 = lyapunov(A + matrix(dA %*% theta, nrow = m, ncol = m)*eps, 
                  Q + matrix(dQ %*% theta, nrow = m, ncol = m)*eps)
    all.equal(dP, (P1 - P)/eps, scale = mean(abs(out$J)), tol = 1e-6)
  }
  
  set.seed(NULL)
  
})
