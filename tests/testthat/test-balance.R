# test-balance.R
# 

# create a statespace model with a given number of 
# controllable and observable states (s[1])
# controllable but not observable states (s[2])
# not controllable and not observable states (s[3])
# not controllable but observable states (s[4])
test_stsp2 = function(dim, s) {
  m = dim[1]
  n = dim[2]
  ss = sum(s)
  if ( ss == 0 ) {
    return(test_stsp(dim = c(m,n), s = ss))
  }
  
  err = TRUE
  sd = 2
  i.trial = 0
  while ((err) & (i.trial < 100)) {
    A = matrix(rnorm(ss*ss, sd = sd), nrow = ss, ncol = ss)
    B = matrix(rnorm(ss*n, sd = sd), nrow = ss, ncol = n)
    C = matrix(rnorm(m*ss, sd = sd), nrow = m, ncol = ss)
    D = matrix(rnorm(m*n), nrow = m, ncol = n)
    
    # non-controllabe
    i = iseq(s[1]+s[2]+1, ss)
    j = (1:ss)[-i]
    B[i,] = 0
    A[i,j] = 0
    
    # non-observable
    i = iseq(s[1]+1,s[1]+s[2]+s[3])
    j = (1:ss)[-i]
    C[,i] = 0
    A[j,i] = 0
    
    err = (max(abs(eigen(A, only.values = TRUE)$values)) >= 1)
    if (m == n) {
      Ainv = try(A - B %*% solve(D, C))
      err = (err || inherits(Ainv,'try-error') || (max(abs(eigen(Ainv, only.values = TRUE)$values)) >= 1) )
    }
    sd = sd / 1.1
    i.trial = i.trial +1 
  }
  if (err) {
    print(i.trial)
    stop('could not generate model')
  }
  x = stsp(A,B,C,D)
  # print(eigen(x$A)$values)
  T = matrix(rnorm(ss*ss), nrow = ss, ncol = ss)
  x = state_trafo(x, T)
  # print(eigen(x$A)$values)
  return(x)  
}

set.seed(12345)

n.sim = 100
n.err = c(0,0)
for (i in (1:n.sim)) {

  s = c( sample(1:4,1), sample(0:1,1), sample(0:1,1), sample(0:1,1) )
  x = test_stsp2(dim=c(sample(1:3,1), sample(1:3,1)), s)
  # print(x, digits = 2)
  # 
  gr = grammians(x)
  
  # balance, no truncation, true statespace dimension
  bal = balance(x, gr, s = s[1], truncate = FALSE)
  gr2 = grammians(bal$obj)
  
  testthat::expect_equivalent(bal$P, gr2$P)
  testthat::expect_equivalent(bal$Q, gr2$Q)
  testthat::expect_equivalent(pseries(x), pseries(bal$obj))
  
  # balance, truncation, true statespace dimension
  bal = balance(x, gr, s = s[1], truncate = TRUE)
  gr2 = grammians(bal$obj)
  
  testthat::expect_equivalent(bal$P, gr2$P)
  testthat::expect_equivalent(bal$Q, gr2$Q)
  testthat::expect_equivalent(pseries(x), pseries(bal$obj))
  
  # balance, truncation, estimate statespace dimension
  bal = balance(x, gr, truncate = TRUE)
  gr2 = grammians(bal$obj)
  testthat::expect_equivalent(pseries(x), pseries(bal$obj))
  
  # minimal statespace dimension has not been correctly estaimated
  n.err = n.err + c((dim(bal$obj)[3] > s[1]), (dim(bal$obj)[3] < s[1]))
}
cat('\nNumber of errors (', paste(n.err, collapse = ','), ') in ', n.sim, ' simulations\n', sep = '')

set.seed(NULL)
