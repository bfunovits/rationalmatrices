test_that("munkres works", {
  for (i in (1:10)) {
    m = sample(3:5, 1)
    n = sample(3:5, 1) 

    C = matrix(sample(1:10, m*n, replace = TRUE), nrow = m, ncol = n)
    # we simply test that the function does not hang and does not throw an error
    expect_named(munkres(C), c("a", "c", "C"))
  }
})
