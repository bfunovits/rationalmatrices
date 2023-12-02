context("lpolm(), its methods, and transformations")
library(rationalmatrices)

lpolm_methods = ls(getNamespace("rationalmatrices"), all.names = TRUE)
lpolm_methods[grep("lpolm", lpolm_methods)]

# Some functions cannot be found with the command above
# Group methods, rbind(), cbind()
# get_bwd(), get_fwd(), polm2fwd()

# OK: test_lpolm() Tests for Laurent polynomials ####

# OK: Extract methods ####
(lp = test_lpolm(dim = c(3,2),
                 degree_max = 2,
                 degree_min = -2))
lp %>% print(format = "c", digits = 2)

# Case 1: Everything missing
lp[]          # returns lp
lp[,]         # returns lp

# Case 2: One argument
lp[c(1,2,4)] %>% print(digits = 2)  # returns a "vector" with the (1,1), (2,1) and (2,2) element of lp
lp[c(1,2,4)] %>% print(format = "c", digits = 2)  # returns a "vector" with the (1,1), (2,1) and (2,2) element of lp

# Case 3: 2 arguments, one missing
# Matrices (vectors with additional dim and class attribute) are treated as vectors
lp %>% print(digits = 2)
lp %>% print(format = "c", digits = 2)
lp[lower.tri(matrix(0, nrow = 3, ncol = 2))]
lp[lower.tri(matrix(0, nrow = 3, ncol = 2))] %>% print(format = "c", digits = 2)        # returns the first row of lp

lp %>% print(digits = 2)
lp %>% print(format = "c", digits = 2)
lp[1:2,]        # returns the first row of lp
lp[1:2,] %>% print(format = "c", digits = 2)        # returns the first row of lp

lp %>% print(digits = 2)
lp %>% print(format = "c", digits = 2)
lp[,2]        # returns the second column of lp
lp[,2]  %>% print(format = "c", digits = 2)       # returns the second column of lp

# Case 3: 2 arguments, none missing
lp %>% print(digits = 2)
lp %>% print(format = "c", digits = 2)
lp[c(TRUE,FALSE),c(FALSE, TRUE)] # returns a 1 by 1 Laurent polynomial matrix
lp[c(TRUE,FALSE),c(FALSE, TRUE)] %>% print(format = "c", digits = 2)

lp %>% print(digits = 2)
lp %>% print(format = "c", digits = 2)
lp[c(1,1),c(2,1)] # returns a 2 by 2 matrix
lp[c(1,1),c(2,1)] %>% print(format = "c", digits = 2)

# OK: "Replace" methods ####

(a = test_lpolm(dim = c(3,2), degree_max = 1, degree_min = -1))
a[FALSE] = 0   # no items to replace, a is not changed
a

(a = test_lpolm(dim = c(3,2), degree_max = 1, degree_min = -1))
a[lower.tri(matrix(0, nrow = 3, ncol = 2))] = 0 # set elements below the diagonal equal to zero
a

(a = test_lpolm(dim = c(3,2), degree_max = 1, degree_min = -1))
a[3,1] = c(1,-1) # set (3,1) element
a

(a = test_lpolm(dim = c(3,2), degree_max = 1, degree_min = -1))
a[1:2, 2:1] = c(0,1)
a

(a = test_lpolm(dim = c(3,2), degree_max = 1, degree_min = -1))
a[2, ] = test_lpolm(dim = c(1,2), degree_max = 4, degree_min = -2)
a

(a = test_lpolm(dim = c(3,2), degree_max = 1, degree_min = -1))
a[, 1] = test_lpolm(dim = c(3,1), degree_max = 4, degree_min = -2) # this gives a warning
a

# OK: "as.lpolm.polm" ####
#
# Comments:
# Returns a lpolm object with min_deg = 0 (this is not possible with lpolm which returns a polm objct when min_deg = 0)

(p = test_polm(dim = c(2,2),
               degree = 3))
lp = as.lpolm(p)
expect_true(is.lpolm(lp))

# OK: "as.polm.lpolm" throws error if min_deg < 0  ####
(lp = test_lpolm(dim = c(1,1),
                degree_max = 2,
                degree_min = -2))
expect_error(as.polm(lp))

# OK: "derivative.lpolm" not implemented ####
(lp = test_lpolm(dim = c(4,3),
                 degree_max = 2,
                 degree_min = -1))
expect_error(derivative(lp))

# OK: "bind" ####
(lp1 = test_lpolm(dim = c(2,2),
                  degree_max = 1,
                  degree_min = -1))
(lp2 = test_lpolm(dim = c(2,2),
                  degree_max = 1,
                  degree_min = -1))

expect_equal(c(cbind(lp1, lp2)), 
             c(dbind(d=2, unclass(lp1), unclass(lp2))))

expect_equal(c(rbind(lp1, lp2)),
             c(dbind(d=1, unclass(lp1), unclass(lp2))))

(lp1 = test_lpolm(dim = c(1,1),
                  degree_max = 1,
                  degree_min = -2))
(lp2 = test_lpolm(dim = c(1,1),
                  degree_max = 2,
                  degree_min = -1))
cbind(lp1, lp2)
rbind(lp1, lp2)

(lp1 = test_lpolm(dim = c(1,1),
                  degree_max = 2,
                  degree_min = -2))
(lp2 = test_lpolm(dim = c(1,1),
                  degree_max = 1,
                  degree_min = -1))
cbind(lp1, lp2)
rbind(lp1, lp2)

(lp1 = test_lpolm(dim = c(1,1),
                  degree_max = 1,
                  degree_min = -4))
(lp2 = test_lpolm(dim = c(1,1),
                  degree_max = 4,
                  degree_min = -1))
cbind(lp1, lp2)
rbind(lp1, lp2)


# OK: "dim.lpolm" ####
(lp = test_lpolm(dim = c(4,3),
                 degree_max = 2,
                 degree_min = -1))
dim(lp)


# OK: "get_bwd" and "get_fwd" ####

(lp = test_lpolm(degree_max = 2, degree_min = -2))
get_fwd(lp)
get_bwd(lp)

# OK: "Ht.lpolm" ####

(lp = test_lpolm(degree_max = 2, degree_min = -2))
Ht(lp)

# OK: "is.lpolm" ####

(lp = test_lpolm(dim = c(2,2),
                 degree_max = 1,
                 degree_min = -1))
(p = test_polm(dim = c(2,1), degree = 1))

is.lpolm(p)
is.lpolm(lp)

# OK: "lpolm" ####
#
# Comments:
# min_deg can be positive too
# lpolm returns a polm object when min_deg is non-negative!
(lp = lpolm(1:5, min_deg = -7))
(lp = lpolm(1:5, min_deg = 0))
(lp = lpolm(1:5, min_deg = 2))

# OK: "polm2fwd() ####
(p = test_polm(dim = c(1,1), degree = 3))
polm2fwd(p)

# OK: "print.lpolm" ####
#
# Comments:
# Attribute min_deg is printed as well. Delete it?

# (2 x 3) polynomial matrix a(z) = a0 (degree is zero)
(lp = lpolm(diag(1, nrow = 2, ncol = 3), min_deg = 1))
lp %>% str()
lp %>% unclass() %>% str()
print(lp, digits = NULL,
      format = c("i|jz"))
print(lp, digits = NULL,
      format = c("i|zj"))
print(lp, digits = NULL,
      format = c("iz|j"))
print(lp, digits = NULL,
      format = c("i|j|z"))
print(lp, digits = NULL,
      format = c("character"))

# Minimal degree = -1

(lp = lpolm(array(c(matrix(rnorm(8), 2, 4), diag(1, nrow = 2, ncol = 2)), dim = c(2,2,3)), min_deg = -2))
lp %>% class()
lp %>% str()
lp %>% unclass()
lp %>% unclass() %>% str()
print(lp, digits = NULL,
      format = c("i|jz"))
print(lp, digits = NULL,
      format = c("i|zj"))
print(lp, digits = NULL,
      format = c("iz|j"))
print(lp, digits = NULL,
      format = c("i|j|z"))
print(lp, digits = NULL,
      format = c("character"))

# OK: "pseries.lpolm": not possible, not implemented ####
(lp = test_lpolm(dim = c(2,1),
                 degree_max = 2,
                 degree_min = -2))
expect_error(pseries(lp))

# OK: "str.lpolm" ####
# Should have same first line as print!

(lp = test_lpolm(dim = c(2,1),
                 degree_max = 2,
                 degree_min = -2))
str(lp)
print(lp)


# OK: "t.lpolm" ####

# Transpose
t(lp)

# OK: "test_lpolm" ####
#
# Comments:
# Not really an essential function
# Overwriting of column_start_matrix, value_at_0, column_end_matrix and consistency of degrees should be tested
(lp = test_lpolm(dim = c(1,1), degree_max = 1, degree_min = -2))
(lp = test_lpolm(dim = c(3,3), degree_max = c(0,1,2), degree_min = -2))
(lp = test_lpolm(dim = c(3,3), degree_max = 1, degree_min = c(0,-1,-2)))

(lp = test_lpolm(dim = c(1,1), degree_max = 0, degree_min = 0))
(lp = test_lpolm(dim = c(1,1), degree_max = 0, degree_min = -1))
(lp = test_lpolm(dim = c(1,1), degree_max = 0, degree_min = -2))
(lp = test_lpolm(dim = c(1,1), degree_max = 1, degree_min = 0))
(lp = test_lpolm(dim = c(1,1), degree_max = 1, degree_min = -1))
(lp = test_lpolm(dim = c(1,1), degree_max = 1, degree_min = -2))
(lp = test_lpolm(dim = c(1,1), degree_max = 2, degree_min = 0))
(lp = test_lpolm(dim = c(1,1), degree_max = 2, degree_min = -1))
(lp = test_lpolm(dim = c(1,1), degree_max = 2, degree_min = -2))

# max degree cannot be smaller than min degree:
expect_error(test_lpolm(dim = c(1,1), degree_max = -1, degree_min = 0))

(lp = test_lpolm(dim = c(1,1), degree_max = -1, degree_min = -1))
(lp = test_lpolm(dim = c(1,1), degree_max = -1, degree_min = -2))
(lp = test_lpolm(dim = c(2,2), degree_max = 0, degree_min = 0))
(lp = test_lpolm(dim = c(2,2), degree_max = 0, degree_min = -1))
(lp = test_lpolm(dim = c(2,2), degree_max = 0, degree_min = -2))
(lp = test_lpolm(dim = c(2,2), degree_max = 1, degree_min = 0))
(lp = test_lpolm(dim = c(2,2), degree_max = 1, degree_min = -1))
(lp = test_lpolm(dim = c(2,2), degree_max = 1, degree_min = -2))
(lp = test_lpolm(dim = c(2,2), degree_max = 2, degree_min = 0))
(lp = test_lpolm(dim = c(2,2), degree_max = 2, degree_min = -1))
(lp = test_lpolm(dim = c(2,2), degree_max = 2, degree_min = -2))

# max degree cannot be smaller than min degree:
expect_error(lp = test_lpolm(dim = c(2,2), degree_max = -1, degree_min = 0))

(lp = test_lpolm(dim = c(2,2), degree_max = -1, degree_min = -1))
(lp = test_lpolm(dim = c(2,2), degree_max = -1, degree_min = -2))
(lp = test_lpolm(dim = c(3,3), degree_max = c(0,1,2), degree_min = 0))
(lp = test_lpolm(dim = c(3,3), degree_max = c(0,1,2), degree_min = -1))
(lp = test_lpolm(dim = c(3,3), degree_max = c(0,1,2), degree_min = -2))
(lp = test_lpolm(dim = c(3,3), degree_max = 1, degree_min = c(0,-1,-2)))

# min_degree can be zero
(lp = test_lpolm(dim = c(3,2), degree_max = 3, degree_min = 0))

# min_degree cannot be positive
expect_error(test_lpolm(dim = c(3,2), degree_max = 3, degree_min = 1))

# col_start_matrix, value_at_0, col_end_matrix

# degrees as matrices




# OK: "zeroes.lpolm" ####
(lp = test_lpolm(dim = c(2,2),
                 degree_max = 2,
                 degree_min = -2))
zeroes(lp)

# OK: "zvalue.lpolm" ####
(lp = test_lpolm(dim = c(2,2),
                 degree_max = 2,
                 degree_min = -2))
zvalue(lp, z = 1.2)

# OK: "zvalues.lpolm" ####

zvalues(lp, n.f = 5)








