# Create an All-Pass Rational Matrix

Create a square, all-pass rational matrix in statespace form with given
\\A,B\\.

## Usage

``` r
make_allpass(A, B)
```

## Arguments

- A:

  square \\(s,s)\\ matrix

- B:

  \\(s,m)\\ matrix

## Value

A
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
object which represents the all-pass rational matrix.

## Notes

The function is intended as an internal helper function and thus does
not check the inputs. See also
[`reflect_poles`](https://bfunovits.github.io/rationalmatrices/reference/reflect_poles.md)
and
[`reflect_zeroes`](https://bfunovits.github.io/rationalmatrices/reference/reflect_zeroes.md).

If \\A,C\\ is given we can use `t(make_allpass(t(A), t(C)))`.

## Examples

``` r
m = 2
s = 6
A = matrix(rnorm(s*s), nrow = s, ncol = s)
B = matrix(rnorm(s*m), nrow = s, ncol = m)

K = make_allpass(A, B)
all.equal(cbind(A,B), cbind(K$A, K$B))
#> [1] TRUE
# check that K(z) is all-pass
Kf = zvalues(K)
all.equal(zvalues(polm(diag(m))) + complex(imaginary = 0), 
          Kf %r% Ht(Kf) + complex(imaginary = 0))
#> [1] TRUE

# if (A,C) is given, we proceed as follows:
C = matrix(stats::rnorm(m*s), nrow = m, ncol = s)
K = t(make_allpass(t(A), t(C)))
all.equal(rbind(A,C), rbind(K$A, K$C))
#> [1] TRUE
# check that K(z) is all-pass
Kf = zvalues(K)
all.equal(zvalues(polm(diag(m))) + complex(imaginary = 0), 
          Kf %r% Ht(Kf) + complex(imaginary = 0))
#> [1] TRUE
```
