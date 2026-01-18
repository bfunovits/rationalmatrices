# Munkres Assignment Algorithm

This algorithm solves the *assignment problem*: Given a set of \\m\\
"agents" and \\n\\ "tasks" with agent-task specific costs. Find an
assignment of tasks to agents such that the total cost is *minimal*,
given the following restrictions. A task may be assigned to at most one
agent and an agent may be assigned to at most one task. In total
\\k=\min(m,n)\\ assignments are required. (E.g. if there are less agents
than tasks then each agent gets a "job", however, \\n-m\\ tasks remain
unaccomplished. )

## Usage

``` r
munkres(C)
```

## Arguments

- C:

  (\\(m,n)\\ numeric matrix) \\C\_{ij}\\ represents the cost for
  assigning the \\j\\-th "job" to the \\i\\-th "agent".

## Value

List with slots

- a:

  (\\(k,2)\\ dimensional (integer) matrix where \\k=\min(m,n)\\. This
  matrix represents the optimal assignment. For each \\i=1,\ldots,k\\
  the task `a[i,2]` is assigned to agent `a[i,1]`.

- c:

  Total cost for this (optimal) assignment
  `c = C[a[1,1], a[1,2]] + ... + C[a[k,1], a[k,2]]`.

- C:

  The cost matrix.

## Details

The original references are (Munkres 1956) and (Bourgeois and Lassalle
1971) for the non-square case.

## References

Munkres J (1956). “Algorithms for Assignment and Transportation
Problems.” *Journal of the Society for Industrial and Applied
Mathematics*, **5**(1).

Bourgeois F, Lassalle J (1971). “An extension of the Munkres algorithm
for the assignment problem to rectangular matrices.” *Commun. ACM*,
**14**, 802-804.

## See also

This helper function is mainly needed to match poles (zeroes).

## Examples

``` r
C = matrix(c (0, 0, 0, 0, 0, 1, 3, 3, 0, 5, 5, 9, 0, 1, 3, 7), 
           nrow = 4, ncol = 4, byrow = TRUE)
out = munkres(C)
print(out)
#> $a
#>      [,1] [,2]
#> [1,]    1    4
#> [2,]    2    3
#> [3,]    3    1
#> [4,]    4    2
#> 
#> $c
#> [1] 4
#> 
#> $C
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    1    3    3
#> [3,]    0    5    5    9
#> [4,]    0    1    3    7
#> 
```
