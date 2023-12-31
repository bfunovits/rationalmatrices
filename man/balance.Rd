% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsp_methods.R
\name{balance}
\alias{balance}
\title{Balanced Realization and Balanced Truncation}
\usage{
balance(
  obj,
  gr,
  tol = 10 * sqrt(.Machine$double.eps),
  s0 = NULL,
  truncate = TRUE
)
}
\arguments{
\item{obj}{(\code{\link{stsp}} object) rational matrix in statespace form.}

\item{gr}{a list with two components \code{P} and \code{Q} which contain two Grammians. 
These Grammians may e.g. be computed with \code{\link{grammians}}.}

\item{tol}{a tolerance parameter used for the determination of the rank \eqn{PQ}. 
If \code{s0 = NULL} then the procedure estimates the rank of \eqn{PQ} 
and sets \code{s0} equal to this estimate. If \code{s0} is given, 
then \code{tol} is ignored.}

\item{s0}{determines the  size of the two diagonal blocks of transformed 
Grammians \eqn{P}, \eqn{Q}, respectively (for \code{truncate = TRUE} 
the statespace dimension of the balanced truncated statespace realization.}

\item{truncate}{(boolean) If true then the balanced truncated model is returned.}
}
\value{
A list with components
\item{obj}{(\code{\link{stsp}} object) represents the balanced (truncated) statespace 
           realization.}
\item{T,Tinv}{the state transformation and the inverse state transformation matrix. 
              Note that for the case \code{truncate=TRUE} the matrix \code{T} is 
              an \code{s0 x s} matrix and \code{Tinv} is of dimension \code{s x s0}.}
\item{P,Q}{the transformed (and truncated) Grammians.}
\item{sigma}{the vector of singular values of the matrix \eqn{(P^{1/2}Q^{1/2})}.}
}
\description{
Compute a balanced realization or a balanced truncated statespace realization
of a rational matrix in statespace form.
}
\details{
Let \eqn{P,Q} denote the controllability and the observability Grammian of a statespace 
realization 
\deqn{K(z) = C(Iz^{-1} - A)^{-1}B +D} 
of a rational matrix \eqn{K(z)} with statespace dimension \eqn{s}. 
The matrix \eqn{PQ} is diagonalizable and the eigenvalues are real and non negative. 
Let \eqn{\sigma_1^2 \geq \cdots \geq \sigma_s^2 \geq 0}{\sigma[1]^2 \ge \dots \ge \sigma[s]^2 \ge 0}
denote the eigenvalues of \eqn{PQ} and suppose that \eqn{PQ} has rank \eqn{k\leq s}{k \le s}, 
i.e. \eqn{\sigma_i=0}{\sigma[i]=0} for \eqn{i>k}. 
Then there exists a statespace transformation \eqn{T} which renders both Grammians into diagonal form 
and furthermore we have \eqn{p_{ii}=q_{ii}=\sigma_i}{p[i,i]=q[i,i]=\sigma[i]} for 
the first \eqn{k} diagonal entries of the (transformed) Grammians.

The (non zero) \eqn{\sigma_i}{\sigma[i]}'s are equal to the (non zero) singular values of the 
Hankel matrix of the impulse response coefficients of \eqn{K} and hence are called 
\emph{Hankel singular values}. These singular values are also the singular values of the 
product \eqn{P^{1/2}Q^{1/2}} where \eqn{P^{1/2}} and \eqn{Q^{1/2}} denote "square roots" 
of the Grammians \eqn{P} and \eqn{Q} respectively.
 
The procedure \code{balance(obj,...)} computes a somewhat simplified "balanced" form where 
the transformed Grammians \eqn{P,Q} are block diagonal with two diagonal blocks of dimension 
\eqn{s_0}{s[0]} and \eqn{s_1 = s - s_0}{s[1] = s - s[0]} respectively. The two 
upper, left blocks of \eqn{P} and \eqn{Q} are diagonal and the diagonal entries are  
\eqn{p_{ii}=q_{ii}=\sigma_i}{p[i,i]=q[i,i]=\sigma[i]} for  \eqn{i=1,\ldots,s_0}{i=1,...,s[0]}.
If \eqn{s_0=k}{s[0]=k} is equal to the rank of \eqn{PQ} 
then the product of the two lower, right blocks of \eqn{P} and \eqn{Q} is zero 
(up to numerical errors). 

Note that \eqn{k=\mathrm{rk}(PQ)}{k = rk(PQ)} is equal to the minimal statespace dimension, i.e. 
there exists a statespace realization with statespace dimension \eqn{k} and 
any statespace realization of \eqn{K} has a statespace dimension \eqn{geq k}{\ge k}. 
Such a \emph{minimal} statespace realization now be constructed by "truncating" the balanced 
realization:  
\deqn{K(z) = C_{1}(Iz^{-1} - A_{11})^{-1} B_{1} + D}{K(z) = C[1](Iz^{-1} - A[1,1])^{-1} B[1] + D}
where \eqn{A_{11}}{A[1,1]} is the left, upper, \eqn{(k,k)} dimensional block of \eqn{A}, 
\eqn{B_1}{B[1]} is the upper, \eqn{(k,n)} dimensional block of \eqn{B} and 
\eqn{C_1}{C[1]} is the left, \eqn{(m,k)} dimensional block of \eqn{C}. 

This balanced truncated statespace realization is returned if the 
optional parameter \code{truncate} is \code{TRUE}. Note that in the case where 
\code{s0} is less than the rank of \eqn{PQ} the truncated realization 
is just an approximate realization of \eqn{K(z)}. The approximation error 
depends on the size of "neglected" Hankel singular values. Note also that the 
truncated realization is not balanced (if \eqn{s_0<k}{s[0]<k}). 

If the "target" statespace dimension \eqn{s_0}{s[0]} is not given then the 
procedure sets \eqn{s_0}{s[0]} equal to an estimate of the rank of 
\eqn{PQ}. This estimate is computed by inspecting the singular values 
\eqn{\sigma_i}{\sigma_[i]} of the product of the square roots 
\eqn{P^{1/2}} and \eqn{Q^{1/2}}. 

The above discussion deals with the controllability and the observability 
Grammian of the statespace realization. However, one may also use 
other pairs of Grammians, e.g. for invertible matrices \eqn{K} 
one may use the controllabaility Grammian of the statespace realizatton 
and the observability Grammian of the statespace realization of the 
inverse \code{K^{-1}(z)}.
}
\examples{
# example A ############################################################

# "obj" is a (1 by 1) rational matrix in statespace form, 
# with stespace dimension s = 2. 

obj = stsp(A = c(0,0.2,1,-0.5),
           B = c(1,1), C = c(1,0))
gr = grammians(obj, 'lyapunov')
bal = balance(obj, gr)

print(cbind(bal$P, bal$Q, diag(bal$sigma, nrow = 2, ncol = 2)))
all.equal(grammians(bal$obj), bal[c('P','Q')])

# example B (non minimal statespace realization #########################

# The "rbind" operation below returns a statespace realization with 
# statespace dimension s = 4. However the minimal statespace dimensions 
# is s0 = 2. 
obj = rbind(obj, obj)
gr = grammians(obj, 'lyapunov')
bal = balance(obj, gr, s0 = 2, truncate = FALSE)

# the upper (2 by 2) block of the (transformed) controllability 
# Grammian is diagonal, the lower (2 by 2) block is "zero". 
# This shows that the (balanced) realization is not controllable. 
print(bal$P)  
              
# the upper (2 by 2) block of the (transformed) observability 
# Grammian is diagonal and equal to the upper block of bal$P.
print(bal$Q) 

# the product of the (transformed) controllability and observability 
# Grammians is (approximately) diagonal and the diagonal entries are  
# the squares of the Hankel singular values.
print(bal$P \%*\% bal$Q)
print(bal$sigma^2)    

all.equal(grammians(bal$obj), bal[c('P','Q')])

# we may construct a minimal realization by 'balanced truncation'.
# note that we let the procedure determine the minimal statespace dimension 
trunc = balance(obj, gr)  
print(trunc$obj)
# compare with the above balanced realization 
print(bal$obj)
# check 
all.equal(pseries(obj), pseries(trunc$obj))

# example C (balanced truncation) ##########################

# construct a random rational matrix with statespace dimension s=10
obj = test_stsp(dim = c(2,2), s = 10, bpoles = 1, bzeroes = 1)
# compute an approximate realization with s0 = 8
gr = grammians(obj, 'minimum phase')
trunc = balance(obj, gr, s0 = 5)
print(trunc$sigma)

max(abs(unclass(pseries(obj, lag.max = 25)) - 
        unclass(pseries(trunc$ob, lag.max = 25))))
plot(pseries(obj, lag.max = 25), x_list= list(pseries(trunc$obj, lag.max = 25)), 
     type = c('l','p'), legend = c('s=10', 's=5')) 
}
