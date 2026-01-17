// rationalmatrices_lyapunov.h
// Header-only implementation of the Lyapunov equation solver
// This allows direct C++ inclusion by dependent packages (e.g., RLDM)
// bypassing Rcpp's SEXP serialization which breaks pass-by-reference semantics.

#ifndef RATIONALMATRICES_LYAPUNOV_H
#define RATIONALMATRICES_LYAPUNOV_H

#include <RcppArmadillo.h>

namespace rationalmatrices {

// Core loop for solving Lyapunov equation via Schur method
// This helper function performs one iteration of the recursive solution algorithm
inline void solve_lyapunov_core_loop(int i, arma::cx_mat& cQ, const arma::cx_mat& cS, arma::cx_mat& cA) {
  // P22 = S22 * P22 * conj(S22) + Q22
  // Q22 <- P22 = Q22 / ( 1- S22*conj(S22)) )
  cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));

  // [Q11, Q12] <- [Q11, Q12] + S12 * P22 * [conj(S21), conj(S22)]
  cQ.submat(0,0,i-1,i) = cQ.submat(0,0,i-1,i) + ( cS.submat(0,i,i-1,i) * cQ(i,i) ) * cS.submat(0,i,i,i).t();

  // (I - S11 *conj(S22))
  cA.submat(0,0,i-1,i-1) = (-conj(cS(i,i))) * cS.submat(0,0,i-1,i-1);
  cA.submat(0,0,i-1,i-1).diag() += 1;

  // Q12 <- P12 = (I - S11 B22)^-1 * Q12
  cQ.submat(0,i,i-1,i) = solve(trimatu(cA.submat(0,0,i-1,i-1)), cQ.submat(0,i,i-1,i));
  // Q21 <- P21 = conj(P12)
  cQ.submat(i,0,i,i-1) = cQ.submat(0,i,i-1,i).t();

  // S11 * P12 * conj(S12)
  cA.submat(0,0,i-1,i-1) = ( cS.submat(0,0,i-1,i-1) * cQ.submat(0,i,i-1,i) ) * ( cS.submat(0,i,i-1,i).t() );

  // Q11 <- Q11 + S11 * P12 * conj(S12) + conj(S11 * P12 * conj(S12))
  cQ.submat(0,0,i-1,i-1) = cQ.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1).t();
}

// Solve Lyapunov equation P = A * P * A' + Q via Schur decomposition
// Returns true if A is stable (all eigenvalues inside unit circle), false otherwise
// If stop_if_non_stable is true and A is not stable, returns early without computing P
//
// Parameters:
//   A - system matrix
//   Q - noise covariance matrix
//   P - [output] solution matrix
//   lambda_r, lambda_i - [output] real and imaginary parts of eigenvalues of A
//   stop_if_non_stable - if true, return false early when A is not stable
inline bool lyapunov_impl(const arma::mat& A, const arma::mat& Q, arma::mat& P,
                          arma::vec& lambda_r, arma::vec& lambda_i, bool stop_if_non_stable) {

  int m = A.n_rows;

  arma::cx_mat cA = arma::cx_mat(A, arma::mat(m,m,arma::fill::zeros));
  arma::cx_mat cQ = arma::cx_mat(Q, arma::mat(m,m,arma::fill::zeros));
  arma::cx_mat cU = cQ;
  arma::cx_mat cS = cQ;

  // Schur decomposition: A = U*S*U' where S is upper triangular
  bool ok = arma::schur(cU, cS, cA);
  if (!ok) {
    Rcpp::stop("RcppArmadillo \"schur\" algorithm failed");
  }
  bool is_stable = (max(abs(cS.diag())) < 1);
  lambda_r = real(cS.diag());
  lambda_i = imag(cS.diag());

  if ( (stop_if_non_stable)  &&  (!is_stable) ) {
    return false;
  }

  // Transform Q to Schur basis
  cQ = cU.t() * cQ * cU;

  // Recursive solution: cQ is overwritten with the solution P
  int i;
  for (i = (m-1); i>0; i--) {
    solve_lyapunov_core_loop(i, cQ, cS, cA);
  }
  // For i=0, only the diagonal element update is needed (no submatrices exist)
  cQ(0,0) = cQ(0,0) / (1.0 - cS(0,0) * conj(cS(0,0)));

  // Transform back and ensure Hermitian
  cQ = cU * cQ * cU.t();
  cQ = (cQ + cQ.t())/2;
  P = real(cQ);

  return is_stable;
}

} // namespace rationalmatrices

#endif // RATIONALMATRICES_LYAPUNOV_H
