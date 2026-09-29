// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
inline double corr_idx_fast(const arma::mat& A, arma::uword ia,
                            const arma::mat& B, arma::uword ib,
                            const arma::vec& alpha, const arma::vec& theta) {
  const arma::uword p = A.n_cols;
  double logc = 0.0;
  
  for (arma::uword k = 0; k < p; ++k) {
    const double d = A(ia, k) - B(ib, k);
    const double invtheta = 1.0 / (theta[k] * theta[k]);      // per-dim theta
    logc += -alpha[k] * std::log1p((d * d) * invtheta);      // per-dim alpha
  }
  return std::exp(logc);
}


// Return the correlation matrix R
// [[Rcpp::export]]
arma::mat R_fun_rcpp(const arma::mat& X,
                     const arma::vec& alpha, const arma::vec& theta)
{
  const arma::uword n = X.n_rows;
  arma::mat R;
  
  // Full build: n x n
  R.eye(n, n);
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const double c = corr_idx_fast(X, i, X, j, alpha, theta);
      R(i, j) = c;
      R(j, i) = c;
    }
  }
  return R;
  
  // // Nugget + Cholesky
  // arma::mat M = R;
  // M.diag() += 1e-6;
  // 
  // arma::mat U;
  // const bool ok = arma::chol(U, M); // upper
  // if (!ok) stop("Cholesky failed: matrix not positive definite (even with nugget).");
  // return U.t(); // lower
}


// Return r(x)
// [[Rcpp::export]]
arma::mat r_fun_rcpp(const arma::mat& X, const arma::mat& qmc,
                     const arma::vec& alpha, const arma::vec& theta) {
  const arma::uword n = qmc.n_rows;     // # test points
  const arma::uword m = X.n_rows;       // # design points 
  arma::mat new_R;
  new_R.set_size(n, m);
  
  for (arma::uword j = 0; j < m; ++j) {
    for (arma::uword i = 0; i < n; ++i) {
      new_R(i, j) = corr_idx_fast(qmc, i, X, j, alpha, theta);
    }
  }
  return new_R;
}

// nllh_fun for MIM kernel, optimize for theta and beta at the same time
// [[Rcpp::export]]
double nllh_fun(const arma::vec& param,
                const arma::mat& X, const arma::vec& y, 
                const int p,
                const double delta, const arma::vec& t) {
  const double beta = param(0);
  const arma::vec theta = param.tail(p);
  const arma::vec alpha = -arma::log(1-beta*t) / arma::log1p(delta*delta/arma::square(theta));
  arma::mat R = R_fun_rcpp(X, alpha, theta);
  const double n = R.n_rows;
  
  // cholesky decomposition
  arma::mat M = R;
  M.diag() += 1e-6;
  arma::mat U;
  const bool ok = arma::chol(U, M); // upper
  arma::mat L =  U.t(); // lower
  
  arma::vec ones = arma::ones<arma::vec>(n);
  arma::vec a = arma::solve(arma::trimatl(L), ones, arma::solve_opts::fast);
  arma::vec b = arma::solve(arma::trimatl(L), y, arma::solve_opts::fast);
  
  double mu = arma::dot(a,b) / arma::dot(a,a);
  double sigma2 = (arma::dot(b,b) - mu*mu*arma::dot(a,a)) / n;
  
  arma::vec d = L.diag();
  double nllh_1 = n*std::log(sigma2) + 2*arma::accu(arma::log(d));
  return nllh_1;
}


// nllh_fun for MIM kernel, optimize for theta and alpha at the same time
// [[Rcpp::export]]
double nllh_fun_alpha(const arma::vec& param,
                      const arma::mat& X, const arma::vec& y, 
                      const int p) {
  const arma::vec alpha = param.head(p);
  const arma::vec theta = param.tail(p);
  arma::mat R = R_fun_rcpp(X, alpha, theta);
  const double n = R.n_rows;
  
  // cholesky decomposition
  arma::mat M = R;
  M.diag() += 1e-12;
  arma::mat U;
  const bool ok = arma::chol(U, M); // upper
  arma::mat L =  U.t(); // lower
  
  arma::vec ones = arma::ones<arma::vec>(n);
  arma::vec a = arma::solve(arma::trimatl(L), ones, arma::solve_opts::fast);
  arma::vec b = arma::solve(arma::trimatl(L), y, arma::solve_opts::fast);
  
  double mu = arma::dot(a,b) / arma::dot(a,a);
  double sigma2 = (arma::dot(b,b) - mu*mu*arma::dot(a,a)) / n;
  
  arma::vec d = L.diag();
  double nllh_1 = n*std::log(sigma2) + 2*arma::accu(arma::log(d));
  return nllh_1;
} 


// [[Rcpp::export]]
Rcpp::List gp_predict(const arma::mat& X, const arma::vec& y, 
                     const arma::mat& new_X,
                     const arma::vec& alpha, const arma::vec& theta) {
  arma::mat R = R_fun_rcpp(X, alpha, theta);
  
  // cholesky decomposition
  arma::mat M = R;
  M.diag() += 1e-12;
  arma::mat U;
  const bool ok = arma::chol(U, M); // upper
  arma::mat L =  U.t(); // lower

  const arma::uword n = L.n_rows;
  arma::vec ones = arma::ones<arma::vec>(n);
  
  arma::vec a = arma::solve(arma::trimatl(L), ones, arma::solve_opts::fast);
  arma::vec b = arma::solve(arma::trimatl(L), y, arma::solve_opts::fast);
  double mu = arma::dot(a,b) / arma::dot(a,a);
  arma::vec c = b - mu*a;
  double sigma2 = arma::dot(c,c) / n;
  arma::vec d = arma::solve(arma::trimatu(L.t()), c, arma::solve_opts::fast);
  
  arma::mat new_R = r_fun_rcpp(X, new_X, alpha, theta);
  const arma::uword m = new_X.n_rows;
  arma::vec new_ones = arma::ones<arma::vec>(m);
  arma::vec pred_mean = mu*new_ones + new_R*d;

  arma::mat V = arma::solve(arma::trimatl(L), new_R.t(), arma::solve_opts::fast);
  arma::rowvec rRinvr = arma::sum(arma::square(V), 0);
  arma::vec pred_var = sigma2 * (1.0 - rRinvr.t());
  pred_var.elem(arma::find(pred_var <= 1e-6)).zeros();
  arma::vec pred_sd = arma::sqrt(pred_var);

  return Rcpp::List::create(
    Rcpp::Named("mean") = pred_mean,
    Rcpp::Named("sd") = pred_sd
  );
}


