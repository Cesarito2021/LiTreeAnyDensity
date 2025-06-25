#include <RcppArmadillo.h>
#include <set>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Fit a circle using algebraic least squares
// [[Rcpp::export]]
NumericVector fitCircle(const NumericVector &x, const NumericVector &y) {
  int n = x.size();
  if (n < 3) stop("At least 3 points are required");
  
  arma::mat A(n, 3);
  arma::vec B(n);
  for (int i = 0; i < n; i++) {
    double xi = x[i], yi = y[i];
    A(i, 0) = xi;
    A(i, 1) = yi;
    A(i, 2) = 1.0;
    B[i] = -(xi * xi + yi * yi);
  }
  
  arma::vec sol;
  try {
    sol = arma::solve(A, B);
  } catch (...) {
    return NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  }
  
  double D = sol[0], E = sol[1], F = sol[2];
  double a = -D / 2.0, b = -E / 2.0;
  double r2 = a * a + b * b - F;
  double r = (r2 >= 0) ? std::sqrt(r2) : NA_REAL;
  return NumericVector::create(a, b, r);
}

// Helper: check if 3 points are colinear
bool isColinear(double x1, double y1, double x2, double y2, double x3, double y3) {
  double area = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
  return std::abs(area) < 1e-8;
}

// [[Rcpp::export]]
List rc_cpp(NumericMatrix data, int n, int k, double t, int d) {
  int N = data.nrow();
  NumericVector X = data(_, 0), Y = data(_, 1), Z = data(_, 2);
  NumericVector best_par = NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  double best_err = R_PosInf;
  
  for (int iter = 0; iter < k; iter++) {
    IntegerVector samp = sample(N, n, false);
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) idx[i] = samp[i] - 1;
    
    if (n == 3 && isColinear(X[idx[0]], Y[idx[0]], X[idx[1]], Y[idx[1]], X[idx[2]], Y[idx[2]]))
      continue;
    
    NumericVector sx(n), sy(n);
    for (int i = 0; i < n; i++) {
      sx[i] = X[idx[i]];
      sy[i] = Y[idx[i]];
    }
    
    NumericVector model = fitCircle(sx, sy);
    if (NumericVector::is_na(model[0])) continue;
    
    double a = model[0], b = model[1], r = model[2];
    NumericVector dist(N);
    for (int i = 0; i < N; i++) {
      dist[i] = std::sqrt((X[i] - a) * (X[i] - a) + (Y[i] - b) * (Y[i] - b));
    }
    
    std::vector<int> inliers;
    for (int i = 0; i < N; i++) {
      if (std::abs(dist[i] - r) < t) inliers.push_back(i);
    }
    
    if ((int)inliers.size() > d) {
      std::vector<int> all_idx = idx;
      for (int j : inliers) {
        if (std::find(all_idx.begin(), all_idx.end(), j) == all_idx.end())
          all_idx.push_back(j);
      }
      
      NumericVector ax(all_idx.size()), ay(all_idx.size());
      for (size_t i = 0; i < all_idx.size(); i++) {
        ax[i] = X[all_idx[i]];
        ay[i] = Y[all_idx[i]];
      }
      
      NumericVector refined = fitCircle(ax, ay);
      if (NumericVector::is_na(refined[0])) continue;
      
      double a2 = refined[0], b2 = refined[1], r2 = refined[2];
      double err = 0;
      for (int j : inliers) {
        double d2 = std::sqrt((X[j] - a2) * (X[j] - a2) + (Y[j] - b2) * (Y[j] - b2));
        err += std::abs(d2 - r2);
      }
      err /= inliers.size();
      
      if (err < best_err) {
        best_err = err;
        best_par = NumericVector::create(a2, b2, r2);
      }
    }
  }
  
  if (NumericVector::is_na(best_par[0])) return R_NilValue;
  double meanZ = mean(Z);
  return List::create(Named("par") = best_par, Named("height") = meanZ);
}

// [[Rcpp::export]]
List rmc_cpp(NumericMatrix data, int n, int k, double t, int d, int max_circles) {
  NumericMatrix rem = clone(data);
  List circles;
  
  while (rem.nrow() >= n && circles.size() < max_circles) {
    SEXP model = rc_cpp(rem, n, k, t, d);
    if (TYPEOF(model) == NILSXP) break;
    Rcpp::List modelList(model);
    NumericVector par = modelList["par"];
    
    double a = par[0], b = par[1], r = par[2];
    
    NumericVector dist(rem.nrow());
    for (int i = 0; i < rem.nrow(); i++) {
      double xi = rem(i, 0), yi = rem(i, 1);
      dist[i] = std::sqrt((xi - a) * (xi - a) + (yi - b) * (yi - b));
    }
    
    std::vector<int> inliers;
    for (int i = 0; i < rem.nrow(); i++) {
      if (std::abs(dist[i] - r) < t) inliers.push_back(i);
    }
    
    circles.push_back(model);
    
    if (!inliers.empty()) {
      std::set<int> to_remove(inliers.begin(), inliers.end());
      int new_n = rem.nrow() - to_remove.size();
      NumericMatrix new_rem(new_n, 3);
      int row_idx = 0;
      for (int i = 0; i < rem.nrow(); i++) {
        if (to_remove.find(i) == to_remove.end()) {
          new_rem(row_idx, 0) = rem(i, 0);
          new_rem(row_idx, 1) = rem(i, 1);
          new_rem(row_idx, 2) = rem(i, 2);
          row_idx++;
        }
      }
      rem = new_rem;
    }
  }
  
  return circles;
}
