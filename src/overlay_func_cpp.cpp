#include <Rcpp.h>
#include <Rinternals.h>
using namespace Rcpp;

// Forward declaration: implemented in vectorized_custom_fun_bis.cpp
NumericMatrix vectorized_custom_fun_bis(NumericMatrix x,
                                        NumericMatrix y,
                                        int valx, int valy);

// Scalar wrapper: call the matrix-based function with 1x1 matrices
inline double vectorized_custom_fun_bis_scalar(double x, double y, int valx, int valy) {
  NumericMatrix mx(1, 1), my(1, 1);
  mx(0, 0) = x;
  my(0, 0) = y;
  NumericMatrix out = vectorized_custom_fun_bis(mx, my, valx, valy);
  return out(0, 0);
}

// [[Rcpp::export]]
List overlay_func_cpp(List data,
                      int number_layers,
                      IntegerVector values_x,
                      IntegerVector values_y,
                      IntegerVector j1) {
  
  const int n_layers = std::min<int>(number_layers, data.size());
  List output_list(n_layers);
  
  const int n_idx = std::min({ (int)data.size(),
                             (int)values_x.size(),
                             (int)values_y.size(),
                             (int)j1.size(),
                             n_layers });
  
  for (int i = 0; i < n_idx; ++i) {
    const int valx = values_x[i];
    const int valy = values_y[i];
    
    // R is 1-based, C++ is 0-based
    const int j2 = j1[i] - 1;
    if (j2 < 0 || j2 >= data.size()) continue;
    
    // NULL checks: use Rf_isNull (proxy has no .isNULL())
    RObject xi = data[i];
    RObject yj = data[j2];
    if (Rf_isNull(xi) || Rf_isNull(yj)) continue;
    
    // Coerce to matrices (throws if incompatible)
    NumericMatrix mat1 = as<NumericMatrix>(xi);
    NumericMatrix mat2 = as<NumericMatrix>(yj);
    
    if (mat1.nrow() != mat2.nrow() || mat1.ncol() != mat2.ncol()) {
      // different shapes: skip safely
      continue;
    }
    
    const int nrow = mat1.nrow();
    const int ncol = mat1.ncol();
    
    NumericMatrix result(nrow, ncol);
    
    for (int r = 0; r < nrow; ++r) {
      for (int c = 0; c < ncol; ++c) {
        const double x = mat1(r, c);
        const double y = mat2(r, c);
        result(r, c) = vectorized_custom_fun_bis_scalar(x, y, valx, valy);
      }
    }
    
    if (Rf_isNull(output_list[i])) {
      output_list[i] = result;
    } else {
      NumericMatrix existing = as<NumericMatrix>(output_list[i]);
      if (existing.nrow() == nrow && existing.ncol() == ncol) {
        for (int r = 0; r < nrow; ++r) {
          for (int c = 0; c < ncol; ++c) {
            existing(r, c) = vectorized_custom_fun_bis_scalar(existing(r, c), mat2(r, c), valx, valy);
          }
        }
        output_list[i] = existing;
      } else {
        // shape mismatch – overwrite to keep consistent shape
        output_list[i] = result;
      }
    }
  }
  
  return output_list;
}
