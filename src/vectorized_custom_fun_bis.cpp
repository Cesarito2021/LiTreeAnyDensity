#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix vectorized_custom_fun_bis(NumericMatrix x, NumericMatrix y, int valx, int valy) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix out(nrow, ncol);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      out(i, j) = (x(i,j) == valx && y(i,j) == valy) ? x(i,j) + 1 : 0;
    }
  }
  return out;
}

// [[Rcpp::export]]
List overlay_func(List data, int number_layers, IntegerVector values_x, IntegerVector values_y, IntegerVector j1) {
  List output_list(number_layers);
  
  for (int i = 0; i < std::min(data.size(), values_x.size()); i++) {
    int valx = values_x[i];
    int valy = values_y[i];
    int j2 = j1[i] - 1;

    if (!Rf_isNull(data[i]) && !Rf_isNull(data[j2])) {
      NumericMatrix m1 = as<NumericMatrix>(data[i]);
      NumericMatrix m2 = as<NumericMatrix>(data[j2]);

      // Skip matrices that are empty (0 rows or cols)
      if (m1.nrow() == 0 || m1.ncol() == 0 || m2.nrow() == 0 || m2.ncol() == 0) {
        continue;
      }

      NumericMatrix result = vectorized_custom_fun_bis(m1, m2, valx, valy);

      if (Rf_isNull(output_list[i])) {
        output_list[i] = result;
      } else {
        NumericMatrix prev = as<NumericMatrix>(output_list[i]);
        output_list[i] = vectorized_custom_fun_bis(prev, m2, valx, valy);
      }
    }
  }

  return output_list;
}
