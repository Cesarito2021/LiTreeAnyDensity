#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List overlay_func_cpp(List data, 
                      int number_layers, 
                      IntegerVector values_x, 
                      IntegerVector values_y, 
                      IntegerVector j1) {
  
  List output_list(number_layers);
  
  for (int i = 0; i < std::min(data.size(), values_x.size()); i++) {
    int valx = values_x[i];
    int valy = values_y[i];
    int j2   = j1[i] - 1; // R is 1-based, C++ is 0-based
    
    if (!data[i].isNULL() && !data[j2].isNULL()) {
      NumericMatrix mat1 = data[i];
      NumericMatrix mat2 = data[j2];
      
      int nrow = mat1.nrow();
      int ncol = mat1.ncol();
      
      NumericMatrix result(nrow, ncol);
      
      for (int r = 0; r < nrow; r++) {
        for (int c = 0; c < ncol; c++) {
          double x = mat1(r, c);
          double y = mat2(r, c);
          result(r, c) = vectorized_custom_fun_bis(x, y, valx, valy);  // già definita altrove
        }
      }
      
      if (output_list[i].isNULL()) {
        output_list[i] = result;
      } else {
        NumericMatrix existing = output_list[i];
        for (int r = 0; r < nrow; r++) {
          for (int c = 0; c < ncol; c++) {
            existing(r, c) = vectorized_custom_fun_bis(existing(r, c), mat2(r, c), valx, valy);
          }
        }
        output_list[i] = existing;
      }
    }
  }
  
  return output_list;
}
