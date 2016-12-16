// Unit vector version
#include <Rcpp.h>
#include "kernels.h"
#include <vector>

using namespace Rcpp;

//#define PI (3.141592653589793)

// [[Rcpp::export]]
NumericVector c_pcf_biv_2d_box(NumericMatrix x, int n1,
                     NumericVector r,
                     double bw, NumericMatrix bbox, int correction) {
  int i, j, k, l;
  double d, d2, w;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();

  nr = r.length();
  NumericVector val(nr);
  double rmax = r(nr-1);
  double upper_limit = rmax + bw;
  upper_limit *= upper_limit;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  // int rm, rM;
  //
  // double dr = r(1)-r(0);
  //
  // int dri = (int)( 1.1*dr/bw );
  //
  // main loop
  for(i=0; i < n1; i++) {
      for(j=n1; j < n; j++)  {
        // distance^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        if(d2 < upper_limit){
          d = sqrt(d2);
          if(d > 0){
            // border correction
            w = 1;
            if(correction>0){
              for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
            }
            // add
            for(k = 0; k < nr; k++){
              val(k) += epa_kernel((r(k)-d)/bw)/(bw * w * d);
            }
          }
        }
      }
  }
  return val;
}




// Univariate

// [[Rcpp::export]]
NumericVector c_pcf_2d_box(NumericMatrix x,
                               NumericVector r,
                               double bw, NumericMatrix bbox, int correction,
                               int kern = 1) {
  int i, j, k, l;
  double d, d2, w;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();

  double (*kernel)(double);
  if(kern == 1) kernel = &epa_kernel;
  if(kern == 0) kernel = &box_kernel;

  nr = r.length();
  NumericVector val(nr);
  double rmax = r(nr-1);
  double upper_limit = rmax + bw;
  upper_limit *= upper_limit;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  // int rm, rM;
  //
  // double dr = r(1)-r(0);
  //
  // int dri = (int)( 1.1*dr/bw );
  //
  // main loop
  for(i=0; i < n; i++) {
    for(j=i+1; j < n; j++)  {
      // distance^2
      d2=0;
      for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
      if(d2 < upper_limit){
        d = sqrt(d2);
        if(d > 0){
          // border correction
          w = 1;
          if(correction>0){
            for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
          }
          // add
          for(k = 0; k < nr; k++){
            val(k) += kernel((r(k)-d)/bw)/(bw * w * d);
          }
        }
      }
    }
  }
  return val;
}
