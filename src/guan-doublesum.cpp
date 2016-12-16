// Unit vector version
#include <Rcpp.h>
#include <vector>
#include "kernels.h"

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector c_guan_doublesum_2d_box(NumericMatrix x,
                               NumericMatrix bbox,
                               NumericVector bw) {
  int i, j, l;
  double d2, w;

  double (*kernel)(double);
  kernel = &box_kernel;

  int bwn = bw.size();

  int n = x.nrow();
  int dim = x.ncol();
  NumericVector bw2 = bw*bw;
  double bw2max = max(bw2);

  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  NumericVector S(bwn);
  // main loop
  for(i=0; i < n-1; i++) {
    for(j=i+1; j < n; j++) {
      // distance^2
      d2=0;
      for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
      if(d2 < bw2max ){
        w = 1;
        for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
        for(l=0;l<bwn;l++) {
          if(d2 < bw2(l)) S(l)+=1/w;
        }
      }
    }
  }
  return S;
}
