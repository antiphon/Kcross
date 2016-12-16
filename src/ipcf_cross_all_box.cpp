#include <Rcpp.h>
#include <vector>
#include "kernels.h"

using namespace Rcpp;

//#define PI (3.141592653589793)


// double epa_kernel(double d){
//   if(d < 0) d = -d;
//   if(d > 1) return 0.0;
//   return 0.75 * (1 - d*d);
// }



// [[Rcpp::export]]
List c_ipcf_cross_2d_box(NumericMatrix x,
                        NumericMatrix bbox,
                        int ntypes,
                        NumericVector types,
                        NumericVector intensities,
                        NumericVector r,
                        double adjust,
                        int correction) {

  int i, j, k, l, t1, t2, idx1, idx2;
  double d, d2, w;
  double bw1, bw2;

  int nr = r.length();
  int n = x.nrow();
  int dim = x.ncol();

  std::vector<std::vector <double> > pcfs(nr);

  for(i = 0; i < nr; i++) pcfs.at(i).resize(ntypes * ntypes);

  // max lambda
  double lmin = min(intensities);
  double bwf = 0.15 * adjust;
  double bwm = bwf / sqrt(lmin);
  double rmax = r(nr-1);
  double dr;
  double kr;
  double upper_limit = (rmax + bwm) * (rmax + bwm);

  NumericVector sqintensities = sqrt(intensities);

  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  // main loop
  for(i=0; i < n; i++) {
    t1 = types(i) - 1;
      for(j=i+1; j < n; j++) {
        //if(i != j)
          {
          t2 = types(j) - 1;
          // distance^2
          d2=0;
          for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
          if(d2 < upper_limit ){
            d = sqrt(d2);
            if(d > 0){
              bw1 = bwf / sqintensities(i);
              bw2 = bwf / sqintensities(j);
              //if(d < fmin(bw1+rmax, bw2+rmax))
                {
                // border correction
                w = 1;
                if(correction>0){
                  for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
                }
                // add to correct position:
                idx1 = t1 + ntypes * t2;
                idx2 = t2 + ntypes * t1;
                w = w * intensities(i) * intensities(j) * PI;
                for(k = 0; k < nr; k++){
                  dr = fmax(d, r(k));
                  kr = r(k)-d;
                  pcfs.at(k).at(idx1) += epa_kernel(kr/bw2) / ( bw2 * w * dr );
                  pcfs.at(k).at(idx2) += epa_kernel(kr/bw1) / ( bw1 * w * dr );
                }
              }
            }
          }
        }
      }
  }

  // compile result
  return wrap(pcfs);
}



// [[Rcpp::export]]
NumericVector c_ipcf_2d_box(NumericMatrix x,
                            NumericMatrix bbox,
                            NumericVector intensities,
                            double bw,
                            NumericVector r,
                            double adjust,
                            int correction) {

  int i, j, k, l;
  double d, d2, w;
  double dr;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();

  // bandwidth
  nr = r.length();
  NumericVector pcf(nr);

  double bwf = 0.15 * adjust;
  double rmax = r(nr-1);
  double upper_limit = (rmax + bw) * (rmax + bw);

  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  // main loop
  for(i=0; i < n; i++) {
    for(j=0; j < n; j++) {
      if(i != j){
        // distance^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        if(d2 < upper_limit ){
          d = sqrt(d2);
          if(d > 0){
            // border correction
            w = 1;
            if(correction>0){
              for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
            }
            bw = bwf / sqrt(intensities(j));
            w *= bw * intensities(i) * intensities(j) * PI;
            for(k = 0; k < nr; k++){
              dr = fmax(d, r(k));
              pcf.at(k) += epa_kernel((r(k)-d)/bw) / ( w * dr );
            }
          }
        }
      }
    }
  }
  return pcf;
}
