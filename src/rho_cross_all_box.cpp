// Unit vector version
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
List c_rho_cross_2d_box(NumericMatrix x,
                        NumericMatrix bbox,
                        NumericVector types,
                        NumericVector intensities,
                        NumericVector r,
                        double adjust,
                        int correction,
                        bool kernel_correction,
                        int verb = 0 ) {

  int i, j, k, l, t1, t2, idx1, idx2;
  double d, d2, w;
  int nr;

  int n = x.nrow();
  int ntypes = intensities.length();
  int dim = x.ncol();


  // bandwidths
  NumericVector bws = adjust * 0.15 / sqrt(intensities);

  nr = r.length();

  std::vector<std::vector <double> > pcfs(nr);

  for(i = 0; i < nr; i++) pcfs.at(i).resize(ntypes * ntypes);

  void (*verbose)(int i, int n) = print_none;
  if(verb > 0) verbose = &print_i_n;

  double rmax = r(nr-1);
  double dr;
  NumericVector upper_limit = (rmax + bws) * (rmax + bws);

  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  NumericMatrix kern_edge_weight(r.size(), ntypes);
  for(i =0; i < ntypes; i ++)
  for(k = 0; k < nr; k++) {
    if(kernel_correction){
      w = r(k)/bws(i);
      kern_edge_weight(k, i) = epa_cumu(w);
    }
    else{
      kern_edge_weight(k, i) = 1;
    }
  }



  // main loop
  for(i=0; i < n-1; i++) {
    t1 = types(i) - 1;
    verbose(i+1, n);
      for(j=i+1; j < n; j++) {
        t2 = types(j) - 1;
        // distance^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        if(d2 < fmax(upper_limit(t1), upper_limit(t2)) ){
          d = sqrt(d2);
          if(d > 0){
            // border correction
            w = 1;
            if(correction>0){
              for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
            }
            // add to correct position:
            idx1 = t1 + ntypes * t2;
            idx2 = t2 + ntypes * t1;
            for(k = 0; k < nr; k++){
              dr = d;//fmax(d, r(k));
              pcfs.at(k).at(idx1) += epa_kernel((r(k)-d)/bws(t2))/(bws(t2) * w * dr * kern_edge_weight(k, t2));
              pcfs.at(k).at(idx2) += epa_kernel((r(k)-d)/bws(t1))/(bws(t1) * w * dr * kern_edge_weight(k, t1));
            }
          }
        }
      }
  }
  if(verb > 0) Rprintf("\n");
  // compile result
  return wrap(pcfs);
}



// [[Rcpp::export]]
NumericVector c_rho_2d_box(NumericMatrix x,
                           NumericMatrix bbox,
                           NumericVector r,
                           double bw,
                           int correction,
                           int kernel_correction,
                           int kern = 1) {

  int i, j, k, l;
  double d, d2, w;
  double dr;
  int nr;

  double (*kernel)(double);
  if(kern == 1) kernel = &epa_kernel;
  if(kern == 0) kernel = &box_kernel;


  int n = x.nrow();
  int dim = x.ncol();

  nr = r.length();

  NumericVector pcf(nr);

  double rmax = r(nr-1);
  double upper_limit = (rmax + bw) * (rmax + bw);

  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  NumericVector kern_edge_weight(r.size());
  for(k = 0; k < nr; k++) {
      if(kernel_correction){
        w = r(k)/bw;
        kern_edge_weight(k) = epa_cumu(w);
      }
      else{
        kern_edge_weight(k) = 1;
      }
  }

  // main loop
  for(i=0; i < n-1; i++) {
    for(j=i+1; j < n; j++) {
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
          w *= bw;
          for(k = 0; k < nr; k++){
            dr = d;//fmax(d, r(k));
            pcf.at(k) += kernel((r(k)-d)/bw) / ( w * dr * kern_edge_weight(k));
          }
        }
      }
    }
  }
  return pcf;
}
