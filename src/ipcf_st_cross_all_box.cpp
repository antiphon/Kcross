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

// V <- c_ipcf_st_cross_2d_box(xy, bbox, ntypes, types,
//                             int,
//                             r, t,
//                             sigmas,
//                             as.integer( do_correction) )
//

// [[Rcpp::export]]
List c_ipcf_st_cross_2d_box(NumericMatrix x,
                            NumericMatrix bbox,
                            int ntypes,
                            NumericVector types,
                            NumericVector intensities,
                            NumericVector r,
                            NumericVector t,
                            NumericVector sigmas,
                            int correction) {

  int i, j, k, l, m1, m2, idx1, idx2;
  double d, d2, w;

  int nr = r.length();
  int nt = t.length();
  int n = x.nrow();
  int dim = x.ncol();
  // time should dim-1'th, mark dim'th

  // result object
  std::vector<std::vector <double> > pcfs(nr * nt);
  for(i = 0; i < nr*nt; i++) pcfs.at(i).resize(ntypes * ntypes);

  double bwx = sigmas(0);
  double bwt = sigmas(1);
  double dr, dt;
  double kr, kt;
  double vr, vt; // contributions
  double d2_upper_limit = pow(max(r) + 2*bwx, 2);
  double dt_upper_limit = max(t) + 2*bwt;
  int idx_t = dim - 2; // which column is time
  int dimx  = dim - 2; // number of space dimensions

  double bd = pow(PI, dimx/2.0) / R::gammafn(dimx/2.0+1.0) * 2.0; // last div from i<j
  //Rprintf("%f\n", kd);

  double (*kernel)(double) = epa_kernel;

  NumericVector boxlen(dimx+1);
  for(i=0; i < dimx+1; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
  // main loop
  for(i=0; i < n; i++) {
    m1 = types(i) - 1;
    for(j=i+1; j < n; j++) {
      // temporal distance
      dt = fabs(x(i, idx_t) - x(j, idx_t));
      if(dt < dt_upper_limit)
      {
        // spatial distance^2
        d2 = 0;
        for(l=0; l < dimx; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        if(d2 < d2_upper_limit ){ // this ij pair contributes
          m2 = types(j) - 1;
          d = sqrt(d2);
          if(d > 0){ // avoid artifacts
            // border correction?
            w = 1;
            if(correction>0){
              for(l=0; l < dim-1; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
            }
            // add to correct position in the result array:
            idx1 = m1 + ntypes * m2;
            idx2 = m2 + ntypes * m1;
            w *= intensities(i) * intensities(j) * bd;
            for(k = 0; k < nr; k++){
              dr = fmax(d, r(k));
              kr = r(k) - d;
              vr = kernel(kr / bwx);
              if(vr > 0) {
                vr /= bwx * pow(dr, dimx-1) ;
                for(l = 0; l < nt; l++) {
                  kt = t(l) - dt;
                  vt = kernel(kt / bwt) / bwt;
                  dr = vt * vr / w;
                  pcfs.at(k * nt + l).at(idx1) += dr;
                  pcfs.at(k * nt + l).at(idx2) += dr;
                }
              }
            }
          }
        }
      }
    }
  }
  // Rprintf("dd ");
  // compile result
  return wrap(pcfs);
}
