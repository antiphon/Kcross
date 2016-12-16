//
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//#define PI (3.141592653589793)
  
  // [[Rcpp::export]]
NumericMatrix c_K_cross_partial_box(NumericMatrix x,
                                    IntegerVector types,
                                    int ntypes,
                                    NumericVector counts,
                                    NumericVector intensities,
                                    NumericMatrix bbox, 
                                    NumericVector r) {
  int i, j, l, s, t, t1, t2;
  double d2, w;
  int nr;
  
  int n = x.nrow();
  int dim = x.ncol();
  NumericVector r2 = r*r;
  nr = r.length();
  
  NumericMatrix Nindiv(n, ntypes * nr);
  
  double rmax = r(nr-1);
  double rmax2 = rmax*rmax;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
  
  // main loop
  for(i=0; i < n-1; i++) {
    t1 = types[i]-1;
    for(j=i+1; j < n; j++) {
      t2 = types[j]-1;
      // distance^2
      d2=0;
      for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
      if(d2 < rmax2){
        // border correction
        w = 1;
        for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
        //
        for(t=0; t < nr; t++){
          if(d2 < r2[t]){
            Nindiv(i, t + t2 * nr) +=1/w;
            Nindiv(j, t + t1 * nr) +=1/w;
          }
        }
      }
    }
  }
  // products
  NumericMatrix Ks(nr*ntypes, nr*ntypes);
  for(i=0; i < n; i++) {
    t1 = types[i]-1;
    for(l = 0; l < ntypes; l++){
      for(t = 0; t < nr; t++) {
        for(s = 0; s < nr; s++) {
          w = sqrt( Nindiv(i, t + nr*l) * Nindiv(i, s + nr*l) ); 
          w = w/(intensities(l)*counts(t1));
          Ks(t + nr * t1, s + nr * l) += w;
        }
      }
    }
  }
  // compile result
  return Ks;
}

