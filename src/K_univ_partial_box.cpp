// Unit vector version
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

#define PI (3.141592653589793)

// [[Rcpp::export]]
NumericMatrix c_K_partial_box(NumericMatrix x,
                     NumericMatrix bbox,
                     NumericVector r) {

  int i, j, l, s, t;
  double d2, w;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();
  NumericVector r2 = r*r;
  nr = r.length();


  NumericMatrix Nindiv(n,nr);

  double rmax = r(nr-1);
  double rmax2 = rmax*rmax;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
  double V = 1;
  for(i=0; i < dim; i++) V *= boxlen[i];
  double lambda = (n-1)/V;
  NumericVector pir2 = PI * r2 * lambda;

  // main loop
  for(i=0; i < n-1; i++) {
    for(j=i+1; j < n; j++) {
      // distance^2
      d2=0;
      for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
      if(d2 < rmax2 ){
        // border correction
        w = 1;
        for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
        for(t=0; t < nr; t++){
          if(d2 < r2[t]){
            Nindiv(i,t) +=1/w;
            Nindiv(j,t) +=1/w;
          }
        }
      }
    }
  }
  // products
  NumericMatrix Ks(nr, nr);
  for(i=0; i < n; i++) {
    for(t = 0; t < nr; t++) {
      for(s = 0; s < nr; s++) {
        w = 0.5 * ( (Nindiv(i,t)-pir2[t]) + (Nindiv(i,s)-pir2[s])  );
        Ks(t,s) += w;

      }
    }
  }

  // compile result
  return Ks;
}

