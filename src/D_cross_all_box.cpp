#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//#define PI (3.141592653589793)

// [[Rcpp::export]]
List c_D_cross_2d_box(NumericMatrix x, 
                      NumericMatrix bbox, 
                      NumericVector types,
                      NumericVector intensities,
                      NumericVector r) {
  
  int i, j, l, t1, t2, idx1;
  double d2, w;
  int nr;
  
  int n = x.nrow();
  int ntypes = intensities.length();
  int dim = x.ncol();
  
  nr = r.length();
  NumericVector r2 = r*r;
  std::vector<std::vector <double> > Ds(nr);
  
  std::vector<double> nnd(ntypes);
  for(i = 0; i < nr; i++) Ds.at(i).resize(ntypes * ntypes);
  
  NumericMatrix counts(ntypes, ntypes);
  
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
  
  double rmax = r(nr-1);
  double rmax2 = rmax*rmax;
  
  // main loop
  for(i=0; i < n; i++) {
    t1 = types(i) - 1;
    // distance to border
    w = boxlen(0);
    for(l=0; l < dim; l++) {
      w = fmin(w, fmin(bbox(1,l)-x(i,l), x(i,l)-bbox(0,l)));
    }
    w *= w; //distances are squared
    // nearest neighbours per type
    fill(nnd.begin(), nnd.end(), rmax2 + 1);
    for(j=0; j < n; j++){
      if(i!=j){
        t2 = types(j) - 1;
        // dist^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        if(nnd.at(t2) > d2) nnd.at(t2) = d2;
      }
    }
    // then we add the contributions 
    for(t2=0; t2 < ntypes; t2++) {
      if(nnd.at(t2) < w){
        counts(t1, t2) += 1.0;
        idx1 = t1 + t2 * ntypes;
        for(l=0; l < nr; l++) {
          if(nnd.at(t2) < r2(l)) Ds.at(l).at(idx1) += 1.0;
        }
      }
    }
  }
  // average
  for(t1 = 0; t1 < ntypes; t1++){
    for(t2 = t1; t2 < ntypes; t2++){
      for(l = 0; l < nr; l++){
        if(counts(t1,t2)>0) Ds.at(l).at(t1+t2*ntypes) /= counts(t1,t2);
        if(t1 != t2) if(counts(t2,t1)>0) Ds.at(l).at(t2+t1*ntypes) /= counts(t2,t1);
      }
    }
  }
  // compile result
  return wrap(Ds);
}

