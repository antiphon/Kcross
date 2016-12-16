// Unit vector version
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//#define PI (3.141592653589793)

// [[Rcpp::export]]
List c_K_cross_2d_box(NumericMatrix x,
                        NumericMatrix bbox,
                        NumericVector types,
                        int ntypes,
                        NumericVector r,
                        int correction) {

  int i, j, k, l, t1, t2, idx1, idx2;
  double d2, w;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();
  NumericVector r2 = r*r;
  nr = r.length();

  std::vector<std::vector <double> > Ks(nr);

  for(i = 0; i < nr; i++) Ks.at(i).resize(ntypes * ntypes);

  double rmax = r(nr-1);
  double rmax2 = rmax*rmax;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++){
    boxlen[i] = bbox(1,i) - bbox(0,i);
  }

  // main loop
  for(i=0; i < n-1; i++) {
    t1 = types(i) - 1;
      for(j=i+1; j < n; j++) {
        t2 = types(j) - 1;
        // distance^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        if(d2 < rmax2 ){
          // border correction
          w = 1;
          if(correction>0){
            for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
          }
          // add to correct position:
          idx1 = t1 + ntypes * t2;
          idx2 = t2 + ntypes * t1;
          for(k = nr-1; k > -1; k--){
            if(d2 > r2.at(k)) break;
              Ks.at(k).at(idx1) += 1/w;
              Ks.at(k).at(idx2) += 1/w;
          }
        }
      }
  }

  // compile result
  return wrap(Ks);
}





// [[Rcpp::export]]
NumericVector c_K_2d_box(NumericMatrix x,
                      NumericMatrix bbox,
                      NumericVector r,
                      int correction) {

  int i, j, k, l;
  double d2, w;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();
  NumericVector r2 = r*r;
  nr = r.length();

  NumericVector Ks(nr);


  double rmax = r(nr-1);
  double rmax2 = rmax*rmax;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);

  // main loop
  for(i=0; i < n-1; i++) {
    for(j=i+1; j < n; j++) {
      // distance^2
      d2=0;
      for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
      if(d2 < rmax2 ){
        // border correction
        w = 1;
        if(correction>0){
          for(l=0; l < dim; l++) w *= boxlen(l) - fabs(x(i,l)-x(j,l));
        }
        for(k = nr-1; k > -1; k--){
          if(d2 > r2.at(k)) break;
          Ks(k) += 1/w;
        }
      }
    }
  }

  // compile result
  return Ks;
}


// [[Rcpp::export]]
List c_K_cross_2d_box_just_one_from(NumericMatrix x,
                                    NumericMatrix bbox,
                                    NumericVector types,
                                    int from,
                                    int ntypes,
                                    NumericVector r,
                                    int correction) {

  int i, j, k, l, t1, t2, idx1, idx2;
  double d2, w;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();
  NumericVector r2 = r*r;
  nr = r.length();

  std::vector<std::vector <double> > Ks(nr);

  for(i = 0; i < nr; i++) Ks.at(i).resize(ntypes * ntypes);

  double rmax = r(nr-1);
  double rmax2 = rmax*rmax;
  NumericVector boxlen(dim);
  for(i=0; i < dim; i++){
    boxlen[i] = bbox(1,i) - bbox(0,i);
  }

  // main loop
  bool t1ok;
  for(i=0; i < n-1; i++) {
      t1ok = (types(i) == from);
      for(j=i+1; j < n; j++) {
        if( t1ok | (types(j) == from) ){
          t1 = types(i) - 1;
          t2 = types(j) - 1;
          // distance^2
          d2=0;
          for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
          if(d2 < rmax2 ){
            // border correction
            w = 1;
            if(correction>0){
              for(l=0; l < dim; l++) w *= boxlen(l) - std::abs(x(i,l)-x(j,l));
            }
            // add to correct position:
            idx1 = t1 + ntypes * t2;
            idx2 = t2 + ntypes * t1;
            for(k = nr-1; k > -1; k--){
              if(d2 > r2.at(k)) break;
              Ks.at(k).at(idx1) += 1/w;
              Ks.at(k).at(idx2) += 1/w;
            }
          }
        }
      }
  }

  // compile result
  return wrap(Ks);
}

