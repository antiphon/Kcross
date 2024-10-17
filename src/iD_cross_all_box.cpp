#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List c_iD_cross_2d_box(NumericMatrix x,
                       NumericMatrix bbox,
                       NumericVector types,
                       int ntypes,
                       NumericVector intensities,
                       NumericVector min_intensities,
                       NumericVector bdist,
                       NumericVector r,
                       int correction) {

  int i, j, k, l, t1, t2, idx1, idx2;
  double d2, d;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();
  nr = r.length();

  std::vector<std::vector <double> > Ds(nr);
  std::vector<std::vector <double> > ww(nr);
  std::vector<std::vector <double> > nn(nr);
  IntegerVector ni(ntypes);

  for(i = 0; i < nr; i++) {
    Ds.at(i).resize(ntypes * ntypes);
    ww.at(i).resize(ntypes);
    for(j = 0; j < ntypes; j++) ww.at(i).at(j) = 1.0;
    nn.at(i).resize(ntypes);
  }
  double rmax = r(nr-1);
  double rmax2 = rmax * rmax;

  // main loop
  for(i=0; i < n; i++) {
    t1 = types(i) - 1;
    for(j=0; j < n; j++) {
      if(i != j){
        t2 = types(j) - 1;
        ni(t2)++;
        // distance^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        //
        if(d2 < rmax2){
          d = sqrt(d2);
          for(k = 0; k < nr; k++){
            if(bdist(i) > r(k) && d < r(k)){ // hit
              ww.at(k).at(t2) *= 1.0 - min_intensities(t2)/intensities(j);
            }
          }
        }
      }
    }
    // add to sums
    for(k = 0; k < nr; k++) {
      if(bdist(i) > r(k)){
        nn.at(k).at(t1) += 1.0/intensities(i); // Hamilton principle
        for(t2 = 0; t2 < ntypes; t2++) {
          if(ni(t2)>0){
            idx1 = t1 + ntypes * t2;
            Ds.at(k).at(idx1) += ww.at(k).at(t2) / intensities(i);
          }
        }
      }
    }
    for(k = 0; k < nr; k++)
    for(t2 = 0; t2 < ntypes; t2++) { // reset for next point
      ww.at(k).at(t2) = 1.0;
      ni(t2)=0;
    }

  }

  // averaging, weighted by intensities
  for(k = 0; k < nr; k++) {
    for(t1 = 0; t1 < ntypes; t1++)
        for(t2 = t1; t2 < ntypes; t2++){
          idx1 = t1 + ntypes * t2;
          idx2 = t2 + ntypes * t1;
          if(nn.at(k).at(t1)>0) Ds.at(k).at(idx1) /= nn.at(k).at(t1);
          if(t1<t2) if(nn.at(k).at(t2)>0) Ds.at(k).at(idx2) /= nn.at(k).at(t2);
        }
  }

  // compile result
  return wrap(Ds);
}



// [[Rcpp::export]]
List c_iD_cross_2d_box_from(NumericMatrix x,
                            NumericMatrix bbox,
                            NumericVector types,
                            int ntypes,
                            NumericVector intensities,
                            NumericVector min_intensities,
                            NumericVector bdist,
                            NumericVector r,
                            int correction,
                            int from) {

  int i, j, k, l, t1, t2, idx1, idx2;
  double d2, d;
  int nr;

  int n = x.nrow();
  int dim = x.ncol();
  nr = r.length();

  std::vector<std::vector <double> > Ds(nr);
  std::vector<std::vector <double> > ww(nr);
  std::vector< double> nn(nr);
  IntegerVector ni(ntypes);

  for(i = 0; i < nr; i++) {
    Ds.at(i).resize(ntypes);
    ww.at(i).resize(ntypes);
    for(j = 0; j < ntypes; j++) ww.at(i).at(j) = 1.0;
  }
  double rmax = r(nr-1);
  double rmax2 = rmax * rmax;

  // main loop
  for(i=0; i < n; i++) {
    t1 = types(i) - 1;
    if((t1 + 1) == from) {
    for(j=0; j < n; j++) {
      if(i != j){
        t2 = types(j) - 1;
        ni(t2)++;
        // distance^2
        d2=0;
        for(l=0; l < dim; l++)  d2 += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        //
        if(d2 < rmax2){
          d = sqrt(d2);
          for(k = 0; k < nr; k++){
            if(bdist(i) > r(k) && d < r(k)){ // hit
              ww.at(k).at(t2) *= 1.0 - min_intensities(t2)/intensities(j);
            }
          }
        }
      }
    }
    // add to sums
    for(k = 0; k < nr; k++) {
      if(bdist(i) > r(k)){
        nn.at(k) += 1.0/intensities(i); // Hamilton principle
        for(t2 = 0; t2 < ntypes; t2++) {
          if(ni(t2)>0){
            Ds.at(k).at(t2) += ww.at(k).at(t2) / intensities(i);
          }
        }
      }
    }
    for(k = 0; k < nr; k++)
      for(t2 = 0; t2 < ntypes; t2++) { // reset for next point
        ww.at(k).at(t2) = 1.0;
        ni(t2)=0;
      }
    }
  }
  // averaging, weighted by intensities
  for(k = 0; k < nr; k++) {
    for(t1 = 0; t1 < ntypes; t1++){
        if(nn.at(k)>0) Ds.at(k).at(t1) /= nn.at(k);
    }
  }

  // compile result
  return wrap(Ds);
}



