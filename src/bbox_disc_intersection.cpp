#include <Rcpp.h>
using namespace Rcpp;

void f(double d1, double d2, double r, double *ra){
  double r2 = r*r, d3, w;
  if(d1 < r) {
    d3 = sqrt(r2 - d1*d1);
    if(d3 < d2) {
      w = atan2(d2,d1) - atan2(d3, d1);
      ra[0] += 0.5 * (d1 * d3 + r2 * w);
      ra[1] += w*r;
    }
    else ra[0] += 0.5 * d1 * d2;
  }
  else{
    w = atan2(d2,d1);
    ra[0] += 0.5 * r2 * w;
    ra[1] += w*r;
  }
}


// [[Rcpp::export]]
NumericMatrix c_bbox_disc_intersection(NumericMatrix x,
                                     NumericMatrix bbox,
                                     double r) {
  NumericMatrix out(x.nrow(),2);
  int i;
  double a0 = bbox(0,0), a1 = bbox(1,0);
  double b0 = bbox(0,1), b1 = bbox(1,1);
  double ra[2];

  for(i = 0; i < x.nrow(); i++) {
    ra[0] = 0;
    ra[1] = 0;
    f(a1-a0-x(i,0),b1-b0-x(i,1),r,ra);
    f(a1-a0-x(i,0),x(i,1)-b0,r,ra);
    f(x(i,1)-b0,a1-a0-x(i,0),r,ra);
    f(x(i,1)-b0,x(i,0)-b0,r,ra);
    f(x(i,0)-b0,x(i,1)-b0,r,ra);
    f(x(i,0)-b0,b1-b0-x(i,1),r,ra);
    f(b1-b0-x(i,1),x(i,0)-b0,r,ra);
    f(b1-b0-x(i,1),a1-a0-x(i,0),r,ra);
    out(i,0) = ra[0];
    out(i,1) = ra[1];
  }
  return out;
}
