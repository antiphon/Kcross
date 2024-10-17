#include <R.h>
#include "kernels.h"

double epa_kernel(double d){
  if(d < 0) return epa_kernel(-d);
  if(d > 1) return 0.0;
  return 0.75 * (1 - d*d);
}

double epa_cumu(double u){
  if(u < -1) return 0;
  if(u > 1) return 1;
  return  0.75 * (u - u*u*u/3.0 + 2.0/3.0);
}


double box_kernel(double d){
  if(d < -1) return 0.0;
  if(d > 1) return 0.0;
  return 0.5;
}


void print_i_n(int i, int n){ // only every 100th
  if(i % 100 == 0) Rprintf("           \r[%i/%i]", i, n);
}
void print_none(int i, int n){}

double gauss_kernel(double d){
  return exp(- 0.5 * d * d) / 2.50663;
}
