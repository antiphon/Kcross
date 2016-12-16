#include "kernels.h"

double epa_kernel(double d){
  if(d < 0) return epa_kernel(-d);
  if(d > 1) return 0.0;
  return 0.75 * (1 - d*d);
}

double box_kernel(double d){
  if(d < -1) return 0.0;
  if(d > 1) return 0.0;
  return 0.5;
}
