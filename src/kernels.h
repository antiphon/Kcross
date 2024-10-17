#include <Rcpp.h>

#define PI (3.141592653589793)


double epa_kernel(double d);
double epa_cumu(double u);
double box_kernel(double d);

void print_i_n(int i, int n);
void print_none(int i, int n);
