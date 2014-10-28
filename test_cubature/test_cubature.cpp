#include "cubature.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>

using namespace std;

int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
  double sigma = *((double *) fdata);
  double sum = 0;
  unsigned i;
  for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
  // compute the output value: note that fdim should == 1 from below
  fval[0] = exp(-sigma * sum);
  return 0;
}


int f2(unsigned ndim, unsigned npts, const double *x, void *fdata, unsigned fdim, double *fval) 
{
  double sigma = *((double *) fdata);
  unsigned i, j;
  //#pragma omp parallel for
  for (j = 0; j < npts; ++j) { // evaluate the integrand for npts points
    double sum = 0;
    for (i = 0; i < ndim; ++i) sum += x[j*ndim+i] * x[j*ndim+i];
    fval[j] = exp(-sigma * sum);
    int tmp;
    for (int k = 0; k < 1e20; k++) {
      tmp += k;
    }
  }
  return 0; // success
}


int main()
{
  double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5, val, err;
  //hcubature(1, f, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
  hcubature_v(1, f2, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
  printf("Computed integral = %0.10g +/- %g\n", val, err);  
  cout << "finished algMosaic: " << clock()*1e-6 << endl;
  return 0;
}
