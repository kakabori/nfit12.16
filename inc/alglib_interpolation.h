#ifndef GUARD_ALGLIB_INTERPOLATION_H
#define GUARD_ALGLIB_INTERPOLATION_H

#include "alglib/src/stdafx.h"
#include "alglib/src/interpolation.h"
#include <vector>

class Alglib_CubicSpline2D {
public:
  void initialize(double (*f)(double, double, void*), void *p) {function = f; 
                                                                parameter = p;}
  void buildInterpolant(const std::vector<double>&, const std::vector<double>&,
                        const std::vector<double>&);
  void buildInterpolant(const std::vector<double>&, const std::vector<double>&);
  double evaluate(double, double);
private:
  double xMin, xMax, yMin, yMax;
  void *parameter;
  // The function to be interpolated takes two doubles for the independent 
  // variables and one void pointer that points to the parameters to evaluate 
  // the function
  double (*function)(double, double, void*); 
  alglib::real_1d_array xarray, yarray, farray;
  alglib::spline2dinterpolant spline;
};

#endif
