#include <math.h>
#include <cmath>
#include <iostream>
#include <queue>
#include <vector>
#include <queue>
#include <algorithm>
#include <stdexcept>
#include "alglib/src/stdafx.h"
#include "alglib/src/interpolation.h"
#include "alglib_interpolation.h"

using namespace alglib;
using std::vector; using std::queue; using std::domain_error;
using std::cout; using std::endl;


/******************************************************************************
This function builds 2D cublic spline interpolant given two std::vectors
that contain x and y values at which the function must be computed.

It copies these vectors to arrays because alglib's setcontent function does 
not work with C++ vectors.

Once the interpolant is built, one can call evalute function to get 
an interpolated value at (x,y)
******************************************************************************/
void Alglib_CubicSpline2D::buildInterpolant(const vector<double>& xgrid, 
                                            const vector<double>& ygrid)
{
  vector<double> f;
  typedef vector<double>::size_type sz;
  
  for (sz i = 0; i < ygrid.size(); i++) {
    for (sz j = 0; j < xgrid.size(); j++) {
      f.push_back(function(xgrid[j], ygrid[i], parameter));
    }
  }
  
  buildInterpolant(xgrid, ygrid, f);
}


void Alglib_CubicSpline2D::buildInterpolant(const vector<double>& xgrid,
                                            const vector<double>& ygrid,
                                            const vector<double>& fvalue)
{  
  typedef vector<double>::size_type sz;
  sz sizex = xgrid.size();
  sz sizey = ygrid.size();
  sz sizef = fvalue.size();
  //cout << sizex << " "<< sizey << " " << sizef << endl;
  
  double _x[sizex];
  double _y[sizey];
  double _f[sizef];

  for (sz i = 0; i < sizex; i++) _x[i] = xgrid[i];
  for (sz i = 0; i < sizey; i++) _y[i] = ygrid[i];
  //for (sz i = 0; i < sizef; i++) _f[i] = fvalue[i];
  for (sz i = 0; i < sizef; i++) {
    _f[i] = fvalue[i];
    if (!std::isfinite(fvalue[i])) cout << i << " " << fvalue[i] << endl;
  }

  xarray.setcontent(sizex, _x);
  yarray.setcontent(sizey, _y);
  farray.setcontent(sizef, _f);
  
  try {
    //spline2dbuildbicubicv(xarray, sizex, yarray, sizey, farray, 1, spline);
    spline2dbuildbilinearv(xarray, sizex, yarray, sizey, farray, 1, spline);
  } catch (alglib::ap_error e) {
    cout << e.msg.c_str() << endl;
  }
  xMin = xgrid[0];
  xMax = xgrid.back();
  yMin = ygrid[0];
  yMax = ygrid.back();
}


double Alglib_CubicSpline2D::evaluate(double vx, double vy)
{
  if (vx > xMax || vx < xMin) {
    cout << "vx: " << vx << " xMin: " << xMin << " xMax: " << xMax << endl;
    throw domain_error("The requested x value is outside of the range.");
  } else if ( vy > yMax || vy < yMin) {
    cout << "vy: " << vy << " yMin: " << yMin << " yMax: " << yMax << endl;
    throw domain_error("The requested y value is outside of the range.");
  }
  
  return spline2dcalc(spline, vx, vy);
}

