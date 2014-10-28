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








/**************************************************************************************************
Check whether left, middle, and right points form a straight line within specified tolerances
**************************************************************************************************/
bool Alglib_CubicSpline2D::onStraightLine(double dx, double dy, double y2)
{
  if( dx > maxdx) return false;
  return ( dx < mindx || dy < abserr || dy < fabs(y2*relerr) );
}

/**************************************************************************************************
Determines the grid points at which the function will get evaluated
Returns the number of grid points contained in oneDimGrid object
  initLeft -> smallest grid point
  initRight -> largest grid point
  line -> the value at which the grid line gets built
  oneDimGrid -> a list of grid points along one direction
**************************************************************************************************/
int Alglib_CubicSpline2D::buildGridPoints(double initLeft, double initRight, double line, vector<double>& oneDimGrid)
{
  double left, right, vleft, vright;
  double middle, vmiddle;
  queue<double> recursionQueue;

  double initVleft = function(initLeft, line, parameter);
  double initVright = function(initRight, line, parameter);

  oneDimGrid.push_back(initLeft);

  /*Start by setting up a queue with the initial right and left end points and
    their corresponding function values.
  */
  recursionQueue.push(initLeft);
  recursionQueue.push(initRight);
  recursionQueue.push(initVleft);
  recursionQueue.push(initVright);

  while (! recursionQueue.empty()) {
    left = recursionQueue.front();
    recursionQueue.pop();
    right = recursionQueue.front();
    recursionQueue.pop();
    vleft = recursionQueue.front();
    recursionQueue.pop();
    vright = recursionQueue.front();
    recursionQueue.pop();
    middle = (left + right) / 2;
    vmiddle = function(middle, line, parameter);

    /*If the function is not linear between two end points within
      specified tolerances, split the range in half, and add two sets of
      points and function values to the back of the queue.
    */
    if ( !onStraightLine(middle-left, fabs(vleft+vright-vmiddle-vmiddle) , vmiddle)) {
      recursionQueue.push(left);
      recursionQueue.push(middle);
      recursionQueue.push(vleft);
      recursionQueue.push(vmiddle);

      recursionQueue.push(middle);
      recursionQueue.push(right);
      recursionQueue.push(vmiddle);
      recursionQueue.push(vright);

    /*If the function is linear within tolerances over the given range, add the
      middle and right points to the list of interpolation points. Moreover,
      add the function values corresponding to the middle and right points to
      the list of calculated points.
    */
    } else {
      oneDimGrid.push_back(middle);
      oneDimGrid.push_back(right);
    }

  
  }

  /*The points are not guaranteed to be in sorted order, so x values must be
    sorted. Finally corresponding function values are added to the vector of y
  values.*/
  sort(oneDimGrid.begin(), oneDimGrid.end());
    return oneDimGrid.size();
}

/**
  Sort input vector x in increasing order, and then sort y accordingly.
  The underlining idea was taken from Accelerated C++ (textbook).
**/
/*
void Alglib_CubicSpline2D::sortXYvectors(vector<double>& x, vector<double>& y)
{
  ;
}
*/

void Alglib_CubicSpline2D::buildInterpolant(double xleft, double xright, double ybottom, double ytop)
{
  vector<double> xgrid, ygrid;
  buildGridPoints(xleft, xright, ybottom, xgrid);
  buildGridPoints(ybottom, ytop, xleft, ygrid);

  typedef vector<double>::size_type sz;
  sz sizex = xgrid.size();
  sz sizey = ygrid.size();
  double _x[sizex];
  double _y[sizey];
  double _f[sizex*sizey];

  for (sz i = 0; i < sizex; i++) _x[i] = xgrid[i];
  for (sz i = 0; i < sizey; i++) _y[i] = ygrid[i];
  for (sz i = 0; i < sizey; i++) {
    for (sz j = 0; j < sizex; j++) {
      _f[i*sizex+j] = function(_x[j], _y[i], parameter);
    }
  }

  xarray.setcontent(sizex, _x);
  yarray.setcontent(sizey, _y);
  farray.setcontent(sizex*sizey, _f);
  spline2dbuildbicubicv(xarray, sizex, yarray, sizey, farray, 1, spline);
}
