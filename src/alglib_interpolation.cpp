#include <math.h>
#include <iostream>
#include <queue>
#include <vector>
#include <algorithm>
#include "alglib/src/stdafx.h"
#include "alglib/src/interpolation.h"

using namespace alglib;


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

  vector<double>::size_type sizex = xgrid.size();
  vector<double>::size_type sizey = ygrid.size();
  double _x[sizex];
  double _y[sizey];
  double _f[sizex*sizey];

  for (int i = 0; i < sizex; i++) _x[i] = xgrid[i];
  for (int i = 0; i < sizey; i++) _y[i] = ygrid[i];
  for (int i = 0; i < sizey; i++) {
    for (int j = 0; j < sizex; j++) {
      _f[i*sizex+j] = function(_x[j], _y[i], parameter);
    }
  }

  xarray.setcontent(sizex, _x);
  yarray.setcontent(sizey, _y);
  farray.setcontent(sizex*sizey, _f);
  spline2dbuildbicubicv(xarray, size_x, yarray, size_y, farray, 1, spline);
}

double Alglib_CubicSpline2D::evaluate(double vx, double vy)
{
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

  double initVleft = Function(initLeft, line, para);
  double initVright = Function(initRight, line, para);

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
    vmiddle = Function(middle, line, para);

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

    return oneDimGrid.size();
  }

  /*The points are not guaranteed to be in sorted order, so x values must be
    sorted. Finally corresponding function values are added to the vector of y
  values.*/
  sort(oneDimGrid.begin(), oneDimGrid.end());
}

/**
  Sort input vector x in increasing order, and then sort y accordingly.
  The underlining idea was taken from Accelerated C++ (textbook).
**/
void Alglib_CubicSpline2D::sortXYvectors(vector<double>& x, vector<double>& y)
{
  ;
}


void Alglib_CubicSpline2D::buildInterpolant(double xleft, double xright, double ybottom, double ytop)
{
  vector<double> xgrid, ygrid;
  buildGridPoints(xleft, xright, ybottom, xgrid);
  buildGridPoints(ybottom, ytop, xleft, ygrid);

  vector<double>::size_type sizex = xgrid.size();
  vector<double>::size_type sizey = ygrid.size();
  double _x[sizex];
  double _y[sizey];
  double _f[sizex*sizey];

  for (int i = 0; i < sizex; i++) _x[i] = xgrid[i];
  for (int i = 0; i < sizey; i++) _y[i] = ygrid[i];
  for (int i = 0; i < sizey; i++) {
    for (int j = 0; j < sizex; j++) {
      _f[i*sizex+j] = Function(_x[j], _y[i], para);
    }
  }

  xarray.setcontent(sizex, _x);
  yarray.setcontent(sizey, _y);
  farray.setcontent(sizex*sizey, _f);
  spline2dbuildbicubicv(xarray, size_x, yarray, size_y, farray, 1, spline);
}
