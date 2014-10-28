#ifndef GUARD_MOSAICSTRUCTUREFACTOR_H
#define GUARD_MOSAICSTRUCTUREFACTOR_H

//#include <gsl/gsl_integration.h>
#include <math.h>
#include <vector>
#include <string>
//#include "utable.h"
//#include "funsupport.h"
#include "alglib_interpolation.h"
#include "Para.h"
#include "tvds.h"
#include "nfit.h"
#include "BareStructureFactor.h"

#define PI 3.1415926535897932384626433832795

class MosaicStructureFactor : public BareStructureFactor {
public:
  void setpara(double, int);
  void init(double, double, double, double);
  double evalMosaicSF(double, double);
  double calcMosaicConvSF(double, double);
  double mosaicDist(double);
  static int s_f(unsigned, const double *, void *, unsigned, double *);
  static int s_f2(unsigned, unsigned, const double *, void *, unsigned, double *);
  void getMosaicSF(double, double, double, double, double, double, vector<double>&, vector<double>&, vector<double>&);
private:
  void set_mosaic(double a) {mosaic = PI * a / 180;}
  double q, theta;
  Alglib_CubicSpline2D interp2dMosaicSF;
  double mosaic;
};

double getThetaPrime(double, double, double);

#endif
