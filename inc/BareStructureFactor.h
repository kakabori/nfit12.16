#ifndef GUARD_BARESTRUCTUREFACTOR_H
#define GUARD_BARESTRUCTUREFACTOR_H

#include <gsl/gsl_integration.h>
#include <math.h>
#include <vector>
#include <string>
#include "utable.h"
#include "funsupport.h"
#include "alglib_interpolation.h"
#include "Para.h"
#include "tvds.h"
#include "nfit.h"

#define PI 3.1415926535897932384626433832795

//double bessel_J0(double x);
enum Var stringToVarEnum(const std::string &);

class BareStructureFactor {
public:
  BareStructureFactor();
  //This constructor provides a means of cloning
  BareStructureFactor(const BareStructureFactor&);
  ~BareStructureFactor() {cleanup();}
  void cleanup();
  
  // setters and getters
  void setModelParameter(double, const std::string&);
  void setModelParameter(double, const char *);
  void setpara(double, int);
  void paraSet(Para *p);

  // Utable utilities
  void read_in_utable(const char *f){utable.readUtableFile(f);}
  
  // FunSupport utilities
  void resetTOL();
  
  // Functions to build the structuer factor interpolant
  void init(double, double, double, double);
  void buildInterpForStrFct(double, double, double, double);
  void qrSlice(double, double, double);
  void setSliceParameter(double);
  
  // Evaluate interpolated structure factor
  double evalStrFct(double, double); 
  void getBareStrFct(double, double, double, double, double, double, 
                     std::vector<double>&, std::vector<double>&, std::vector<double>&);
  
  // Functions to calculate the structure factor
  static double s_StrFctWrapper(double, void *);
  double StrFct(double);
  static double s_StrFctIntegrand(double, void *);
  static double s_sumnWrapper(double, void *);
  double sumn(double);
  static double s_HrWrapper(double, void *);
  double Hr(double);
  static double s_HrIntegrand(double, void *);
  double Pr(double);
  double finiteSizeFactor(double);
  double Hz(int);

private:
  double Kc, B, D, T, avgLr, avgMz, edisp;
  Utable utable;  
  
  // should be called after parameters are changed. considered as updater
  void avgLrChanged();
  void avgMzChanged();
  void set_eta_lambda();
  
  // cutoff for r integration and sum over n
  double cutoff_r, cutoff_n, currqr, curr_r;
  double qz;
  
  // gsl workspace
  gsl_integration_workspace *workspace;
  // interpolation support object for sum over n
  FunSupport interp1dsumn;
  // interpolation support objects for the structure factor
  FunSupport interp1dStrFct;
  Alglib_CubicSpline2D interp2dStrFct;
  // interpolation support object for effective finite size factor in r
  FunSupport interp1dHr;
};

#endif
