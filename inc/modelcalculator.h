#ifndef GUARD_MODELCALCULATOR_H
#define GUARD_MODELCALCULATOR_H

#include <math.h>
#include <gsl/gsl_integration.h>
//#include <tcl.h>
#include <vector>
#include <string>
#include "Para.h"
#include "tvds.h"
#include "nfit.h"
#include "funsupport.h"
#include "utable.h"
#include "alglib_interpolation.h"
//#include "alglib/src/stdafx.h"
//#include "alglib/src/interpolation.h"

#define PI 3.1415926535897932384626433832795


double bessel_J0(double x);
enum Var stringToVarEnum(const std::string &);

/******************************************************************************
ModelCalculator class: 
This is the class that calculates the theoretical structure factor.
FuncLmdif should bind to this class object for chi square minimization
******************************************************************************/
class ModelCalculator {
public:
  ModelCalculator();
  //This constructor provides a means of cloning a ModelCalculator object
  ModelCalculator(const ModelCalculator& mc);
  ~ModelCalculator(){cleanup();}
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

	// functions to build the interpolating function for the structure factor
	void qxSlice(double, double);
	void buildInterpForRotatedStrFct(double, double);
	void init(double, double, double, double);
	void buildInterpForMosaicStrFct(double, double, double, double);
	void buildInterpForStrFct(double, double, double, double);
	void qrSlice(double, double);
	void setSliceParameter(double qz);
	double getCCDStrFct(double); 

	// core functions that calculate the structure factor
	double beamConvolutedStrFct(double);
	static double s_beamConvIntegrand(double, void *);
	double beamProfile(double);
	static double s_rotatedWrapper(double, void *);
	double rotated(double);
	static double s_rotatedIntegrand(double, void *);
	static double s_convolveMosaicWrapper(double, double, void *);
	double convolveMosaic(double, double);
	static double s_convolveMosaicIntegrand(double, void *);
	double mosaicDist(double);
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
  double getLowerLimit(double, double);
  double getUpperLimit(double, double);
  
  // a bunch of functions that are useful in testing and debugging
  void get_spStrFctPoints(std::vector<double>&, std::vector<double>&);
  void get_spsumnPoints(std::vector<double>&, std::vector<double>&);
  void get_spHrPoints(std::vector<double>&, std::vector<double>&);
  void get_spRotatedPoints(std::vector<double>&, std::vector<double>&);
  void get2DStrFct(double, double, double, std::vector<double>&, std::vector<double>&);
  void getStrFct(double, double, double, std::vector<double>&, std::vector<double>&);
  void getSumn(double, double, double, std::vector<double>&, std::vector<double>&);
  void getHr(double, double, std::vector<double>&, std::vector<double>&);
  void getRotatedStrFct(double, double, double, std::vector<double>&, std::vector<double>&);
  void getMosaicStrFct(double, double, double, std::vector<double>&, std::vector<double>&);
  void getMosaicStrFct(double, double, double, double, std::vector<double>&, 
                       std::vector<double>&, std::vector<double>&);
  void getCCDStrFct(double, double, double, std::vector<double>&, std::vector<double>&);
  void read_struct_factor(const char *);
  double evalStrFct(double, double);
  void getBareStrFct(double, double, double, double, double, double, 
                     std::vector<double>&, std::vector<double>&, std::vector<double>&);
private:
  double Kc, B, dspacing, T;
  double avgLr, avgMz;
  double edisp;
  double mosaic;
  double wavelength, pixelSize, bFWHM, sdistance, beamSigma;

  // cutoff for r integration and sum over n
  double cutoff_r, cutoff_n;
  // qz value for the qx slice
  double qz;
  // gsl workspace
  gsl_integration_workspace *workspace;
  // interpolation support object for sum over n
  FunSupport spsumn;
  // interpolation support objects for the structure factor
  FunSupport spStrFct;
  Alglib_CubicSpline2D algStrFct;
  // interpolation support object for effective finite size factor in r
  FunSupport spHr;
	// interpolation support object for rotated structure factor
  FunSupport spRotated;
  // interpolation support object for mosaic-convolved structure factor
  Alglib_CubicSpline2D algMosaic;
  
  Utable utable;
  double curr_r, currtheta, currqx, currq, currqr;
  
  // should be called after parameters are changed. considered as updater
  void avgLrChanged();
  void avgMzChanged();
  void set_eta_lambda();
  void set_beamSigma();
  // mosaic is in radian, but the program takes the input in degrees
  void set_mosaic(double a) {mosaic = PI * a / 180;}
};

#endif
