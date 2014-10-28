/* 
 * CCD structure factor is in (qx, qz) space, map calculated slice by slice
 * rotated structure factor is in (qx, qz) space, map calculated slice by slice
 * mosaic-convolved structure factor is in (q, theta) space, full 2D map
 * pure structure factor is in (qr, qz) space, full 2D map
 */

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
//#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <errno.h>
#include <stdexcept>
#include <algorithm>
#include <time.h>
#include <interp2d.h>
#include <interp2d_spline.h>

#include "globalVariables.h"
#include "BareStructureFactor.h"
#include "fileTools.h"

using std::string; using std::ifstream; using std::getline; 
using std::istringstream;
using std::vector; using std::cout; using std::endl;
using std::cerr;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001
#define WORKSPACE_SIZE 10000
#define KEY 6
#define ARB_SCALE 100000
const double QMIN = 0.02; // in invrse-Angstrom
const double QSTEP = 0.001;
const double THETAMAX = PI / 2.5; 
const double THETASTEP = 0.1 * PI / 180; // in radian
const double QRSTEP = 0.001;
const double QZSTEP = 0.001;
const double QZMIN = 0.001;


BareStructureFactor::BareStructureFactor()
{
  //Turn off the error handler in gsl_integration_qag
  gsl_set_error_handler_off();
  //For qag numerical integration routine
  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  //Interpolation support functions
  interp1dsumn.init(s_sumnWrapper, this);
  interp1dHr.init(s_HrWrapper, this);
  interp1dStrFct.init(s_StrFctWrapper, this);
  resetTOL();
}


/******************************************************************************
Cloning constructor that initializes necessary structures and copies
variables to provide a duplicate of a BareStructureFactor object already
in use.
******************************************************************************/
BareStructureFactor::BareStructureFactor(const BareStructureFactor& mc)
{
  // Copy necessary variables
  Kc = mc.Kc;
  B = mc.B;
  D = mc.D;
  T = mc.T;
  avgLr = mc.avgLr;
  avgMz = mc.avgMz;
  cutoff_r = mc.cutoff_r;
  cutoff_n = mc.cutoff_n;
  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
  utable = mc.utable;
  interp1dsumn.init(s_sumnWrapper, this);
  interp1dHr.init(s_HrWrapper, this);
  interp1dStrFct.init(s_StrFctWrapper, this);
  resetTOL();

  // call the updaters
  avgLrChanged();
  avgMzChanged();
}


/******************************************************************************
Clean up resouces, should be called before the object is destructed
******************************************************************************/
void BareStructureFactor::cleanup()
{
  gsl_integration_workspace_free(workspace);
}


/*******************************************************************************
set tolerance for interpolation
*******************************************************************************/
void BareStructureFactor::resetTOL()
{
  interp1dStrFct.setol(g_interp1dStrFctAbserr, g_interp1dStrFctRelerr, g_interp1dStrFctMindx, g_interp1dStrFctMaxdx);
  interp1dsumn.setol(g_interp1dsumnAbserr, g_interp1dsumnRelerr, g_interp1dsumnMindx, g_interp1dsumnMaxdx);
  interp1dHr.setol(g_interp1dHrAbserr, g_interp1dHrRelerr, g_interp1dHrMindx, g_interp1dHrMaxdx);
}


/******************************************************************************

******************************************************************************/
void BareStructureFactor::init(double qrMin, double qrMax, 
                               double qzMin, double qzMax)
{
  if (qrMin < 0) {
    cout << "Input qrMin < 0. Setting qrMin = 0." << endl;
    qrMin = 0;
  }
  if (qzMin < QZMIN) {
    cout << "Input qzMin is less than the allowed value. Setting qzMin = " 
         << QZMIN << "." << endl;
    qzMin = QZMIN;
  }
  buildInterpForStrFct(qrMin, qrMax, qzMin, qzMax);
}


/******************************************************************************

******************************************************************************/
void BareStructureFactor::buildInterpForStrFct(double qrMin, double qrMax, 
                                           double qzMin, double qzMax)
{
  vector<double> qrvec;
  vector<double> qzvec;
  vector<double> qrqzStrFctvec;

  for (double qr = qrMin; qr < qrMax+0.5*QRSTEP;) {
    qrvec.push_back(qr);
    qr += QRSTEP;
  }
  for (double qz = qzMin; qz < qzMax+0.5*QZSTEP;) {
    qzvec.push_back(qz);
    qz += QZSTEP;
  }
  
  typedef vector<double>::size_type sz;
  for (sz i = 0; i < qzvec.size(); i++) {
    qrSlice(qrvec.front(), qrvec.back(), qzvec[i]);
    for (sz j = 0; j < qrvec.size(); j++) {
      qrqzStrFctvec.push_back(exp(interp1dStrFct.val(log(fabs(qrvec[j])+SMALLNUM))));
    }
  }
  
  interp2dStrFct.buildInterpolant(qrvec, qzvec, qrqzStrFctvec);  
}


/******************************************************************************
Build interpolant for a qr slice for qrMin <= qr <= qrMax at qz
******************************************************************************/
void BareStructureFactor::qrSlice(double qrMin, double qrMax, double qz)
{
  if (qrMin < 0) {
    throw domain_error("ERROR: qrMin < 0 not allowed");   
  } else if (qrMax < qrMin) {
    throw domain_error("ERROR: qrMax < qrMin not allowed");
  } else if (qz < QZMIN) {
    throw domain_error("ERROR: qz < QZMIN not allowed");
  }
  setSliceParameter(qz);
  interp1dStrFct.findPoints( log(qrMin+SMALLNUM), log(qrMax+SMALLNUM) );
}


/******************************************************************************
Return structure factor at (qr, qz) using the 2D interpolant
******************************************************************************/
double BareStructureFactor::evalStrFct(double qr, double qz)
{
  return interp2dStrFct.evaluate(fabs(qr), qz);
}


/******************************************************************************
Static wrapper for structure factor along qr direction.
The input is log(qr), base e.
Returns log(S(qr,qz)).
This wrapper is necessary for FunSupport object.
******************************************************************************/
double BareStructureFactor::s_StrFctWrapper(double logqr, void *ptr)
{
  BareStructureFactor *p = (BareStructureFactor *)ptr;
  double tmpQr = fabs(exp(logqr) - SMALLNUM);
  double ret = p->StrFct(tmpQr);
  if (ret < 0) {
    cout << "\nNegative structure factor was obtained at" << endl
         << "qr: " << p->currqr << " qz: " << p->qz << " Value: " << ret << endl;
    cout << p->Kc << " " << p->B << " " << p->avgLr << " " << p->avgMz << " "
         << p->D << " " << p->T << " " << endl;
    cout << "Recalculate the structure factor at +/- 0.0005 of the current qr\n"
         << "value and take the mean.\n" << endl;
    double tmp1 = p->StrFct(tmpQr-0.0005);
    double tmp2 = p->StrFct(tmpQr+0.0005);
    ret = (tmp1 + tmp2) / 2;
  }
  return log(ret);
}


/******************************************************************************
This function calculates the structure factor as a function of qr for a fixed 
value of qz specified by qz variable
******************************************************************************/
double BareStructureFactor::StrFct(double qr)
{
  currqr = qr;
  double result, abserr;
  gsl_function F;
  F.function = &s_StrFctIntegrand;
  F.params = this;
  double lowerLimit = 0;
  double upperLimit = cutoff_r;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
This calculates f_2(r,qr) integrand
It is equal to r * H_r(r) * J_0(qr*r) * f_1(r)
Two functions are their corresponding interpolation functions,
H_r(r) and f_1(r) are both one dimensional interpolation using FunSupport class
J_0(qr*r) is the zeroth order Bessel function of the first kind
******************************************************************************/
double BareStructureFactor::s_StrFctIntegrand(double r, void *ptr)
{
	BareStructureFactor *p = (BareStructureFactor *)ptr;
  return r * p->interp1dHr.val(r) * bessel_J0(p->currqr*r) *
         p->interp1dsumn.val(log(r+SMALLNUM)) / ARB_SCALE;
}


/******************************************************************************
static function that wraps BareStructureFactor::sumn, so that it can be passed to
FunSupport::init
FunSupport evaluates this function at constant log(r) step. This function
converts log(r) to r. SMALLNUM is used to avoid Nan coming from exp(log(0)).
When compiled under g++, as of 7/3/2013, log(0) yields -inf, but exp(log(0)) 
yeilds Nan, instead of 0.
SMALLNUM is added to r, like log(r+SMALLNUM), when an interpolated value gets
retrieved.
******************************************************************************/
double BareStructureFactor::s_sumnWrapper(double logr, void *ptr)
{
  BareStructureFactor *p = (BareStructureFactor *)ptr;
  double tmpR = fabs(exp(logr) - SMALLNUM);
  return p->sumn(tmpR);
}


/******************************************************************************
Calculate the sum over layer index, n.
It is equal to f_1(r)
******************************************************************************/
double BareStructureFactor::sumn(double r)
{
  double sum = 0;
  double un, t1, t2;
  double qz2 = qz * qz;
  double qzsig2 = (qz*edisp) * (qz*edisp);
  double lambda = 1e16 * sqrt(Kc/B) / D;
  double eta = 2.17 * (T+273.15) / sqrt(Kc*B) / D /D;

  //Calculate n=0 term separately because it has 1/2 weight relatively to n>0 terms
  int n = 0;
  un = utable.getHeightDiffFunc(n, r, lambda, eta, D);
  t1 = 1 + qzsig2*un;
  t2 = n * D;
  sum += Hz(n) * sqrt(1/t1) * exp( -0.5*(qz2*un+t2*t2*qzsig2)/t1 ) * cos(t2*qz/t1);
  for(n = 1; n < cutoff_n; n++) {
    un = utable.getHeightDiffFunc(n, r, lambda, eta, D);
    t1 = 1 + qzsig2*un;
    t2 = n * D;
    sum += 2 * Hz(n) * sqrt(1/t1) * exp( -0.5*(qz2*un+t2*t2*qzsig2)/t1 ) * cos(t2*qz/t1);
  }
  return sum;
}


/******************************************************************************
Static function that wrapps BareStructureFactor::Hr
******************************************************************************/
double BareStructureFactor::s_HrWrapper(double r, void *ptr)
{
  BareStructureFactor *p = (BareStructureFactor *)ptr;
  return p->Hr(r);
}


/******************************************************************************
Effective finite size factor in r direction, H_r(r)
The upper limit is chosen to be 10 times the cutoff in r integration.
The factor of 10 is somewhat arbitrary.
******************************************************************************/
double BareStructureFactor::Hr(double r)
{
  curr_r = r;
  double result, abserr;
  gsl_function F;
  F.function = &s_HrIntegrand;
  F.params = this;
  double lowerLimit = r;
  double upperLimit = 10 * cutoff_r;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return PI * result;
}


/******************************************************************************
The integrand in H_r(r) integration
******************************************************************************/
double BareStructureFactor::s_HrIntegrand(double Lr, void *params)
{
  BareStructureFactor *p = (BareStructureFactor *)params;
  return p->Pr(Lr) * Lr * Lr * p->finiteSizeFactor(p->curr_r/Lr);
}


/******************************************************************************
In-plane domain size distribution
Can be changed to another function easily
******************************************************************************/
double BareStructureFactor::Pr(double Lr)
{
  return exp(-Lr / avgLr) / avgLr;
}


/******************************************************************************
Finite size factor in r direction, F_r(r/Lr)
******************************************************************************/
double BareStructureFactor::finiteSizeFactor(double x)
{
  if (x > 1) {
    return 0;
  } else {
    return acos(x) - x*sqrt(1-x*x);
  }
}


/******************************************************************************
Effective finite size factor in z direction, H_z(nD)
******************************************************************************/
double BareStructureFactor::Hz(int n)
{
  return exp(-n/avgMz);
}


/******************************************************************************
given a slice with qz = tqz, update necessary temporary variables
and tables.
******************************************************************************/
void BareStructureFactor::setSliceParameter(double _qz)
{
  qz = _qz;
  interp1dsumn.findPoints( log(0+SMALLNUM), log(cutoff_r+SMALLNUM) );
}


/******************************************************************************
update domain size parameters
When they are modified, the cutoffs and the interpolant must also be modified
******************************************************************************/
void BareStructureFactor::avgLrChanged()
{
  cutoff_r = 50 * avgLr;
  interp1dHr.findPoints(0, cutoff_r);
}


void BareStructureFactor::avgMzChanged()
{
	cutoff_n = 20 * avgMz;
}


/******************************************************************************
This function sets a model parameter to an input value.
A model parameter is specified by std::string name.
See stringToVarEnum function for a list of strings that can be
passed as an input name.

value: an input value
name: an input name specifying which parameter to be set to an input value
******************************************************************************/
void BareStructureFactor::setModelParameter(double value, const char *_name)
{
	string name(_name);
	setModelParameter(value, name);
}


void BareStructureFactor::setModelParameter(double value, const string& name)
{
  setpara(value, stringToVarEnum(name));
}


/******************************************************************************
set parameters for Sccd calculation

a : the value to be set
idx: the index for the corresponding parameter
******************************************************************************/
void BareStructureFactor::setpara(double a, int idx)
{
  switch(idx) {
		case Var_Kc:                 Kc = a;  break;
		case Var_B:                   B = a;  break;
		case Var_D:            D = a;  break;
		case Var_T:                   T = a;  break;
		case Var_avgLr:           avgLr = a;  avgLrChanged(); break;
		case Var_avgMz:           avgMz = a;  avgMzChanged(); break;
		case Var_edisp:           edisp = a;                  break;
		default: ;
  }
}


/******************************************************************************
make the BareStructureFactor in sync with parameter set p
******************************************************************************/
void BareStructureFactor::paraSet(Para *p)
{
  Kc = p->Kc;
  B = p->B;
  D = p->D;
  T = p->T;
  avgLr = p->avgLr;
  avgMz = p->avgMz;
  edisp = p->edisp;
  avgLrChanged();
  avgMzChanged();
}
                                   

// Calculate bare structure factor, without mosaic spread or resolution function
// convolution. Return S(qr,qz) in qrv, qzv, and sfv vectors.
// sfv.size() = qrv.size() * qzv.size()
void BareStructureFactor::getBareStrFct(double qrmin, double qrmax, double deltaqr,
                                    double qzmin, double qzmax, double deltaqz,
                                    vector<double>& qrv,
                                    vector<double>& qzv,
                                    vector<double>& sfv)
{
  if (qzmin < QZMIN) {
    qzmin = QZMIN;
  }
  init(qrmin, qrmax+deltaqr, qzmin, qzmax);
  
  // adding 0.5 make sure to include qrmax
  for (double qr = qrmin; qr < qrmax + 0.5*deltaqr; ) {
    qrv.push_back(qr);
    qr += deltaqr;
  }
  for (double qz = qzmin; qz < qzmax + 0.5*deltaqz; ) {
    qzv.push_back(qz);
    qz += deltaqz;
  }  
  
  typedef vector<double>::size_type vec_sz;
  for (vec_sz i = 0; i < qzv.size(); i++) {
    for (vec_sz j = 0; j < qrv.size(); j++) {
      sfv.push_back(interp2dStrFct.evaluate(abs(qrv[j]), qzv[i]));
    }
  }
}

