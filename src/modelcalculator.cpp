/* 
 * CCD structure factor is in (qx, qz) space, map calculated slice by slice
 * rotated structure factor is in (qx, qz) space, map calculated slice by slice
 * mosaic-convolved structure factor is in (q, theta) space, full 2D map
 * pure structure factor is in (qr, qz) space, full 2D map
 */

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
//#include "alglib/src/stdafx.h"
//#include "alglib/src/interpolation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <errno.h>
#include <stdexcept>
#include <algorithm>
#include <time.h>

#include "globalVariables.h"
#include "modelcalculator.h"
#include "fileTools.h"


using std::string; using std::ifstream; using std::getline; 
using std::istringstream;
using std::vector; using std::cout; using std::endl;
using std::cerr;

//using namespace alglib;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001
#define WORKSPACE_SIZE 10000
#define KEY 6
#define ARB_SCALE 100000
const double NEGLIGIBLE_MOSAIC = 0.001 * PI / 180; // in radian 
const double NEGLIGIBLE_BEAM_SIZE = 0.1; // in pixel
const double QMIN = 0.02; // in invrse-Angstrom
const double QSTEP = 0.001;
const double THETAMAX = PI / 2.5; 
const double THETASTEP = 0.1 * PI / 180; // in radian
const double QRSTEP = 0.001;
const double QZSTEP = 0.001;
const double QZMIN = 0.001;


ModelCalculator::ModelCalculator()
{
  //Turn off the error handler in gsl_integration_qag
  gsl_set_error_handler_off();
  //For qag numerical integration routine
  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  //Interpolation support functions
  spsumn.init(s_sumnWrapper, this);
  spHr.init(s_HrWrapper, this);
  spStrFct.init(s_StrFctWrapper, this);
  spRotated.init(s_rotatedWrapper, this);
  algMosaic.initialize(s_convolveMosaicWrapper, this);
  resetTOL();
}


/******************************************************************************
Cloning constructor that initializes necessary structures and copies
variables to provide a duplicate of a ModelCalculator object already
in use.
******************************************************************************/
ModelCalculator::ModelCalculator(const ModelCalculator& mc)
{
  // Copy necessary variables
  Kc = mc.Kc;
  B = mc.B;
  dspacing = mc.dspacing;
  T = mc.T;
  avgLr = mc.avgLr;
  avgMz = mc.avgMz;
  edisp = mc.edisp;
  mosaic = mc.mosaic;
  wavelength = mc.wavelength;
  pixelSize = mc.pixelSize;
	bFWHM = mc.bFWHM;
  sdistance = mc.sdistance;
	beamSigma = mc.beamSigma;
  cutoff_r = mc.cutoff_r;
  cutoff_n = mc.cutoff_n;

  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  utable = mc.utable;

  spsumn.init(s_sumnWrapper, this);
  spHr.init(s_HrWrapper, this);
  spStrFct.init(s_StrFctWrapper, this);
  spRotated.init(s_rotatedWrapper, this);
  algMosaic.initialize(s_convolveMosaicWrapper, this);
  resetTOL();

	// call the updaters
	avgLrChanged();
	avgMzChanged();
	set_beamSigma();
}


/******************************************************************************
Clean up resouces, should be called before the object is destructed
******************************************************************************/
void ModelCalculator::cleanup()
{
  gsl_integration_workspace_free(workspace);
}


/*******************************************************************************
set tolerance for interpolation
*******************************************************************************/
void ModelCalculator::resetTOL()
{
  spStrFct.setol(g_spStrFctAbserr, g_spStrFctRelerr, g_spStrFctMindx, g_spStrFctMaxdx);
  spsumn.setol(g_spsumnAbserr, g_spsumnRelerr, g_spsumnMindx, g_spsumnMaxdx);
  spHr.setol(g_spHrAbserr, g_spHrRelerr, g_spHrMindx, g_spHrMaxdx);
  spRotated.setol(g_spRotatedAbserr, g_spRotatedRelerr, g_spRotatedMindx, g_spRotatedMaxdx);
}

/******************************************************************************
This method builds a qx-slice (along constant qz) for qx = 0 to qxMax
using the (mosaic-convolve) 2D structure factor map built by init function.
Specifically, it builds an interpolant that takes care of so-called CCD 
integration, which accounts for the sample rotation during x-ray exposure.

Call this method before calling getCCDStrFct method.

Call init(qxMax, qzMax) before caling this method.
******************************************************************************/
void ModelCalculator::qxSlice(double qxMax, double _qz)
{
	if (qxMax < 0 || _qz < 0)
		throw domain_error("ERROR: qxMax and qz must both take positive values");
  buildInterpForRotatedStrFct(qxMax, _qz);
}


void ModelCalculator::buildInterpForRotatedStrFct(double qxMax, double _qz)
{
  qz = _qz;
  //cout << "qz: " << qz << endl;
  qxMax = qxMax + 20*beamSigma;  
  spRotated.findPoints( log(0+SMALLNUM), log(qxMax+SMALLNUM) );
}


/******************************************************************************
Initiates building of mosaic-convolved structure factor.

The inputs, qxMin, qxMax, qzMin, and qzMax define the size of a desired CCD 
structure factor map.

This method adds 20*beamSigma to qxMax to account for an extra
convolution due to the horizontal beam wdith. Then, it computes qrMax, which is
larger than qxMax because of the qy integration due to the sample rotation.

The mosaic-convolved structure factor is a map in (qr, qz) space, but
the integration is done along q. This map will extend from qzmin to qzmax 
and from qrmin to qrmax. qrmin is assumed to be zero. 

It is necessary to call this function before calling qxSlice.

If mosaic is negligible, the program will bypass building of the mosaic-
convolved structure factor and only build a 2D map of the pure structre 
factor.
******************************************************************************/
void ModelCalculator::init(double qxMin, double qxMax, 
                           double qzMin, double qzMax)
{
  double qrMin = 0;
  qxMax = qxMax + 20*beamSigma;
  double qyMax = max(fabs(getLowerLimit(qxMax, qzMax)), 
                     fabs(getUpperLimit(qxMax, qzMax)));
  double qrMax = sqrt(qxMax*qxMax + qyMax*qyMax);
  
  if (mosaic > NEGLIGIBLE_MOSAIC) { 
    buildInterpForMosaicStrFct(qrMin, qrMax, qzMin, qzMax);
  } else {
    buildInterpForStrFct(qrMin, qrMax, qzMin, qzMax);
    cout << "finished StrFct: " << clock()*1e-6 << " " << qrMax << " " 
         << qzMax << endl;
  }
}


/******************************************************************************
Build mosaic-convolved structure factor.

This function builds a 2D interpolant for the mosaic-convolved structure 
factor via algMosaic, which is an Alglib_CubicSpline2D object.
******************************************************************************/
void ModelCalculator::buildInterpForMosaicStrFct(double qrMin, double qrMax, 
                                                 double qzMin, double qzMax)
{
  // The structure factor map must span in both qr and qz up to qMax
  double qMax = sqrt(qrMax*qrMax+qzMax*qzMax);
  buildInterpForStrFct(0, qMax+0.001, 0.001, qMax+0.001);
  cout << "finished StrFct: " << clock()*1e-6 << " " << qMax << " " 
       << qMax << endl;
  
  vector<double> qrvec;
  vector<double> qzvec;
  // qrMax+QRSTEP ensures that qrMax is included (avoiding floating point issue)
  for (double qr = qrMin; qr < qrMax+QRSTEP;) {
    qrvec.push_back(qr);
    qr += QRSTEP;
  }
  for (double qz = qzMin; qz < qzMax+QZSTEP; ) {
    qzvec.push_back(qz);
    qz += QZSTEP;
  }
  algMosaic.buildInterpolant(qrvec, qzvec);
  cout << "finished algMosaic: " << clock()*1e-6 << endl;
}


/******************************************************************************
Build 2D interpolant for the structure factor. 

For every qz of 0.001, build an interpolant for the qr slice, evaluate it at 
every qr of 0.001, and store the value in the qrqzStrFct vector. Then, the qr, 
qz, and qrqzStrFct vectors are used to build 2D interpolant, which can 
evaluate the structure factor at any arbitrary set of qr and qz.

The structure factor must be calculated qr-slice by slice because the
calculation of the theory is most efficient along a fixed qz value. Therefore,
wrapping a bunch of calculated qr slices by a 2D interpolation is an
effienet way to build the structure factor in qr, qz coordinates. 
******************************************************************************/
void ModelCalculator::buildInterpForStrFct(double qrMin, double qrMax, 
                                           double qzMin, double qzMax)
{
  vector<double> qrvec;
  vector<double> qzvec;
  vector<double> qrqzStrFctvec;

  for (double qr = qrMin; qr < qrMax+QRSTEP;) {
    qrvec.push_back(qr);
    qr += QRSTEP;
  }
  for (double qz = qzMin; qz < qzMax+QZSTEP;) {
    qzvec.push_back(qz);
    qz += QZSTEP;
  }
  
  typedef vector<double>::size_type sz;
  for (sz i = 0; i < qzvec.size(); i++) {
    qrSlice(qrvec.back(), qzvec[i]);
    for (sz j = 0; j < qrvec.size(); j++) {
      qrqzStrFctvec.push_back(exp(spStrFct.val(log(fabs(qrvec[j])+SMALLNUM))));
    }
  }
  
  algStrFct.buildInterpolant(qrvec, qzvec, qrqzStrFctvec);  
}


/******************************************************************************
Initiate building the interpolants for the pure structure factor.
******************************************************************************/
void ModelCalculator::qrSlice(double qrMax, double _qz)
{
	if (qrMax < 0 || _qz < 0)
		throw domain_error("ERROR: qrMax and qz must both take positive values");
	setSliceParameter(_qz);
  spStrFct.findPoints( log(0+SMALLNUM), log(qrMax+SMALLNUM) );
}


/******************************************************************************
This method returns the theoretical value of the CCD structure factor
at the input qx value. The qz value is determined by the input to the 
qxSlice method.

FuncLmdif::funclmdif calls this method and take the difference between 
the returned value and the data point.

Call the qxSlice method before calling this method.
******************************************************************************/
double ModelCalculator::getCCDStrFct(double qx)
{
	if (bFWHM > NEGLIGIBLE_BEAM_SIZE) {
		return beamConvolutedStrFct(qx);
	} else {
		return exp( spRotated.val(log(qx + SMALLNUM)) );
	}
}


/******************************************************************************
This function takes care of beam FWHM convolution with rotated structure
factor.
******************************************************************************/
double ModelCalculator::beamConvolutedStrFct(double qx)
{
	currqx = qx;
  double result, abserr;
  gsl_function F;
  F.function = &s_beamConvIntegrand;
  F.params = this;
  // qxMax + 20*beamSigma is the largest qx value it can go. To be safe,
  // the integration only goes upto qxMax + 10*beamSigma
  double lowerLimit = qx - 10*beamSigma;
  double upperLimit = qx + 10*beamSigma;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel_low,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
Return the integrand for the integration that performs beam convolution
******************************************************************************/
double ModelCalculator::s_beamConvIntegrand(double qx, void *params)
{
	ModelCalculator *p = (ModelCalculator *)params;
	double tmp = exp( p->spRotated.val(log(fabs(qx) + SMALLNUM)) );
	//cout << qx << " " << p->currqx << " " << tmp << endl;
	return exp( p->spRotated.val(log(fabs(qx) + SMALLNUM)) )
	       * p->beamProfile(qx - p->currqx);
}


/******************************************************************************
Return the beam profile in q-space. It is assumed to be Gaussian.
******************************************************************************/
double ModelCalculator::beamProfile(double qx)
{
	return exp(-qx * qx / 2 / beamSigma / beamSigma) / sqrt(2 * PI) / beamSigma;
}


/******************************************************************************
Static function that returns the logarithmic value of the rotated structure 
factor at the input log(qx). (Note that it is base e, not base 10)
It uses the absolute value of qx. This is to deal with qx = 0 point, where
floating point arithmetic could potentially yield a very slightly negative 
value that represents zero in doulbe precision floating point.
******************************************************************************/
double ModelCalculator::s_rotatedWrapper(double logqx, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
	double tmpQx = fabs(exp(logqx) - SMALLNUM);
	double tmp = log(p->rotated(tmpQx));
	//cout << logqx << " " << tmpQx << " " << tmp << endl;
	return tmp;
}


/******************************************************************************
This function performs the integration over the structure factor along qy 
direction at the input qx value
The lower and upper integration limits depend on qx, so their values get 
retrieved by ModelCalculator::getLowerLimit() and 
ModelCalculator::getUpperLimit() functions
*******************************************************************************/
double ModelCalculator::rotated(double qx)
{
  currqx = qx;
	double result, abserr;
	gsl_function F;
	F.function = &s_rotatedIntegrand;
	F.params = this;
	double lowerLimit = getLowerLimit(qx, qz);
	double upperLimit = getUpperLimit(qx, qz);
	gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel_low,
	                    WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
	if (result < 0) cout << "rotated method returning a negative value" << endl;
	return result;
}


/******************************************************************************
This function returns the integrand for the qy integration, which is
the interpolated value of the mosaic convolved structure factor at given qx, 
qy, and qz.
The interpolation over the structure factor is done in (qr,qz) tuple.
Calculates qr corresponding to the input qy and currQx.
This wrapper is necessary for GSL qag.
******************************************************************************/
double ModelCalculator::s_rotatedIntegrand(double qy, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
	double qr = sqrt(p->currqx*p->currqx + qy*qy);
	double tmp;
	if (p->mosaic > NEGLIGIBLE_MOSAIC) {
	  if ((tmp = p->algMosaic.evaluate(qr, p->qz)) == NAN) cout << "NAN";
	  return p->algMosaic.evaluate(qr, p->qz);
	} else {
	  return p->algStrFct.evaluate(qr, p->qz);
	}
}


/****************************************************************************** 
Wrapper needed for AlglibCublicSpline2D object, algMosaic
******************************************************************************/
double ModelCalculator::s_convolveMosaicWrapper(double qr, double qz, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double q = sqrt(qr*qr + qz*qz);
  double theta = atan(qr/qz);
  double ret = p->convolveMosaic(q, theta);
  if (ret < 0) {
    cout << "negative mosaic structure factor obtained at "
         << q << " " << theta << endl;
  }
  return ret;
}


/******************************************************************************
This function performs convolution of the structure factor with the mosaic
distribution.

The upper limit is pi/3 instead of pi/2 b/c qz = 0 is pathological
in the structure factor.
******************************************************************************/
double ModelCalculator::convolveMosaic(double q, double theta)
{
  currq = q;
  currtheta = theta;
  double result, abserr;
  gsl_function F;
  F.function = &s_convolveMosaicIntegrand;
  F.params = this;
  //double lowerLimit = -PI/3;
  //double upperLimit = PI/3;
  double lowerLimit = theta - PI/2.1;
  double upperLimit = PI/2.1;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel_low,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
This function returns the integrand of the mosaic convolution. 
******************************************************************************/
double ModelCalculator::s_convolveMosaicIntegrand(double theta, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double qr = p->currq * sin(p->currtheta-theta);
  double qz = p->currq * cos(p->currtheta-theta);
  double tmp;
  try {
    tmp = p->algStrFct.evaluate(fabs(qr), fabs(qz));
  } catch (domain_error e) {
    cout << e.what() << endl 
         << "Set the return value of algStrFct.evaluate to zero" << endl;
    tmp = 0;
  }
	return tmp * p->mosaicDist(theta);
}


/******************************************************************************

******************************************************************************/
double ModelCalculator::evalStrFct(double qr, double qz)
{
  //return algStrFct.evaluate(fabs(qr), fabs(qz));
  return algMosaic.evaluate(qr, qz);
}


/******************************************************************************
Mosaic spread distribution. 
******************************************************************************/
double ModelCalculator::mosaicDist(double theta)
{
  return 2 * mosaic / PI / (4*theta*theta + mosaic*mosaic);
}


/******************************************************************************
Static wrapper for structure factor along qr direction.
The input is log(qr), base e.
Returns log(S(qr,qz)).
This wrapper is necessary for FunSupport object.
******************************************************************************/
double ModelCalculator::s_StrFctWrapper(double logqr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double tmpQr = fabs(exp(logqr) - SMALLNUM);
  double ret = p->StrFct(tmpQr);
  if (ret < 0) {
    cout << "\nNegative structure factor was obtained at" << endl
         << "qr: " << p->currqr << " qz: " << p->qz << " Value: " << ret << endl;
    cout << p->Kc << " " << p->B << " " << p->avgLr << " " << p->avgMz << " "
         << p->dspacing << " " << p->T << " " << p->wavelength << endl;
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
double ModelCalculator::StrFct(double qr)
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
double ModelCalculator::s_StrFctIntegrand(double r, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
  return r * p->spHr.val(r) * bessel_J0(p->currqr*r) *
         p->spsumn.val(log(r+SMALLNUM)) / ARB_SCALE;
}


/******************************************************************************
static function that wraps ModelCalculator::sumn, so that it can be passed to
FunSupport::init
FunSupport evaluates this function at constant log(r) step. This function
converts log(r) to r. SMALLNUM is used to avoid Nan coming from exp(log(0)).
When compiled under g++, as of 7/3/2013, log(0) yields -inf, but exp(log(0)) 
yeilds Nan, instead of 0.
SMALLNUM is added to r, like log(r+SMALLNUM), when an interpolated value gets
retrieved.
******************************************************************************/
double ModelCalculator::s_sumnWrapper(double logr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double tmpR = fabs(exp(logr) - SMALLNUM);
  return p->sumn(tmpR);
}


/******************************************************************************
Calculate the sum over layer index, n.
It is equal to f_1(r)
******************************************************************************/
double ModelCalculator::sumn(double r)
{
  double sum = 0;
  double un, t1, t2;
  double qz2 = qz * qz;
  double qzsig2 = (qz*edisp) * (qz*edisp);
  double lambda = 1e16 * sqrt(Kc/B) / dspacing;
  double eta = 2.17 * (T+273.15) / sqrt(Kc*B) / dspacing /dspacing;

  //Calculate n=0 term separately because it has 1/2 weight relatively to n>0 terms
  int n = 0;
  un = utable.getHeightDiffFunc(n, r, lambda, eta, dspacing);
  t1 = 1 + qzsig2*un;
  t2 = n * dspacing;
  sum += Hz(n) * sqrt(1/t1) * exp( -0.5*(qz2*un+t2*t2*qzsig2)/t1 ) * cos(t2*qz/t1);
  for(n = 1; n < cutoff_n; n++) {
    un = utable.getHeightDiffFunc(n, r, lambda, eta, dspacing);
    t1 = 1 + qzsig2*un;
    t2 = n * dspacing;
    sum += 2 * Hz(n) * sqrt(1/t1) * exp( -0.5*(qz2*un+t2*t2*qzsig2)/t1 ) * cos(t2*qz/t1);
  }
  return sum;
}


/******************************************************************************
Static function that wrapps ModelCalculator::Hr
******************************************************************************/
double ModelCalculator::s_HrWrapper(double r, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  return p->Hr(r);
}


/******************************************************************************
Effective finite size factor in r direction, H_r(r)
The upper limit is chosen to be 10 times the cutoff in r integration.
The factor of 10 is somewhat arbitrary.
******************************************************************************/
double ModelCalculator::Hr(double r)
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
double ModelCalculator::s_HrIntegrand(double Lr, void *params)
{
  ModelCalculator *p = (ModelCalculator *)params;
  return p->Pr(Lr) * Lr * Lr * p->finiteSizeFactor(p->curr_r/Lr);
}


/******************************************************************************
In-plane domain size distribution
Can be changed to another function easily
******************************************************************************/
double ModelCalculator::Pr(double Lr)
{
  return exp(-Lr / avgLr) / avgLr;
}


/******************************************************************************
Finite size factor in r direction, F_r(r/Lr)
******************************************************************************/
double ModelCalculator::finiteSizeFactor(double x)
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
double ModelCalculator::Hz(int n)
{
  return exp(-n/avgMz);
}


/******************************************************************************
Return the lower limit of the qy integration for the current value of qx and 
qz. The equation is derived in an appropirate report.
******************************************************************************/
double ModelCalculator::getLowerLimit(double qx, double qz)
{
	double theta = 0.5 * asin(sqrt(qx*qx+qz*qz)*wavelength/2/PI);
	return -sqrt((4*PI*sin(theta)/wavelength)*(4*PI*sin(theta)/wavelength) -
	             qx*qx -qz*qz);
}


/******************************************************************************
Return the upper limit of the qy integration for the currenct value of qx and 
qz. The equation is derived in an appropriate report.
******************************************************************************/
double ModelCalculator::getUpperLimit(double qx, double qz)
{
	double omega = 2 * asin(wavelength*qz/4/PI);
	return getLowerLimit(qx, qz)*cos(omega)*
	       (1+tan(omega)*tan(omega)) + qz*tan(omega);
}


/******************************************************************************
given a slice with qz = tqz, update necessary temporary variables
and tables.
******************************************************************************/
void ModelCalculator::setSliceParameter(double _qz)
{
  qz = _qz;
  spsumn.findPoints( log(0+SMALLNUM), log(cutoff_r+SMALLNUM) );
}


/******************************************************************************
update domain size parameters
When they are modified, the cutoffs and the interpolant must also be modified
******************************************************************************/
void ModelCalculator::avgLrChanged()
{
  cutoff_r = 50 * avgLr;
  spHr.findPoints(0, cutoff_r);
}


void ModelCalculator::avgMzChanged()
{
	cutoff_n = 20 * avgMz;
}


/******************************************************************************
Converts beam FWHM in pixels to beam Gaussian sigma in q-space units
******************************************************************************/
void ModelCalculator::set_beamSigma()
{
	double deltaQPerPixel = 2 * PI * pixelSize / wavelength / sdistance;
	// bFWHM / 2.3548 is equal to beam Gaussian sigma
	beamSigma = deltaQPerPixel * bFWHM / 2.3548;
}


/******************************************************************************
This function sets a model parameter to an input value.
A model parameter is specified by std::string name.
See stringToVarEnum function for a list of strings that can be
passed as an input name.

value: an input value
name: an input name specifying which parameter to be set to an input value
******************************************************************************/
void ModelCalculator::setModelParameter(double value, const char *_name)
{
	string name(_name);
	setModelParameter(value, name);
}


void ModelCalculator::setModelParameter(double value, const string& name)
{
  setpara(value, stringToVarEnum(name));
}


/******************************************************************************
set parameters for Sccd calculation

a : the value to be set
idx: the index for the corresponding parameter
******************************************************************************/
void ModelCalculator::setpara(double a, int idx)
{
  switch(idx) {
		case Var_Kc:                 Kc = a;  break;
		case Var_B:                   B = a;  break;
		case Var_D:            dspacing = a;  break;
		case Var_T:                   T = a;  break;
		case Var_avgLr:           avgLr = a;  avgLrChanged(); break;
		case Var_avgMz:           avgMz = a;  avgMzChanged(); break;
		case Var_mosaic:                       set_mosaic(a); break;
		case Var_edisp:           edisp = a;                  break;
		case Var_wavelength: wavelength = a; set_beamSigma(); break;
		case Var_pixelSize:   pixelSize = a; set_beamSigma(); break;
		case Var_bFWHM:           bFWHM = a; set_beamSigma(); break;
		case Var_sdistance:   sdistance = a; set_beamSigma(); break;
		case Var_beamSigma:   beamSigma = a;                  break;
		default: ;
  }
}


/******************************************************************************
make the ModelCalculator in sync with parameter set p
******************************************************************************/
void ModelCalculator::paraSet(Para *p)
{
  Kc = p->Kc;
  B = p->B;
  dspacing = p->D;
  T = p->T;
  avgLr = p->avgLr;
  avgMz = p->avgMz;
  edisp = p->edisp;
  wavelength = p->setup.wavelength;
  pixelSize = p->setup.pz;
  bFWHM = p->bFWHM;
  sdistance = p->setup.s;
  avgLrChanged();
  avgMzChanged();
  set_beamSigma();
  set_mosaic(p->mosaic);
}


void ModelCalculator::get_spStrFctPoints(vector<double>& x, vector<double>& y)
{
	spStrFct.getPoints(x, y);
}


void ModelCalculator::get_spsumnPoints(vector<double>& x, vector<double>& y)
{
	spsumn.getPoints(x, y);
}


void ModelCalculator::get_spHrPoints(vector<double>& x, vector<double>& y)
{
	spHr.getPoints(x, y);
}


void ModelCalculator::get_spRotatedPoints(vector<double>& x,vector<double>& y)
{
	spRotated.getPoints(x, y);
}


void ModelCalculator::getCCDStrFct(double qxlow, double qxhigh, double _qz,
                                   vector<double>& qxv, vector<double>& sfv)
{
  init(0, qxhigh, QZMIN, _qz);
  qxSlice(qxhigh, _qz);
  
  for (double qx = qxlow; qx < qxhigh; ) {
    qxv.push_back(qx);
    sfv.push_back(getCCDStrFct(qx));
    qx += 0.001;
  }
}


void ModelCalculator::getRotatedStrFct(double qxlow, double qxhigh, double _qz, 
                                       vector<double>& qxvec, vector<double>& rsf)
{
	init(0, qxhigh, QZMIN, _qz);
	qxSlice(qxhigh, _qz);

	for (double qx = qxlow; qx < qxhigh; ) {
		qxvec.push_back(qx);
		rsf.push_back( exp( spRotated.val(log(qx+SMALLNUM)) ) );
		qx += 0.001;
	}
}


void ModelCalculator::getMosaicStrFct(double qrlow, double qrhigh, double _qz,
                                      vector<double>& qrv, vector<double>& sf)
{
	init(qrlow, qrhigh, _qz, _qz);
	for (double qr = qrlow; qr < qrhigh;) {
		qrv.push_back(qr);
		sf.push_back(algMosaic.evaluate(qr, _qz));
		qr += 0.001;
	}
}


void ModelCalculator::getMosaicStrFct(double qrlow, double qrhigh, 
                                      double qzlow, double qzhigh,
                                      vector<double>& qrv, vector<double>& qzv,
                                      vector<double>& sfv)
{
  for (double qr = qrlow; qr < qrhigh; ) {
    qrv.push_back(qr);
    qr += 0.001;
  }
  for (double qz = qzlow; qz < qzhigh; ) {
    qzv.push_back(qz);
    qz += 0.001;
  }
  init(0, qrv.back(), QZMIN, qzv.back());
  typedef vector<double>::size_type sz;
  for (sz i = 0; i < qzv.size(); i++) {
    for (sz j = 0; j < qrv.size(); j++) {
      try {
        if (mosaic > NEGLIGIBLE_MOSAIC) {
          sfv.push_back(algMosaic.evaluate(sqrt(qrv[j]*qrv[j]+qzv[i]*qzv[i]), 
                                         atan(qrv[j]/qzv[i])));
        } else {
          sfv.push_back(algStrFct.evaluate(qrv[j], qzv[i]));
        }
      } catch (exception& e) {
        cout << e.what() << endl;
        if (mosaic > NEGLIGIBLE_MOSAIC) {
          cout << "vx: " << sqrt(qrv[j]*qrv[j]+qzv[i]*qzv[i]) << " " 
               << "vy: " << atan(qrv[j]/qzv[i]) << endl;
        } else {
          cout << "vx: " << qrv[j] << " " << "vy: " << qzv[i] << endl;
        }
      }
    }
  }
}
                                      

// Calculate bare structure factor, without mosaic spread or resolution function
// convolution. Return S(qr,qz) in qrv, qzv, and sfv vectors.
// sfv.size() = qrv.size() * qzv.size()
void ModelCalculator::getBareStrFct(double qrmin, double qrmax, double deltaqr,
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
      sfv.push_back(algStrFct.evaluate(fabs(qrv[j]), qzv[i]));
    }
  }
}


// function name is misleading
void ModelCalculator::get2DStrFct(double qrlow, double qrhigh , double _qz,
                                  vector<double>& qrv, vector<double>& sfv)
{
  init(0, qrhigh, QZMIN, _qz);
  
  for (double qr = qrlow; qr < qrhigh; ) {
    qrv.push_back(qr);
    sfv.push_back(algStrFct.evaluate(qr, _qz));
    qr += 0.001;
  }
}
                                  

void ModelCalculator::getStrFct(double qrlow, double qrhigh, double _qz, 
                                vector<double>& qrvec, vector<double>& sf)
{
	qrSlice(qrhigh, _qz);
	
	for (double qr = qrlow; qr < qrhigh+0.001;) {
		qrvec.push_back(qr);
		sf.push_back( exp( spStrFct.val(log(qr + SMALLNUM)) ) );
		qr += 0.001;
	}
}


void ModelCalculator::getSumn(double rlow, double rhigh, double _qz, 
                              vector<double>& rvec, vector<double>& zvec)
{
	setSliceParameter(_qz);
	for (double r = rlow; r < rhigh; ) {
		rvec.push_back(r);
		zvec.push_back( exp( spsumn.val(log(r+SMALLNUM)) ) );
		r += 10;
	}
}


void ModelCalculator::getHr(double rlow, double rhigh, 
                            vector<double>& rvec, vector<double>& hvec)
{
	avgLrChanged();
	for (double r = rlow; r < rhigh; ) {
		rvec.push_back(r);
		hvec.push_back(spHr.val(r));
		r += 100;
	}
}


void ModelCalculator::read_struct_factor(const char *filename)
{
  vector<double> qrvec, qzvec, qrqzStrFctvec;
  ifstream myfile;
  myfile.open(filename);
  string line;
  double tmp;
  
  getline(myfile, line);
  istringstream iss(line);
  while (iss >> tmp) {
    qrvec.push_back(tmp);
  }
  
  getline(myfile, line);
  istringstream iss2(line);
  while (iss2 >> tmp) {
    qzvec.push_back(tmp);
  }
  
  getline(myfile, line);
  istringstream iss3(line);
  while (iss3 >> tmp) {
    qrqzStrFctvec.push_back(tmp);
    //if (tmp < 1e3) cout << tmp << endl;
  }

  algStrFct.buildInterpolant(qrvec, qzvec, qrqzStrFctvec);
}
