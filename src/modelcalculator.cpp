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

#include "globalVariables.h"
#include "modelcalculator.h"

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

ModelCalculator::ModelCalculator()
{
  //Turn off the error handler in gsl_integration_qag
  gsl_set_error_handler_off();
  //For qag numerical integration routine
  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  //Interpolation support functions
  spsumn.init(s_sumnWrapper, this);
  spHr.init(s_HrWrapper, this);
  algStrFct.init(s_StrFctWrapper, this);
  algRotated.init(s_rotatedWrapper, this);
  algMosaic.init(s_convolveMosaicWrapper, this);
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
  in = mc.in;
  out = mc.out;

  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  utable = mc.utable;

  spsumn.init(s_sumnWrapper, this);
  spHr.init(s_HrWrapper, this);
  algStrFct.init(s_StrFctWrapper, this);
  algRotated.init(s_rotatedWrapper, this);
  algMosaic.init(s_convolveMosaicWrapper, this);
  resetTOL();

	// call the updaters
	avgLrChanged();
	avgMzChanged();
	KcBDTChanged();
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
  spsumn.setol(g_spsumnAbserr, g_spsumnRelerr, g_spsumnMindx, g_spsumnMaxdx);
  spHr.setol(g_spHrAbserr, g_spHrRelerr, g_spHrMindx, g_spHrMaxdx);
}


void buildCCDStrFct(double _qxMax, double _qzMax);
{
  double qxMax = _qxMax + 20*beamSigma;
  qMax = getqMax(qxMax, _qzMax);
  
  
  buildInterpForMosaicStrFct(qMax);
  

}


void ModelCalculator::getqMax(double qxMax, double qzMax)
{
  double qy = max(fabs(getLowerLimit(qxMax, qzMax, wavelength)), 
                  fabs(getUpperLimit(qxMax, qzMax, wavelength)));
  double qrMax = sqrt(qxMax*qxMax + qy*qy);
  return max(qrMax, qzMax);
}


void ModelCalculator::buildInterpForMosaicStrFct()
{


}



double ModelCalculator::evalStrFct(double q, double theta)
{
  return algStrFct.evaluate(q*sin(theta), q*cos(theta));
}


/******************************************************************************
Build 2D interpolant for the structure factor. 

For every qz of 0.001, build an interpolant for the qr slice, evaluate it at 
every qr of 0.001, and store the value in the qrqzStrFct vector. Then, the qr, 
qz, and qrqzStrFct vectors are used to build 2D interpolant, which can 
evaluate the structure factor at any arbitrary set of qr and qz.

The structure factor must be calculated qr slice by slice because the
calculation of the theory is most efficient along a fixed qz value. Therefore,
wrapping a bunch of calculated qr slices by a 2D interpolation is an
effienet way to build the structure factor in qr, qz coordinates. 
******************************************************************************/
void ModelCalculator::buildInterpForStrFct(double qrMax, double qzMax)
{
  vector<double> qrvec;
  vector<double> qzvec;
  // Create vectors
  qrvec.push_back();
  qzvec.push_back();
  
  typedef vector<double>::size_type sz;
  for (sz i = 0; i < qzvec.size(); i++) {
    qrSlice(qrMax, qzvec[i]);
    for (sz j = 0; j < qrvec.size(); j++) {
      qrqzStrFct.push_back(exp(spStrFct.val(log(fabs(qrvec[j])+SMALLNUM))));
    }
  }
  
  algStrFct.buildInterpolant(qrvec, qzvec, qrqzStrFct);  
}


/******************************************************************************
Initiate building the interpolants for the pure structure factor, mosaic-
convolved structure factor, and rotated structure factor.
******************************************************************************/
void ModelCalculator::qrSlice(double qrMax, double qz)
{
	if (qrMax < 0 || qz < 0)
		throw domain_error("ERROR: qrMax and qz must both take positive values");
  buildInterpForRotatedStrFct(0, qxhigh + 20*beamSigma);
}


/******************************************************************************
In order to build rotated structure factor interpolant, the program needs
mosaic-convolved structure factor. This function decides the necessary
range of qr for mosaic-convolved structure factor. The upper limit on
the range is slightly larger than qx high limit for rotated structure factor
because qr >= qx for any value of qy.
******************************************************************************/
void ModelCalculator::buildInterpForRotatedStrFct(double qxlow, double qxhigh)
{
  double qy = max(fabs(getLowerLimit(qxhigh, qz, wavelength)), 
                  fabs(getUpperLimit(qxhigh, qz, wavelength)));
  double qrhigh = sqrt(qxhigh*qxhigh + qy*qy);
  buildInterpForMosaicStrFct(qxlow, qrhigh);
  spRotated.findPoints( log(qxlow+SMALLNUM), log(qxhigh+SMALLNUM) );
}


void ModelCalculator::buildInterpForMosaicStrFct(double qrlow, double qrhigh)
{
  buildInterpForStrFct(0, qrhigh+0.101);
  qrUpperLimit = qrhigh + 0.1;
  spMosaic.findPoints( log(qrlow+SMALLNUM), log(qrhigh+SMALLNUM) );
}


void ModelCalculator::buildInterpForStrFct(double qrlow, double qrhigh)
{
  spStrFct.findPoints( log(qrlow+SMALLNUM), log(qrhigh+SMALLNUM) );
}


/******************************************************************************
FuncLmdif::modelCalcThread calls this function.
This function returns the interpolated value of the rotated structure factor
at the input qx. The interpolant returns the log of a rotated structure factor, 
so this function takes the exponential of the interpolated value.
******************************************************************************/
double ModelCalculator::getCCDStrFct(double qx)
{
  // 0.05 hard coded for now. For some test cases, this worked well.
	if (bFWHM > 0.05) {
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
	currQx = qx;
  double result, abserr;
  gsl_function F;
  F.function = &s_convIntegrand;
  F.params = this;
  // set integration limits
  double lowerLimit = qx - 10*beamSigma;
  double upperLimit = qx + 10*beamSigma;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
Return the integrand for the integration that performs beam convolution
******************************************************************************/
double ModelCalculator::s_convIntegrand(double qx, void *params)
{
	ModelCalculator *p = (ModelCalculator *)params;
	return exp( p->spRotated.val(log(fabs(qx) + SMALLNUM)) )
	       * p->beamProfile(qx - p->currQx);
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
	return log(p->rotated(tmpQx));
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
	currQx = qx;
	double result, abserr;
	gsl_function F;
	F.function = &s_qyIntegrandWrapper;
	F.params = this;
	double lowerLimit = getLowerLimit(qx, qz, wavelength);
	double upperLimit = getUpperLimit(qx, qz, wavelength);
	gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
	                    WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
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
double ModelCalculator::s_qyIntegrandWrapper(double qy, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
	double qr = sqrt(p->currQx*p->currQx + qy*qy);
	return p->getMosaicStrFct(qr);
}


/******************************************************************************
Return the structure factor after mosaic spread is treated.
If mosaic > 0.001, it will return an interpolated value of the structure 
factor with mosaic convolution.
If mosaic is less than 0.001, it will return an interpolated value of the 
structure factor without any additional operation done on itself.
******************************************************************************/
double ModelCalculator::getMosaicStrFct(double qr)
{
  // Takes the absolute value of qr to avoid a very small negative value 
  // representing zero
  const double negligible_mosaic = 0.001 * PI / 180;
  if (mosaic > negligible_mosaic) {
    return exp(spMosaic.val(log(fabs(qr)+SMALLNUM)));
  } else {
    return exp(spStrFct.val(log(fabs(qr)+SMALLNUM)));
  }
}


/****************************************************************************** 
Wrapper needed for FunSupport object, spMosaic
******************************************************************************/
double ModelCalculator::s_convolveMosaicWrapper(double logqr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double tmpQr = fabs(exp(logqr) - SMALLNUM);
  return log(p->convolveMosaic(tmpQr));
}


/******************************************************************************
This function performs convolution of the structure factor with the mosaic
distribution.
******************************************************************************/
double ModelCalculator::convolveMosaic(double q, double theta)
{
  curr_q = q;
  curr_theta = theta;
  double result, abserr;
  gsl_function F;
  F.function = &s_mosaicIntegrandWrapper;
  F.params = this;
  double lowerLimit = 0;
  double upperLimit = PI / 2;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
This function returns the integrand of the mosaic convolution. 
******************************************************************************/
double ModelCalculator::s_mosaicIntegrandWrapper(double theta, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
	return p->interpStrFct(p->curr_q, curr_theta-theta) *
	       p->mosaicDist(theta);
}


/******************************************************************************
Mosaic spread distribution. Mosaic angle is approximated by qr/qz, which is
good for small angle.
******************************************************************************/
double ModelCalculator::mosaicDist(double theta)
{
  return 2 * mosaic / PI / (4*theta*theta + mosaic*mosaic);
}


void get_q_theta(double qr, double qz, double& q, double& theta)
{
  q = sqrt(qr*qr + qz*qz);
  theta = atan(qr/qz);
}


void get_qr_qz(double q, double theta, double& qr, double& qz)
{
  qr = q * sin(theta);
  qz = q * cos(theta);
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
    cout << "\nNegative structure factor was obtained at" << endl;
    cout << "qr: " << p->currQr << " qz: " << p->qz << " Value: " << ret << endl;
    cout << p->Kc << " " << p->B << " " << p->avgLr << " " << p->avgMz << " "
         << p->dspacing << " " << p->T << " " << p->wavelength << endl;
  }
  return log(ret);
}


/******************************************************************************
This function calculates the structure factor as a function of qr for a fixed 
value of qz specified by qz variable
******************************************************************************/
double ModelCalculator::StrFct(double qr)
{
  currQr = qr;
  double result, abserr;
  gsl_function F;
  F.function = &s_Fr;
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
double ModelCalculator::s_Fr(double r, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
  return r * p->spHr.val(r) * bessel_J0(p->currQr*r) *
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
  return Pr(Lr, p->avgLr) * Lr * Lr * finiteSizeFactor(p->curr_r/Lr);
}


/******************************************************************************
In-plane domain size distribution
Can be changed to another function easily
******************************************************************************/
double Pr(double Lr, double avgLr)
{
  return exp(-Lr / avgLr) / avgLr;
}


/******************************************************************************
Finite size factor in r direction, F_r(r/Lr)
******************************************************************************/
double finiteSizeFactor(double x)
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
double getLowerLimit(double qx, double qz, double wavelength)
{
	double theta = 0.5 * asin(sqrt(qx*qx+qz*qz)*wavelength/2/PI);
	return -sqrt((4*PI*sin(theta)/wavelength)*(4*PI*sin(theta)/wavelength) -
	             qx*qx -qz*qz);
}


/******************************************************************************
Return the upper limit of the qy integration for the currenct value of qx and 
qz. The equation is derived in an appropriate report.
******************************************************************************/
double getUpperLimit(double qx, double qz, double wavelength)
{
	double omega = 2 * asin(wavelength*qz/4/PI);
	return getLowerLimit(qx, qz, wavelength)*cos(omega)*
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
		case Var_Kc:                 Kc = a;  KcBDTChanged(); break;
		case Var_B:                   B = a;  KcBDTChanged(); break;
		case Var_D:            dspacing = a;  KcBDTChanged(); break;
		case Var_T:                   T = a;  KcBDTChanged(); break;
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
  KcBDTChanged();
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


void ModelCalculator::getCCDStrFct(double qxlow, double qxhigh, double qz,
                                   vector<double>& qxv, vector<double>& sfv)
{
  setSliceParameter(qz);
  QxSlice(qxlow, qxhigh);
  
  for (double qx = qxlow; qx < qxhigh; ) {
    qxv.push_back(qx);
    sfv.push_back(getCCDStrFct(qx));
    qx += 0.001;
  }
}


void ModelCalculator::getRotatedStrFct(double qxlow, double qxhigh, double _qz, 
                                       vector<double>& qxvec, vector<double>& rsf)
{
	setSliceParameter(_qz);
	buildInterpForRotatedStrFct(qxlow, qxhigh);

	for (double qx = qxlow; qx < qxhigh; ) {
		qxvec.push_back(qx);
		rsf.push_back( exp( spRotated.val(log(qx+SMALLNUM)) ) );
		qx += 0.001;
	}
}


void ModelCalculator::getMosaicStrFct(double qrlow, double qrhigh, double _qz,
                                      vector<double>& qr, vector<double>& sf)
{
	if (qrlow < 0 || qrhigh < 0)
		throw domain_error("qrlow and qrhigh must both take positive values");
	setSliceParameter(_qz);
  buildInterpForMosaicStrFct(qrlow, qrhigh);

	for (double q = qrlow; q < qrhigh;) {
		qr.push_back(q);
		sf.push_back(getMosaicStrFct(q));
		q += 0.001;
	}
}


void ModelCalculator::getStrFct(double qrlow, double qrhigh, double _qz, 
                                vector<double>& qrvec, vector<double>& sf)
{
	if (qrlow < 0 || qrhigh < 0)
		throw domain_error("qrlow and qrhigh must both take positive values");
	setSliceParameter(_qz);
	buildInterpForStrFct(qrlow, qrhigh);

	for (double qr = qrlow; qr < qrhigh;) {
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


void saveDoubleColumns(vector<double>& xvec, vector<double>& yvec, 
                       const char *filename)
{
  cout << "xvec.size(): " << xvec.size() << " yvec.size(): " << yvec.size() 
       << endl;
  if (xvec.size() != yvec.size()) {
    cerr << "\nThe size of input vectors must be identical for " 
         << "saveDoubleColumns function to work\n" << endl;
    return;
  }
  ofstream myfile;
  myfile.open(filename);
  typedef vector<double>::size_type vec_sz;
  for (vec_sz i = 0; i < xvec.size(); i++) {
    myfile << xvec[i] << " " << yvec[i] << endl;
  }
}
