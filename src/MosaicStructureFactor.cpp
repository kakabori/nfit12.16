#include <sstream>
#include <errno.h>
#include <stdexcept>
#include <algorithm>
#include <time.h>
#include "cubature.h"
#include "MosaicStructureFactor.h"

using std::string; using std::vector; using std::cout; using std::endl;
using std::cerr;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001
const double QMIN = 0.02; // in invrse-Angstrom
const double QSTEP = 0.001;
const double THETAMAX = PI / 2.5; 
const double THETASTEP = 0.1 * PI / 180; // in radian
const double QRSTEP = 0.001;
const double QZSTEP = 0.001;
const double QZMIN = 0.001;


// Initiate creation of mosaic convolved structure factor
void MosaicStructureFactor::init(double qrMin, double qrMax,
                                 double qzMin, double qzMax)
{
  double qMax = sqrt(qrMax*qrMax + qzMax*qzMax);
  // For now, hard code those limits
  double qrlow = 0;
  double qrhigh = qMax + 0.001;
  double qzlow = 0;
  double qzhigh = qMax + 0.001;
  double qrDelta = 0.001;
  double qzDelta = 0.001;
  cout << "1: " << clock()*1e-6 << endl;
  BareStructureFactor::init(qrlow, qrhigh, qzlow, qzhigh);
  cout << "2: " << clock()*1e-6 << endl;
  
  vector<double> qrvec, qzvec, sfvec;
  for (double qr = qrMin; qr < qrMax + 0.5*qrDelta; ) {
    qrvec.push_back(qr);
    qr += qrDelta;
  }
  for (double qz = qzMin; qz < qzMax + 0.5*qzDelta; ) {
    qzvec.push_back(qz);
    qz += qzDelta;
  }
  typedef vector<double>::size_type vec_sz;
  for (vec_sz i = 0; i < qzvec.size(); i++) {
    for (vec_sz j = 0; j < qrvec.size(); j++) {
      sfvec.push_back(calcMosaicConvSF(qrvec[j], qzvec[i]));
    }
  }
  cout << "3: " << clock()*1e-6 << endl;
  interp2dMosaicSF.buildInterpolant(qrvec, qzvec, sfvec);
  cout << "4: " << clock()*1e-6 << endl;
}


// Return mosaic convolved structure factor
double MosaicStructureFactor::evalMosaicSF(double qr, double qz)
{
  return interp2dMosaicSF.evaluate(fabs(qr), qz);
}


// Calculate mosaic convolved structure factor at a point (qr, qz)
double MosaicStructureFactor::calcMosaicConvSF(double qr, double qz)
{
  // Calculate corresponding (q, theta) from (qr, qz)
  q = sqrt(qr*qr + qz*qz);
  theta = atan(qr/qz);
  
  // Set up cubature
  double xmin[2] = {0, 0}, xmax[2] = {PI/6, 2*PI}, val, err;
  hcubature(1, s_f, this, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
  //hcubature_v(1, s_f2, this, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
    
  return val;
}


double MosaicStructureFactor::mosaicDist(double alpha)
{
  return 2 * mosaic / PI / (4*alpha*alpha + mosaic*mosaic);
}


int MosaicStructureFactor::s_f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
  MosaicStructureFactor *msf = (MosaicStructureFactor *)fdata;
  double q = msf->q;
  double theta = msf->theta;
  double alpha = x[0];
  double beta = x[1];
  double theta_prime = getThetaPrime(theta, alpha, beta);
  double qr = q * sin(theta_prime);
  double qz = q * cos(theta_prime);
  double tmp;
  try {
    tmp = msf->evalStrFct(fabs(qr), fabs(qz)) * msf->mosaicDist(alpha);
  } catch (exception& e) {
    cout << e.what() << endl;
    tmp = 0;
  }
  fval[0] = tmp; 
  return 0;
}


double getThetaPrime(double theta, double alpha, double beta)
{
  return acos(sin(theta)*sin(alpha)*cos(beta)+cos(theta)*cos(alpha));
}


int MosaicStructureFactor::s_f2(unsigned ndim, unsigned npts, const double *x, void *fdata, unsigned fdim, double *fval) 
{
  MosaicStructureFactor *msf = (MosaicStructureFactor *)fdata;
  double q = msf->q;
  double theta = msf->theta;
  for (unsigned j = 0; j < npts; ++j) {
    double alpha = x[j*ndim + 0];
    double beta = x[j*ndim + 1];
    double theta_prime = getThetaPrime(theta, alpha, beta);
    double qr = q * sin(theta_prime);
    double qz = q * cos(theta_prime);
    fval[j] = msf->evalStrFct(fabs(qr), fabs(qz)) * msf->mosaicDist(alpha);
  }
  return 0; // success
}


void MosaicStructureFactor::setpara(double a, int idx)
{
  BareStructureFactor::setpara(a, idx);
  switch(idx) {
    case Var_mosaic: set_mosaic(a); break;
    default: ;
  }
}


// Calculate bare structure factor, without mosaic spread or resolution function
// convolution. Return S(qr,qz) in qrv, qzv, and sfv vectors.
// sfv.size() = qrv.size() * qzv.size()
void MosaicStructureFactor::getMosaicSF(
  double qrmin, double qrmax, double deltaqr,
  double qzmin, double qzmax, double deltaqz,
  vector<double>& qrv, vector<double>& qzv, vector<double>& sfv)
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
      sfv.push_back(interp2dMosaicSF.evaluate(abs(qrv[j]), qzv[i]));
    }
  }
}

