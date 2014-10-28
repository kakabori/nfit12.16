//#include "nfit.h"
#include "BareStructureFactor.h"
#include "MosaicStructureFactor.h"
//#include "utable.h"
#include "fileTools.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001


// Create mosaic convolved structure factor
// dqr and dqz are step size
void createStructureFactor(double Kc, double B, double D, double T, double Lr,
                           double Mz, double mosaic, 
                           double qrmin, double qrmax, double dqr,
                           double qzmin, double qzmax, double dqz)
{
  MosaicStructureFactor mc;
  mc.read_in_utable("../dat/utab_nfit12.15.dat");  
  mc.setModelParameter(Kc, "Kc");
  mc.setModelParameter(B, "B");
  mc.setModelParameter(D, "D");
  mc.setModelParameter(T, "T");
  mc.setModelParameter(Lr, "Lr");
  mc.setModelParameter(Mz, "Mz"); 
  mc.setModelParameter(mosaic, "mosaic");
  mc.setModelParameter(0, "edisp");
  
  vector<double> qrvec, qzvec, sfvec;
  mc.getMosaicSF(qrmin, qrmax, dqr, qzmin, qzmax, dqz, qrvec, qzvec, sfvec);
	
  saveMatrix(qrvec, qzvec, sfvec, "columns.dat");
  saveMatrix(qrvec.size(), qzvec.size(), sfvec, "matrix.dat");  
  saveRowVector(qrvec, "qr.dat");
  saveRowVector(qzvec, "qz.dat");
}

// s is the string that gets displayed
// x is the default value
double getInput(string s, double x)
{
  string input;
  cout << s;
  getline(cin, input);
  if (!input.empty()) {
    istringstream stream(input);
    stream >> x;
  } 
  return x;
}


int main()
{
  double Kc = getInput("Enter Kc [default: 6e-13]: ", 6e-13);
  double B = getInput("Enter B [default: 2e13: ", 2e13);
  double D = getInput("Enter D [default: 62.8]: ", 62.8);
  double T = getInput("Enter T [default: 30]: ", 30);
  double Lr = getInput("Enter Lr [default: 2500]: ", 2500);
  double Mz = getInput("Enter Mz [default: 10]: ", 10);
  double mosaic = getInput("Enter mosaic [default: 0.1]: ", 0.1);
  double qrmin = getInput("Enter qr min [default: -0.1]: ", -0.1);
  double qrmax = getInput("Enter qr max [default: 0.1]: ", 0.1);
  double dqr = getInput("Enter delta qr [default: 0.001]: ", 0.001);
  double qzmin = getInput("Enter qz min (must be >= 0.001) [default: 0.05]: ", 0.05);
  double qzmax = getInput("Enter qz max [default: 0.25]: ", 0.25);
  double dqz = getInput("Enter delta qz [default: 0.001]: ", 0.001);
  
  createStructureFactor(Kc, B, D, T, Lr, Mz, mosaic,
                        qrmin, qrmax, dqr, 
                        qzmin, qzmax, dqz);
  return 0;
}
