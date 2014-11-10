//#include "nfit.h"
#include "BareStructureFactor.h"
//#include "utable.h"
#include "fileTools.h"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001


// dqr and dqz are step size
void createStructureFactor(double Kc, double B, double D, double T, double Lr,
                           double Mz, double qrmin, double qrmax, double dqr,
                           double qzmin, double qzmax, double dqz)
{
  BareStructureFactor mc;
  mc.read_in_utable("../dat/utab_nfit12.15.dat");  
  mc.setModelParameter(Kc, "Kc");
  mc.setModelParameter(B, "B");
  mc.setModelParameter(D, "D");
  mc.setModelParameter(T, "T");
  mc.setModelParameter(Lr, "Lr");
  mc.setModelParameter(Mz, "Mz"); 
  mc.setModelParameter(0, "mosaic");
  mc.setModelParameter(0, "edisp");
  
  vector<double> qrvec, qzvec, sfvec, xvec, yvec;
  mc.getBareStrFct(qrmin, qrmax, dqr, qzmin, qzmax, dqz, qrvec, qzvec, sfvec);
	
  saveMatrix(qrvec, qzvec, sfvec, "columns.dat");
  saveMatrix(qrvec.size(), qzvec.size(), sfvec, "matrix.dat");  
  saveRowVector(qrvec, "qr.dat");
  saveRowVector(qzvec, "qz.dat");
  
  for (unsigned int i = 0; i < qrvec.size(); i++) xvec.push_back(i);
  saveDoubleColumns(xvec, qrvec, "qr.dat");
  for (unsigned int i = 0; i < qzvec.size(); i++) yvec.push_back(i);
  saveDoubleColumns(yvec, qzvec, "qz.dat");
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


void user_interface()
{
  double Kc = getInput("Enter Kc [default: 6e-13]: ", 6e-13);
  double B = getInput("Enter B [default: 0.5e13]: ", 0.5e13);
  double D = getInput("Enter D [default: 62.8]: ", 62.8);
  double T = getInput("Enter T [default: 30]: ", 30);
  double Lr = getInput("Enter Lr [default: 2500]: ", 2500);
  double Mz = getInput("Enter Mz [default: 10]: ", 10);
  double qrmin = getInput("Enter qr min [default: -0.05]: ", -0.05);
  double qrmax = getInput("Enter qr max [default: 0.05]: ", 0.05);
  double dqr = getInput("Enter delta qr [default: 0.001]: ", 0.001);
  double qzmin = getInput("Enter qz min (must be >= 0.001) [default: 0.05]: ", 0.05);
  double qzmax = getInput("Enter qz max [default: 0.45]: ", 0.45);
  double dqz = getInput("Enter delta qz [default: 0.001]: ", 0.001);
  
  createStructureFactor(Kc, B, D, T, Lr, Mz, 
                        qrmin, qrmax, dqr, 
                        qzmin, qzmax, dqz);
}


void automate()
{
  BareStructureFactor mc;
  mc.read_in_utable("../dat/utab_nfit12.15.dat");  
  mc.setModelParameter(6e-13, "Kc");
  mc.setModelParameter(0.5e13, "B");
  mc.setModelParameter(62.8, "D");
  mc.setModelParameter(30, "T");
  mc.setModelParameter(0, "mosaic");
  mc.setModelParameter(0, "edisp");
  
  vector<double> Lr_vec, Mz_vec;
  Lr_vec.push_back(250);
  Lr_vec.push_back(2500);
  Mz_vec.push_back(2);
  Mz_vec.push_back(6);
  Mz_vec.push_back(10);
  
  double qrmin = -0.05;
  double qrmax = 0.05;
  double dqr = 0.001;
  double qzmin = 0.05;
  double qzmax = 0.45;
  double dqz = 0.001;
  
  for (unsigned int i = 0; i < Lr_vec.size(); i++) {
    for (unsigned int j = 0; j < Mz_vec.size(); j++) {
      double Lr = Lr_vec[i];
      double Mz = Mz_vec[j];
      mc.setModelParameter(Lr, "Lr");
      mc.setModelParameter(Mz, "Mz"); 
      
      vector<double> qrvec, qzvec, sfvec, xvec, yvec;
      mc.getBareStrFct(qrmin, qrmax, dqr, qzmin, qzmax, dqz, qrvec, qzvec, sfvec);      
      
      string columns = "columns_" + to_string(int(Lr)) + "_" + to_string(int(Mz)) + ".dat";
      saveMatrix(qrvec, qzvec, sfvec, columns.c_str());
      
      string matrix = "matrix_" + to_string(int(Lr)) + "_" + to_string(int(Mz)) + ".dat";
      saveMatrix(qrvec.size(), qzvec.size(), sfvec, matrix.c_str());  

      string qr = "qr_" + to_string(int(Lr)) + "_" + to_string(int(Mz)) + ".dat";  
      for (unsigned int i = 0; i < qrvec.size(); i++) xvec.push_back(i);
      saveDoubleColumns(xvec, qrvec, qr.c_str());
      
      string qz = "qz_" + to_string(int(Lr)) + "_" + to_string(int(Mz)) + ".dat";          
      for (unsigned int i = 0; i < qzvec.size(); i++) yvec.push_back(i);
      saveDoubleColumns(yvec, qzvec, qz.c_str());    
    }
  }
}


int main()
{
  //user_interface();
  automate();
  
  return 0;
}
