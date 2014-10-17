#include "nfit.h"
#include "modelcalculator.h"
#include "utable.h"
#include "fileTools.h"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001

void createStructureFactor(double Kc, double B, double D, double T, double Lr,
                           double Mz, double edisp, double moasic, double wavelength,
                           double pixelSize, double bFWHM, double s)
{
  ModelCalculator mc;
  mc.read_in_utable("../dat/utab_nfit12.15.dat");
  
  mc.setModelParameter(Kc, "Kc");
  mc.setModelParameter(B, "B");
  mc.setModelParameter(D, "D");
  mc.setModelParameter(T, "T");
  mc.setModelParameter(Lr, "Lr");
  mc.setModelParameter(Mz, "Mz");
  mc.setModelParameter(edisp, "edisp");
  mc.setModelParameter(mosaic, "mosaic");
  mc.setModelParameter(wavelength, "wavelength");
  mc.setModelParameter(pixelSize, "pixelSize");
  mc.setModelParameter(bFWHM, "bFWHM");
  mc.setModelParameter(s, "s");  
}

int main()
{
    
  return 0;
}
