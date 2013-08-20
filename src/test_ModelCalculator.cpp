#include "nfit.h"
#include "modelcalculator.h"
#include "utable.h"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001


void testCCDStrFct()
{
  ModelCalculator mc;
	mc.read_in_utable("../dat/utab_nfit12.15.dat");
	
  // Do not change the following parameters
	mc.setModelParameter(6e-13, "Kc");
	mc.setModelParameter(2e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0.1, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");
	
	vector<double> qvec, sfvec;
	mc.getCCDStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_CCD_struct_factor.dat");		
}


void testRotatedStrFct()
{
  ModelCalculator mc;
	mc.read_in_utable("../dat/utab_nfit12.15.dat");
	
  // Do not change the following parameters
	mc.setModelParameter(6e-13, "Kc");
	mc.setModelParameter(2e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0.1, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");
	
	vector<double> qvec, sfvec;
	mc.getRotatedStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_rotated_struct_factor.dat");		
}


void testMosaicStrFct()
{
	ModelCalculator mc;
	mc.read_in_utable("../dat/utab_nfit12.15.dat");
	
	// Do not change the following parameters
	mc.setModelParameter(6e-13, "Kc");
	mc.setModelParameter(2e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0.1, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");

	vector<double> qvec, sfvec;
	mc.getMosaicStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_mosaic_struct_factor.dat");
}


void testStrFct()
{
	ModelCalculator mc;
	mc.read_in_utable("../dat/utab_nfit12.15.dat");
	  
	mc.setModelParameter(6e-13, "Kc");
	mc.setModelParameter(2e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0.1, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");

	vector<double> qvec, sfvec;
	mc.getStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_struct_factor.dat");
}


// Run to create a new utable and test it
void testUtable()
{
	Utable utab;
	utab.writeUtableFile("utab_nfit12.15.dat", 600, 2000);
	utab.readUtableFile("utab_nfit12.15.dat");
}


int main()
{
  testStrFct();
  testMosaicStrFct();
  testRotatedStrFct();
  testCCDStrFct();

	return 0;
}

