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


void test1DStrFct()
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


void test2DStrFctWithoutMosaic()
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
	mc.setModelParameter(0, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");

	vector<double> qvec, sfvec;
	mc.get2DStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_2D_struct_factor_without_mosaic.dat");
}


void test2DStrFctWithMosaic()
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
	mc.get2DStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_2D_struct_factor_with_mosaic.dat");
}


void create2DStrFctWithMosaic()
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

	vector<double> qrvec, qzvec, sfvec;
	mc.getMosaicStrFct(0, 0.15, 0.05, 0.15, qrvec, qzvec, sfvec);

	saveMatrix(qrvec, qzvec, sfvec, "mosaic_0.1.ssg");
	saveMatrix(qrvec.size(), qzvec.size(), sfvec, "mosaic_0.1.dat");
}


void create2DStrFctWithoutMosaic()
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
	mc.setModelParameter(0.0, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");

	vector<double> qrvec, qzvec, sfvec;
	mc.getMosaicStrFct(0, 0.15, 0.05, 0.15, qrvec, qzvec, sfvec);

	saveMatrix(qrvec, qzvec, sfvec, "no_mosaic.ssg");
	saveMatrix(qrvec.size(), qzvec.size(), sfvec, "no_mosaic.dat");
}


// Run to create a new utable and test it
void testUtable()
{
	Utable utab;
	utab.writeUtableFile("utab_nfit12.15.dat", 600, 2000);
	utab.readUtableFile("utab_nfit12.15.dat");
}


void testConvolveMosaic()
{
  ModelCalculator mc;
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
	
  mc.read_struct_factor("structure_factor.dat");
  
  double qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evalStrFct(qr, 0.061) << " ";
    qr += 0.001;
  }
  cout << endl;
  qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evalStrFct(qr, 0.062) << " ";
    qr += 0.001;
  }
  cout << endl;
  qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evalStrFct(qr, 0.063) << " ";
    qr += 0.001;
  }
  cout << endl;
  qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evalStrFct(qr, 0.064) << " ";
    qr += 0.001;
  }
  cout << endl;

  cout << mc.evalStrFct(0.0011, 0.063) << endl;
  //cout << mc.evalStrFct(0.0010995, 0.0629904) << endl;
  //cout << mc.convolveMosaic(0.063, 1*PI/180) << endl;
}


int main()
{
  //test1DStrFct();
  //test2DStrFctWithoutMosaic();
  //test2DStrFctWithMosaic();
  //testMosaicStrFct();
  //testRotatedStrFct();
  //testCCDStrFct();
  //create2DStrFctWithoutMosaic();
  //create2DStrFctWithMosaic();
  testConvolveMosaic();
	return 0;
}

