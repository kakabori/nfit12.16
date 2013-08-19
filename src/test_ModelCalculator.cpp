#include "nfit.h"
#include "modelcalculator.h"
#include "utable.h"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001


void spsumn(ModelCalculator& mc, double qz)
{
	mc.setSliceParameter(qz);
	mc.QxSlice(0, 0.302);
	
	vector<double> r;
	vector<double> sumn;
	mc.get_spsumnPoints(r, sumn);
	
	vector<double>::iterator iter, iter2;
	iter2 = sumn.begin();
	int count = 0;
	for (iter = r.begin(); iter != r.end(); iter++) {
		cout << *iter << " " << exp(*iter) - SMALLNUM << " " << *iter2 << endl;
		iter2++;
		count++;
	}
	cout << count << endl;
}


void rotated(ModelCalculator& mc, double qz)
{
	vector<double> qxvec;
	vector<double> rotatedSFvec;
	mc.getRotatedSF(0.008, 0.3, qz, qxvec, rotatedSFvec);
	
	ofstream myfile;
	myfile.open("rotatedSF_qz0.45_new.dat");
	for (size_t i = 0; i < rotatedSFvec.size(); i++) {
		myfile << qxvec[i] << " " << rotatedSFvec[i] << endl;
		//cout << qxvec[i] << " " << rotatedSFvec[i] << endl;
	}
	myfile.close();

	vector<double> qx, rotatedSF;
	mc.get_spRotatedPoints(qx, rotatedSF);
	vector<double>::iterator iter, iter2;
	iter2 = rotatedSF.begin();
	int count = 0;
	for (iter = qx.begin(); iter != qx.end(); iter++) {
		cout << *iter << " " << exp(*iter) - SMALLNUM << " " << exp(*iter2) << endl;
		iter2++;
		count++;
	}
	cout << count << endl;

}


void testMosaic()
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

	mc.setSliceParameter(0.3);
	mc.QxSlice(0.00, 0.3);
	vector<double> qvec, sfvec;
	mc.get_spStrFctPoints(qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << exp(qvec[i]) - SMALLNUM << " " << exp(sfvec[i]) << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.3 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
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
	mc.setModelParameter(0, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");

	vector<double> qvec, sfvec;
	mc.getStrFct(0, 0.3, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.3 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_struct_factor.dat");
}


void CCD(ModelCalculator& mc, double qz)
{
	mc.setSliceParameter(qz);
	mc.QxSlice(0.00, 0.3);
	
	for (double qx = 0.00; qx < 0.3;){
		cout << qx << " " << mc.getCCDStrFct(qx) << endl;
		qx += 0.01;
	}
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

	return 0;
}

