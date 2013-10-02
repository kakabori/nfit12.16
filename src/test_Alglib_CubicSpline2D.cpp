#include "alglib_interpolation.h"
#include <vector>
#include <iostream>
using namespace std;
int main()
{
  Alglib_CubicSpline2D mc;
  vector<double> qrv, qzv, sfv;
  qrv.resize(4);
  qzv.resize(4);
  sfv.resize(16);
  qrv[0] = 0.000;
  qrv[1] = 0.001;
  qrv[2] = 0.002;
  qrv[3] = 0.003;
  qzv[0] = 0.061;
  qzv[1] = 0.062;
  qzv[2] = 0.063;
  qzv[3] = 0.064;
  sfv[0] = 4.2033e+09;
  sfv[1] = 2.04038e+07;
  sfv[2] = 3.19547e+06;
  sfv[3] = 1.18481e+06; 
  sfv[4] = 4.27329e+09;
  sfv[5] = 2.09309e+07; 
  sfv[6] = 3.29908e+06; 
  sfv[7] = 1.22898e+06; 
  sfv[8] = 4.35387e+09;
  sfv[9] = 2.15203e+07; 
  sfv[10] = 3.4138e+06;
  sfv[11] = 1.27763e+06;
  sfv[12] = 4.44585e+09;
  sfv[13] = 2.2178e+07;
  sfv[14] = 3.54076e+06;
  sfv[15] = 1.33124e+06;
  mc.buildInterpolant(qrv, qzv, sfv);
  
  double qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evaluate(qr, 0.061) << " ";
    qr += 0.001;
  }
  cout << endl;
  qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evaluate(qr, 0.062) << " ";
    qr += 0.001;
  }
  cout << endl;
  qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evaluate(qr, 0.063) << " ";
    qr += 0.001;
  }
  cout << endl;
  qr = 0.000;
  for (int i = 0; i < 4; i++) {
    cout << mc.evaluate(qr, 0.064) << " ";
    qr += 0.001;
  }
  cout << endl;
  
  cout << mc.evaluate(0.0011, 0.063) << endl;
  
  return 0;
}
