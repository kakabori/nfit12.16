// see fileTools.h for the implementation of genric functions
#include <vector>
#include <iostream>
#include <fstream>

using std::vector; using std::ofstream; using std::cerr;
using std::endl; using std::cout;


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
  myfile.close();
}


void saveMatrix(vector<double>& xv, vector<double>& yv, vector<double>& zv,
                const char *filename)
{
  cout << "xv.size(): " << xv.size() << " yv.size(): " << yv.size() 
       << " zv.size(): " << zv.size() << endl;
  if (xv.size()*yv.size() != zv.size()) {
    cerr << "\nzv.size() is not equal to xv.size()*yv.size()\n" << endl;
    return;
  }
  ofstream myfile;
  myfile.open(filename);
  typedef vector<double>::size_type vec_sz;
  for (vec_sz i = 0; i < yv.size(); i++) {
    for (vec_sz j = 0; j < xv.size(); j++) {
      myfile << xv[j] << " " << yv[i] << " " << zv[j+xv.size()*i] << endl;
    }
  }
  myfile.close();
}


void saveMatrix(unsigned int xs, unsigned int ys, vector<double>& zv,
                const char *filename)
{
  if (xs*ys != zv.size()) {
    cerr << "\nzv.size() is not equal to xs*ys\n" << endl;
    return;
  }
  ofstream myfile;
  myfile.open(filename);
  typedef vector<double>::size_type vec_sz;
  for (unsigned int i = 0; i < ys; i++) {
    for (unsigned int j = 0; j < xs-1; j++) {
      myfile << zv[j+xs*i] << " ";
    }
    myfile << zv[xs-1+xs*i] << endl;
  }
  myfile.close();
}


void saveThreeVectorsToFile(vector<double>& v1, vector<double>& v2,
                            vector<double>& v3, const char *filename)
{
  ofstream myfile;
  myfile.open(filename);
  typedef vector<double>::size_type vec_sz;
  for (vec_sz i = 0; i < v1.size(); i++) {
    myfile << v1[i] << " ";
  }
  myfile << endl;
  for (vec_sz i = 0; i < v2.size(); i++) {
    myfile << v2[i] << " ";
  }
  myfile << endl;
  for (vec_sz i = 0; i < v3.size(); i++) {
    myfile << v3[i] << " ";
  }
  cout << v1.size() << " " << v2.size() << " " << v3.size() << endl;
  myfile.close();
}
