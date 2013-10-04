#ifndef _FuncLmdif_h_
#define _FuncLmdif_h_

#include <tcl.h>
#include <float.h>
#include <vector>
#include <string.h>
#include <limits>

#include "modelcalculator.h"
#include "tvds.h"
#include "tvImg.h"
#include "dataset.h"

/*
  FuncLmdif collects Data, Para and ModelCalculator objects
  and do optimization
*/

class FuncLmdif {
public:
  FuncLmdif() {bestChisq = numeric_limits<double>::infinity();}
  static int WrapperFunclmdif(int, int, double *, double *, void *, void *);
  void setData(Data *indata) {data = indata;}
  void setMC(ModelCalculator *b) {mc = b;}
  void setPara(Para *);
  int npoint() {return data->n;}
  void logBest(double, char *);
  double recoverBestParams(Para *);
  char *getBestChain();
private:
  Data *data;
  Para *para;
  ModelCalculator *mc;
  double bestChisq;
  double bestParams[18];
  char bestChain[256];
  int funclmdif(int m,int n,double *par,double *fvec,void*);
  static void* modelCalcThread(void *args);
};

#endif
