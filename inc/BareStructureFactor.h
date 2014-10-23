#ifndef GUARD_BARESTRUCTUREFACTOR_H
#define GUARD_BARESTRUCTUREFACTOR_H

#include <gsl/gsl_integration.h>
#include <vector>
#include "utable.h"

#define PI 3.1415926535897932384626433832795

class BareStructureFactor {
public:

private:
  double Kc, B, D, T;
  double avgLr, avgMz;
  Utable utable;  
  // should be called after parameters are changed. considered as updater
  void avgLrChanged();
  void avgMzChanged();
  void set_eta_lambda();
};

#endif
