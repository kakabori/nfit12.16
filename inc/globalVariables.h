#ifndef GUARD_GLOBALVARIABLES_H
#define GUARD_GLOBALVARIABLES_H

#include <tcl.h>

extern double xlow;
extern double xhigh;
extern int nocz;
extern int noscale;
extern int kiyomask;
extern double NKfilter;
extern double dupe;
extern double backgroundSigmaSquare;
extern double aFactor;
extern int updateSigma;
extern int startflag;
extern int stopflag;
extern int nthreads;

/*FunSupport parameters*/
extern double g_spRotatedAbserr;
extern double g_spRotatedRelerr;
extern double g_spRotatedMindx;
extern double g_spRotatedMaxdx;
extern double g_spMosaicAbserr;
extern double g_spMosaicRelerr;
extern double g_spMosaicMindx;
extern double g_spMosaicMaxdx;
extern double g_interp1dStrFctAbserr;
extern double g_interp1dStrFctRelerr;
extern double g_interp1dStrFctMindx;
extern double g_interp1dStrFctMaxdx;
extern double g_interp1dsumnAbserr;
extern double g_interp1dsumnRelerr;
extern double g_interp1dsumnMindx;
extern double g_interp1dsumnMaxdx;
extern double g_interp1dHrAbserr;
extern double g_interp1dHrRelerr;
extern double g_interp1dHrMindx;
extern double g_interp1dHrMaxdx;


//For qag routine
//extern int g_workspaceSize; //workspace size
extern double g_epsabs; //absolute error tolerance
extern double g_epsrel; //relative error tolerance
extern double g_epsrel_low;
//extern int g_key; //key

extern Tcl_AsyncHandler g_commandHandler;
extern Tcl_Interp *g_interp;
extern std::string g_string;

#endif
