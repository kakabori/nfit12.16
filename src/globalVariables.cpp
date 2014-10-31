#define PI 3.1415926535897932384626433832795

double xlow = 0;
double xhigh = 2048;
int nocz = 1;
int noscale = 1;
int kiyomask = 0;
double NKfilter = 0.0;
double dupe = 1;
double backgroundSigmaSquare = 10;
double aFactor = 0.2;
int updateSigma = 1;
int startflag = 0;
int stopflag = 0;
int nthreads = 2;

/*FunSupport parameters*/
double g_spRotatedAbserr = 1e-20;
double g_spRotatedRelerr = 1e-3;
double g_spRotatedMindx = 0.001;
double g_spRotatedMaxdx = 10;
double g_spMosaicAbserr = 1e-20;
double g_spMosaicRelerr = 1e-2;
double g_spMosaicMindx = 0.001;
double g_spMosaicMaxdx = 10;
double g_interp1dStrFctAbserr = 1e-20;
double g_interp1dStrFctRelerr = 1e-4;
double g_interp1dStrFctMindx = 0.001;
double g_interp1dStrFctMaxdx = 10;
double g_interp1dsumnAbserr = 1e-20;
double g_interp1dsumnRelerr = 1e-4;
double g_interp1dsumnMindx = 0.001;
double g_interp1dsumnMaxdx = 10;
double g_interp1dHrAbserr = 1e-20;
double g_interp1dHrRelerr = 1e-4;
double g_interp1dHrMindx = 0.01;
double g_interp1dHrMaxdx = 100000;

//For qag routine
//int g_workspaceSize = 10000; //workspace size
double g_epsabs = 0; //absolute error tolerance
double g_epsrel = 1e-4; //relative error tolerance for oscillatory integrand
double g_epsrel_low = 1e-4; // relative error tolerance for non-oscilattor integrand
//int g_key = 6; //key
