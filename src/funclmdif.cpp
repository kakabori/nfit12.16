#include <sstream>
#include <tcl.h>
#include <tk.h>
#include <pthread.h>

#include "nfit.h"
#include "funclmdif.h"
#include "globalVariables.h"
#include "dataset.h"
#include "tcl_utility.h"
#include "fileTools.h"

extern void updatelinks(double xisquare, char *chain);
extern double refine(double theor, double exper);
extern Para g_ParaStruct;
extern Tcl_Interp *NKinterp;

//////////////////////////////////////////////////////////////
// compute I_e - I_{model} to be used for least square optimization
// lmdif which is linked upon compilation - see makefile
//////////////////////////////////////////////////////////////

void FuncLmdif::setPara(Para *p)
{
  int i;
  para = p;
  for(i = 0; i < Para::s_numParams; i ++) {
    bestParams[i] = p->getValue(i);
  }
}


/****************************************************************************************
This function logs the combination of parameters that yields the optimal
chi squared value, along with the corresponding chi squared value.
****************************************************************************************/
void FuncLmdif::logBest(double chisq, char *chain)
{
  int i;
  bestChisq = chisq;
  for(i = 0; i < Para::s_numParams; i ++) {
    bestParams[i] = para->getValue(i);
  }
  sprintf(bestChain, "Best values: %s", chain);
}


/****************************************************************************************
This function retirieves the combination of parameters that yielded the
optimal chi squared value, along with the corresponding chi squared value.
****************************************************************************************/
double FuncLmdif::recoverBestParams(Para *p)
{
  int i;
  for(i = 0; i < Para::s_numParams; i ++){
    p->setValue(bestParams[i], i);
  }
  return bestChisq;
}


/****************************************************************************************
This function returns a string reporting the optimal value for each
parameter.
****************************************************************************************/
char *FuncLmdif::getBestChain()
{
  return bestChain;
}


/**************************************************************************************************
typedef int (*LmdifFunc)(int, int, double *, double *, void *, void *);
Lmdif object gets linked to this function at the beginning of a fit

n < m
m: an input variable set to the number of functions (= number of data points)
n: an input variable set to the number of variables (= number of free parameters)
par: an array of length n. On input par must contain an initial estimate of the
solution vector. On output (output from Lmdif?) par contains the final estimate of the solution 
vector.
fvec: an output array of length m which contains the functions evaluated at the 
output x.
**************************************************************************************************/
int FuncLmdif::WrapperFunclmdif(int m, int n, double *par, double *fvec, void *client, void* ctrl)
{
  return ((FuncLmdif*)client)->funclmdif(m, n, par, fvec, ctrl);
}


/******************************************************************************
This function calculates the sum of squares

m: number of data points
n: number of parameters to be fitted
par: an array that stores initial values of the free parameters
fvec: an array to store the residuals (= weighted difference between the model 
      and a data point)
ctrl: ?
******************************************************************************/
int FuncLmdif::funclmdif(int m, int n, double *par, double *fvec, void* ctrl)
{
  if(stopflag) return 0;
  //printf("funclmdif ");
  appendStringToFile("log.txt", "funclmdif ");
  appendStringToFile("buffer.txt", "funclmdif ");
  char chain[256]={0};
  char cmd[256];
  size_t i;
  double sum1, sum2, sum3, sum4, sum5;
  double scale, bias;  
  double sum = 0;
  double chisq;
  size_t nslice = data->qs.size();
  vector<DataPoint>::iterator it;
  vector<DataPoint>::iterator itend;

  // loop through n free parameters and update
  // corresponding ModelCalculator parameters
  for (int k = 0; k < n; k++) {
    //double tmp;
    //tmp = par[k];

    // for variables other than bc2b, don't accept a negative value
		if (k != Var_bc2b) {
		  par[k] = fabs(par[k]);
		  //printf("%g ", par[k]);
		  sprintf(cmd, "%g ", par[k]);
		  appendStringToFile("log.txt", cmd);
		  appendStringToFile("buffer.txt", cmd);
		  
		  sprintf(chain, "%s %g ", chain, par[k]);
		  *(para->xp[k]) = par[k];
		  mc->setpara(par[k], para->idx[k]);
		// bc2b can be negative when the beam is below the detector edge
		} else {
			//printf("%g ", par[k]);
			sprintf(cmd, "%g ", par[k]);
		  appendStringToFile("log.txt", cmd);
		  appendStringToFile("buffer.txt", cmd);
		  
			sprintf(chain, "%s %g ", chain, par[k]);
			*(para->xp[k]) = par[k];
			mc->setpara(par[k], para->idx[k]);
		}
  }

  fflush(stdout);
  double qzmin = para->setup.getqz(data->qs[0].qz);
  itend = data->qs[nslice-1].sdata.end();
  double qxmax = para->setup.getqr( (itend-1)->qx );
  double qzmax = para->setup.getqz(data->qs[nslice-1].qz);
  //cout << "qzmin: " << qzmin << endl;
  //cout << "qxmax: " << qxmax << endl;
  //cout << "qzmax: " << qzmax << endl;
  // initiate building the 2D structure factor map
  mc->init(0, qxmax, qzmin, qzmax);
  
  // mc.init() calculates up to mosaic convoluted 2D map, so
  // big calculation has been finished. Now, do small one slice by slice
  // for each qz value  
  for (i = 0; i < nslice; i++) {
    it = data->qs[i].sdata.begin();
    itend = data->qs[i].sdata.end();
    double qxmax_tmp = para->setup.getqr((itend-1)->qx);
    double qz_tmp = para->setup.getqz(data->qs[i].qz);
    //cout << "qxmax_tmp: " << qxmax_tmp << endl;
    //cout << "qz_tmp: " << qz_tmp << endl;
    mc->qxSlice(qxmax_tmp, qz_tmp);
    sum1 = sum2 = sum3 = sum4 = sum5 = 0;

    /* Obtain linear parameters c(z) background adjustment and
	     the overall scale s(z) that appear in Eq. 5.1 and Eq. 2.1,
	     given fixed values of all the other parameters that determine Sccd.
	     Compute the summations for determining linear parameters,*/    
    for (; it != itend; it++) {
      double vcal;
      double sigma2;    
      // Notation in Numerical Recipes, chapter 15.2, Eq. (15.2.4)
      // vcal = x_i
      // sigma2 = sigma^2
      // sum1 = S_{xy}
      // sum2 = S_{xx}
      // sum3 = S_x
      // sum4 = S_y
      // sum5 = S
      double qx_tmp = para->setup.getqr(it->qx);
      vcal = mc->getCCDStrFct(qx_tmp);
      //cout << qx_tmp << " " << qz_tmp << " " << vcal << endl;
      sigma2 = it->sigma * it->sigma;
      sum1 += it->inte * vcal / sigma2;
      sum2 += vcal * vcal / sigma2;
      sum3 += vcal / sigma2;
      sum4 += it->inte / sigma2;
      sum5 += 1 / sigma2;          
    }
    
    // obtain 'scale' and 'bias' using sum
    if (nocz) {
      scale = sum1 / sum2;
      bias = 0;
    } else {
      double Delta_NR = sum5*sum2 - sum3*sum3;        // -Delta defined in Numerical Recipes
      scale = (sum1*sum5 - sum4*sum3) / Delta_NR;     // s(z) (or phi(z))
      bias = -(sum2*sum4 - sum3*sum1) / Delta_NR;     // c(z) parameter
      double sigma_scale = sqrt(sum5 / Delta_NR);     // uncertainty in s(z)
      double sigma_bias = sqrt(sum2 / Delta_NR);      // uncertainty in c(z)

      // store uncertainties in Data object
    	data->qs[i].sigma_scale = sigma_scale;
    	data->qs[i].sigma_bias = sigma_bias;
    }
    
    if (noscale) {
      scale = 1;
      bias = 0;
    }
    
    // store 'scale' and 'bias' in Data object
    data->qs[i].scale = scale;
    data->qs[i].bias  = bias;
    
    // compute (I_e - I_{model})/(\sigma_p^2) and write them to fvec
    // also update 'data' with model values
    it = data->qs[i].sdata.begin();    
    for (; it != itend; it++) {
      double tmp;
      it->cal = tmp = mc->getCCDStrFct( para->setup.getqr(it->qx) ) * scale - bias;
      //do not include outstanding points (refinement)
      if (NKfilter != 0 && refine(tmp, it->inte) > NKfilter) {tmp=0; m-=1;}
      else {tmp = ( tmp - it->inte ) / it->sigma;}  //looks like the command used.
      /*      else {tmp = ( tmp - it->inte ) / sqrt(backgroundSigmaSquare);}
	      gives much smaller chi2 b/c backgroundSigmaSquare is fixed at 10.
            it->sigma comes from dataset.cxx  */
      *(fvec++) = tmp;
      sum += tmp * tmp;
    }
  }
  
  //Finalize the chi squared caluclation and report current parameter values.
  chisq = sum/m;
  //printf("Xr: %g  / %d = %g  /  %g \n",sum,m, chisq,backgroundSigmaSquare);
  sprintf(cmd, "Xr: %g  / %d = %g  /  %g \n", sum, m, chisq, backgroundSigmaSquare);
  appendStringToFile("log.txt", cmd);
  appendStringToFile("buffer.txt", cmd);
  evalTclCommand("insert_buffer_file");
  fflush(stdout);
  sprintf(chain, "%s Xr: %g  / %d = %g ", chain, sum, m, chisq);
  updatelinks(chisq, chain);
  if(chisq < bestChisq) logBest(chisq, chain);
  // update associated image with 'data'
  data->writeimg();
  data->writefrm((char *) "frm.dat", para);
  
  return 0;
}
  
