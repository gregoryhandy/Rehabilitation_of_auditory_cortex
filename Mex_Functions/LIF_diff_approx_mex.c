/*
**
**  This file translates the Julia code used in
**  Litwin-Kumar et al., 2016 to mex
**
*/
#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"

/* returns a pseudorandom value drawn from the standard normal distribution */
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	int istim, T, NT, recstart; 
	int Npop, *Ncells, Ntot, m1, m2;
	double dt, *rates; 
	int *pinds;
	
	/* Synaptic parameters */
	double *p0, *J;
	double *tau_s;
	double *tau_m, *EL;
	double *vth, *vre, *tauref;
	
	/* Connectivity matrix */
	int *wind, *wipost;
	double *wstr;
	//	, *wext, *rext;
	
	/* indexing variables */
	int tt, cc, kk, i, j, pp;
	
   /******
    * Import variables from matlab
    * This is messy looking and is specific to mex.
    * Ignore if you're implementing this outside of mex.
    *******/
	T =  (int)mxGetScalar(prhs[0]);
	NT = (int)mxGetScalar(prhs[1]);
	recstart =  (int)mxGetScalar(prhs[2]);
	Npop = (int)mxGetScalar(prhs[3]);
	Ntot = (int)mxGetScalar(prhs[4]);
	
	/*
	  Needed for vectors of integers, which is important
	  since they will be used to access array locations
	  As a result, the mex call has been motified
	*/
	mxArray *temp[1];
	mxArray *lhs[1];
	temp[0] = (mxArray *) prhs[5];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	Ncells=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[5]);
	m2 = mxGetN(prhs[5]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("Ncells should be Npopx1 or 1xNpop");

	dt = mxGetScalar(prhs[6]);

	rates= mxGetPr(prhs[7]);
	m1 = mxGetM(prhs[7]);
	m2 = mxGetN(prhs[7]);
	if(m1!=Ntot && m2!=Ntot)
	    mexErrMsgTxt("rates should be 1xNtot or Ntotx1");

	/* load in synaptic parameters */
	p0 = mxGetPr(prhs[8]);
	m1 = mxGetM(prhs[8]);
	m2 = mxGetN(prhs[8]);
	if(m1!=Npop || m2!=Npop)
	    mexErrMsgTxt("p0 should be NpopxNpop");

	J = mxGetPr(prhs[9]);
	m1 = mxGetM(prhs[9]);
	m2 = mxGetN(prhs[9]);
	if(m1!=Npop || m2!=Npop)
	    mexErrMsgTxt("J should be NpopxNpop");

	tau_s = mxGetPr(prhs[10]);
	m1 = mxGetM(prhs[10]);
	m2 = mxGetN(prhs[10]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("taurise should be 1xNpop or Npopx1");

	tau_m= mxGetPr(prhs[11]);
	m1 = mxGetM(prhs[11]);
	m2 = mxGetN(prhs[11]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("tau_m should be 1xNpop or Npopx1");

	EL= mxGetPr(prhs[12]);
	m1 = mxGetM(prhs[12]);
	m2 = mxGetN(prhs[12]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("EL should be 1xNpop or Npopx1");

	vth= mxGetPr(prhs[13]);
	m1 = mxGetM(prhs[13]);
	m2 = mxGetN(prhs[13]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("vth should be 1xNpop or Npopx1");

	vre= mxGetPr(prhs[14]);
	m1 = mxGetM(prhs[14]);
	m2 = mxGetN(prhs[14]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("vre should be 1xNpop or Npopx1");

	tauref= mxGetPr(prhs[15]);
	m1 = mxGetM(prhs[15]);
	m2 = mxGetN(prhs[15]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("tauref should be 1xNpop or Npopx1");

	/*
	  Needed for vectors of integers, which is important
	  since they will be used to access array locations
	  As a result, the mex call has been motified
	*/
	temp[0] = (mxArray *) prhs[16];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	wind=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[16]);
	m2 = mxGetN(prhs[16]);
	if(m1!=(Ntot+1) && m2!=(Ntot+1))
	    mexErrMsgTxt("wind should be (Ntot+1)x1 or 1x(Ntot+1)");

	temp[0] = (mxArray *) prhs[17];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	wipost=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[17]);
	m2 = mxGetN(prhs[17]);

	wstr= mxGetPr(prhs[18]);
	m1 = mxGetM(prhs[18]);
	m2 = mxGetN(prhs[18]);

	temp[0] = (mxArray *) prhs[19];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	pinds=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[19]);
	m2 = mxGetN(prhs[19]);
	if(m1!=(Npop+1) && m2!=(Npop+1))
	    mexErrMsgTxt("pinds should be (Npop+1)x1 or 1x(Npop+1)");

	unsigned int iseed;
	iseed =  (int)mxGetScalar(prhs[20]);
	
	
	double *mu_vec, *var_vec;
	mu_vec = mxGetPr(prhs[21]);
	m1 = mxGetM(prhs[21]);
	m2 = mxGetN(prhs[21]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("mu_vec should be Npopx1 or 1xNpop");
	
	var_vec = mxGetPr(prhs[22]);
	m1 = mxGetM(prhs[22]);
	m2 = mxGetN(prhs[22]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("var_vec should be Npopx1 or 1xNpop");
	
	/* Simple "srand()" seed: just use "time()" */
	/* iseed = (unsigned int)time(NULL); */
	srand (iseed);
	
	/*	
	Code for mex printing
	mexPrintf("test \n");
	mexEvalString("drawnow;");
	*/
	
	/* NOT IMPORTED FROM MATLAB */
	/* state vectors */
	double *v, *lastSpike;
	int *whichpop;

	v = mxMalloc(Ntot * sizeof(double));
	lastSpike = mxMalloc(Ntot * sizeof(double));
	whichpop = mxMalloc(Ntot * sizeof(int));


	double *I_s;
	I_s = mxMalloc(Ntot * sizeof(double));

	int ns = 0;
	int maxns = round(T*Ntot*0.05);

	int *counts;
	counts = mxMalloc(Ntot * sizeof(int));

	/* Allocate output vector */
	double *rates_official;
	plhs[0] = mxCreateDoubleMatrix(1, Ntot, mxREAL);
	rates_official=mxGetPr(plhs[0]);


	double *times;
	double *tinds;
	plhs[1] = mxCreateDoubleMatrix(1, maxns, mxREAL);
	times=mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(1, maxns, mxREAL);
	tinds=mxGetPr(plhs[2]);
	
	
	/* initial conditions */
	for(i = 0; i<Ntot; i++){
		v[i] = -65;
		lastSpike[i] = -100;

		/* save whichpop neuron i is */
		pp = 1;
		int initial_vT = 0;
		while(initial_vT == 0){
			if(i < pinds[pp]){
				whichpop[i] = pp-1;
				initial_vT = 1;
			}else{
				pp++;
			}

			if(pp > (Npop+1) ){
				mexErrMsgTxt("vT initialization failed");
			}
		}

		I_s[i] = 0; /* for each cell, one set of ODEs for each presynaptic population */
	}
	
	/* spike counts and times */
	ns = 0;
	for(i =0; i< maxns; i++){
		times[i] = 0;
		tinds[i] = 0;
	}
	for(cc = 0; cc < Ntot; cc++){
		counts[cc] = 0;
	}
		
	double sigma_hat[3];
	int pc;
	for(pc = 0; pc<Npop; pc++){
		sigma_hat[pc] = var_vec[pc]*sqrt(tau_m[pc])/tau_s[pc]*sqrt(dt);
	}
	 
	/* simulation------------------------------------------------------------- */
	for(tt = 1; tt <= NT && ns<maxns; tt++){

		double t = dt*tt;

		/* update synaptic/adaptation parameters */
		for(cc=0; cc<Ntot; cc++){
			int pc = whichpop[cc];

			/* Note: all matrices are 1D arrays with conversion [i][j] = i + j*num_rows */
			I_s[cc] = I_s[cc]+dt*(mu_vec[pc]-I_s[cc])/tau_s[pc]+sigma_hat[pc]*randn(0,1);
			
			if(t > (lastSpike[cc]+tauref[pc])){ /* not in refractory period */
				v[cc] += dt*(-(v[cc]-EL[pc])+I_s[cc])/tau_m[pc];
			} /* end if(refractory) */
		}


		/* Deal with spikes */
		for(cc=0; cc<Ntot; cc++){
			int pc = whichpop[cc];

			/* spike occurred */
			if((v[cc] > vth[pc]) && (ns<maxns)){

				/* update v */
				v[cc] = vre[pc];

				/* record spike */
				lastSpike[cc] = t;
				ns = ns+1;
				times[ns] = t;
				tinds[ns] = cc;
				if(t > recstart){
					counts[cc] = counts[cc] + 1;
				}

				/* propagate spike */
				for(kk = wind[cc]; kk < (wind[cc+1]-1); kk++){
					int ipost = wipost[kk];
					int ppost = whichpop[ipost];

					I_s[ipost] += wstr[kk]*tau_m[pc]/tau_s[pc];
				}
			}
		} /* end deal with spikes */
	}/* end of loop over time */

	/* Issue a warning if max number of spikes reached */
	if(ns>=maxns)
	   mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");

	/* Calculate the rates */
	for(i = 0; i < Ntot; i++){
		rates[i] = 1000*(double)counts[i]/ (double)(T-recstart);
		rates_official[i] = rates[i];
	}
		
	/* Free malloc'ed variables */
	mxFree(v);
	mxFree(lastSpike);
	mxFree(whichpop);
	mxFree(counts);
	mxFree(I_s);
	
}