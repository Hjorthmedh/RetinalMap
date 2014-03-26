#include <stddef.h>
#include <math.h>
#include "mex.h"

/* these 2 #define lines help make the later code more readable */
/* Input Arguments */
#define S0  prhs[0]
#define TA0 prhs[1]
#define TB0 prhs[2]

/* square macro */
#define SQR(x) (x)*(x)

/* #define LITTLE_ENDIAN  1 */

static union 
{
  double d;
  struct {
#ifdef LITTLE_ENDIAN
    int j,i;
#else 
    int i,j;
#endif
  } n;
} _eco;

#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)

void mexexp(double*y, double*yp, size_t m) {
  while(m--) {
    *yp++ = EXP(*y++);
  }
}

void laplacian(double *conc, double *out, const mxArray *neighbours, int nSC) {
  /* const char *class_name; */
  int j, k, ik;
  mxArray *pneighbour;
  mwSignedIndex *neighbour;
  /* printf("laplacian nSC %d, %d\n", nSC, mxGetN(neighbours)); */
  /* int nmin =  100000; */
  /* int nmax = -100000; */

  for (j=0; j<nSC; j++) {
    /* printf("j %d ", j); */
    pneighbour = mxGetCell(neighbours, j);
    neighbour = (mwSignedIndex*)mxGetData(pneighbour);
    mxAssert(mxGetM(pneighbour) > 1, "Array not big enough");
    out[j] = 0;
    /* mexPrintf("mxGetM(pneigbour): %d, mxGetN(pneigbour): %d, k:",  */
    /*           mxGetM(pneighbour), mxGetN(pneighbour)); */
    /* class_name = mxGetClassName(pneighbour); */
    /* mexPrintf("Class Name: %s%s\n", class_name, */
    /*             mxIsSparse(pneighbour) ? " (sparse)" : ""); */

    for (k=0; k < mxGetM(pneighbour); k++) {
      ik = (int)neighbour[k] - 1;
      /* mexPrintf(" %d:%d", k, ik); */
      /* if (ik < nmin) nmin = ik;  */
      /* if (ik > nmax) nmax = ik;  */
      out[j] += conc[ik] - conc[j];
    }

    /* mexPrintf("\n"); */
  }
  /* if ((nmin != 0) || (nmax != nSC - 1)) {  */
  /*   mexPrintf("nmin: %d, nmax: %d, nSC: %d\n", nmin, nmax, nSC);   */
  /* } */
}

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray*prhs[] )
{
  int n, i, j;

  /* Check for proper number of arguments */    
  if (nrhs != 13) { 
    mexErrMsgTxt("13 input arguments required."); 
  } else if (nlhs != 5) {
    mexErrMsgTxt("Three output arguments required.");
  }

  /* Input arguments */
  size_t nRGC = (size_t)mxGetM(S0); 
  size_t nSC  = (size_t)mxGetN(S0);
  double *RA  = mxGetPr(prhs[3]);
  double *RB  = mxGetPr(prhs[4]);
  const mxArray *neighbours = prhs[5];
  double alpha = mxGetScalar(prhs[6]);
  double beta  = mxGetScalar(prhs[7]);
  double gamma = mxGetScalar(prhs[8]);
  double kappa = mxGetScalar(prhs[9]);
  double RGCEphAScale = mxGetScalar(prhs[10]);
  double dt    = mxGetScalar(prhs[11]);
  double N     = mxGetScalar(prhs[12]);

  /* Create matricies for the return arguments */ 
  plhs[0] = mxDuplicateArray(S0); 
  double *S = mxGetPr(plhs[0]); 
  plhs[1] = mxDuplicateArray(TA0); 
  double *TA = mxGetPr(plhs[1]);
  plhs[2] = mxDuplicateArray(TB0); 
  double *TB = mxGetPr(plhs[2]);

  mxArray *pIA =  mxCreateDoubleMatrix(nSC, 1, mxREAL);
  plhs[3] = pIA;
  mxArray *pIB =  mxCreateDoubleMatrix(nSC, 1, mxREAL);
  plhs[4] = pIB;
  mxArray *pout = mxCreateDoubleMatrix(nSC, 1, mxREAL);
  double *IA = mxGetPr(pIA);
  double *IB = mxGetPr(pIB);
  double *out = mxGetPr(pout);

  double Psi, Phi;
  double Si;                    /* Sum from each RGC */
  double Sj;                    /* Sum to each SC cell */
  /* Do the actual computation*/
  for (n=0; n<N; n++) {
    /* mexPrintf("%d\n", n); */
    /* mexPrintf("nSC: %d\n", nSC); */
    /* mexPrintf("nRGC: %d\n", nRGC); */
    
    /* Compute induced marker */
    for (j=0; j<nSC; j++) {
      /* mexPrintf("j: %d\n", j); */
      /* Total weight onto each SC cell */
      Sj = 0;
      for (i=0; i<nRGC; i++) {
      /* mexPrintf("i: %d\n", i); */
        /* Column-first arrangement */
        Sj += S[i + j*nRGC]; 
      }
      if (Sj == 0) {
        mexPrintf("j = %d; Sj = 0\n", j);
      }
      /* Post-synaptic divisive normalisation */
      IA[j] = 0;
      IB[j] = 0;
      for (i=0; i<nRGC; i++) {
        IA[j] += S[i + j*nRGC]*RA[i]/Sj; 
        IB[j] += S[i + j*nRGC]*RB[i]/Sj; 
      }
    }
    /* Update target marker */
    /* BUG AROUND HERE */
    /* for (j=0; j<nSC; j++) out[j]=0; */
    laplacian(TA, out, neighbours, nSC);
    for (j=0; j<nSC; j++) {
      TA[j] += dt*(alpha*(1 - RGCEphAScale*IA[j]*TA[j]) + beta*out[j]);
      if (isnan(TA[j])) {
        mexPrintf("IA[%d]=%1.5f; out[%d]=%1.5f\n", j, IA[j], j, out[j]);
      }
      mxAssert(!isnan(TA[j]), "NaN in TA");
    }
    laplacian(TB, out, neighbours, nSC);
    for (j=0; j<nSC; j++) {
      TB[j] += dt*(alpha*(IB[j] - TB[j]) + beta*out[j]);
    }

    /* Update synaptic weights */
    for(i=0; i<nRGC; i++) {
      Si = 0;
      for(j=0; j<nSC; j++) {
        Psi = (SQR(RGCEphAScale*RA[i]*TA[j] - 1) + SQR(RB[i] - TB[j]))/(2*SQR(kappa));
        /* EXP would be Cheaty but fast EXP */
        /* if (Psi > 700) { */
        /*   Psi = 700; */
        /* } */
        /* Phi = EXP(-Psi); */
        Phi = exp(-Psi);
        /* if (Phi < 0) { */
        /*   mexPrintf("Phi = %f < 0: Psi = %f\n", Phi, Psi); */
        /* } */
        S[i + j*nRGC] += dt*gamma*Phi;
        Si += S[i + j*nRGC];
      }
      /* Pre-synaptic divisive normalisation */
      for(j=0; j<nSC; j++) {
        S[i + j*nRGC] /= Si;
      }
    }

}
  return;
} 
