/* This is a MEX-version of the step function 
**
** First argument is the object handle
** Second argument is the number of steps to take
**/

#include "mex.h"
#include <math.h>
#include <string.h>

/* Global variables, to avoid excessive argument passing */

/* Chemical interaction parameters and variables */
int     typeFlag;
double  alphaForwardChem;
double  betaForwardChem;
double  alphaReverseChem;
double  betaReverseChem;
double  alphaServoChem;
double  betaServoChem;
double* RGCEphA;
double* RGCEphB;
double* SCephrinA;
double* SCephrinB;
double* RGCephrinA;
double* RGCephrinB;
double* SCEphA;
double* SCEphB;
double  servoExp;

/* Competition parameters and variables */
double  AComp;
double  BComp;
double  DComp;
double  EComp;
double  alphaComp;
double  betaComp;
double  deltaComp;
double  epsilonComp;

/* Activity parameters and variables */
double *CAct;
mxArray *UAct;
double gammaAct;
double actScaling;

/* General parameters */
int *presynapticConnections;
double *presynapticWeight;
double *numPresynapticConnections;
double *totalWeightRGC;
double *totalWeightSC;

/* mxArray *neighbourRGC; */
mxArray *neighbourSC;
/* int *nNeighbourRGC; */
int *nNeighbourSC;
int maxConnections;
int nRGC, nSC;

int totalSynapsesInSystem;
double *time;

/* Function definitions */

int getRandomSCWeighted();

int getRandomNewRGC();

int getRandomExistingConnection(int SCidx);

int randDecision(double deltaEnergy);

void addSynapse(int SCidx, int RGCidx);

void removeSynapse(int SCidx, int RGCidx);

mxArray *setupOutputDouble(const mxArray *inputPtr, const char *name);
mxArray *setupOutputInt(const mxArray *inputPtr, const char *name);

double calculateSynapseAdditionEnergy(int SCidx, int RGCidx);
double calculateSynapseRemovalEnergy(int SCidx, int RGCidx);

void step(int nSteps);

void printParameterDump();
void testAddRemovalSymmetry();

/**************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] ) {

  int nSteps, i; 

  /*** Set up the output arguments data structures ***/

  if(nlhs != 6) {
    mexErrMsgTxt("stepFast: Incorrect number of return arguments");
  }
  
  /* First LHS argument, synapse connections */
  plhs[0] = setupOutputInt(prhs[0],"presynapticConnections");
  presynapticConnections = (int*) mxGetData(plhs[0]);
  
  /* Second LHS argument, synapse weight */
  plhs[1] = setupOutputDouble(prhs[0],"presynapticWeight");
  presynapticWeight = mxGetPr(plhs[1]);
  
  /* Third LHS argument, number of connections */
  plhs[2] = setupOutputDouble(prhs[0],"numPresynapticConnections");
  numPresynapticConnections = mxGetPr(plhs[2]);
  
  /* Fourth LHS argument, number of targets */
  plhs[3] = setupOutputDouble(prhs[0],"totalWeightRGC");
  totalWeightRGC = mxGetPr(plhs[3]);
  
  /* Fifth LHS argument, total weight */
  plhs[4] = setupOutputDouble(prhs[0],"totalWeightSC");
  totalWeightSC = mxGetPr(plhs[4]);
  
  /* Sixth LHS argument, time */
  plhs[5] = setupOutputDouble(prhs[0],"time");
  time = mxGetPr(plhs[5]);

  /*** Read in the input parameters (prhs) ***/
  
  if(nrhs == 0) {
    mexErrMsgTxt("stepFast: First argument object handle, second argument optional, nSteps");
  }
  
  if(!mxIsClass(prhs[0],"RetinalMap")) {
    mexErrMsgTxt("stepFast: First argument should be a RetinalMap object");
  }
  
  if(nrhs > 1) {
    nSteps = (int) mxGetScalar(prhs[1]);
  } else {
    nSteps = (int) (mxGetScalar(mxGetProperty(prhs[0],0,"nSteps")) 
		    - mxGetScalar(mxGetProperty(prhs[0],0,"curStep")));
  }
  
  /*** Read additional data from the RetinalMap object ***/
  
  /* neighbourRGC = mxGetProperty(prhs[0],0,"neighbourRGC"); */
  neighbourSC = mxGetProperty(prhs[0],0,"neighbourSC");
  
  /* nNeighbourRGC = mxGetData(mxGetProperty(prhs[0],0,"nNeighbourRGC")); */
  nNeighbourSC = mxGetData(mxGetProperty(prhs[0],0,"nNeighbourSC"));
  
  maxConnections = (int) mxGetScalar(mxGetProperty(prhs[0],0,"maxConnections"));
  nRGC = (int) mxGetScalar(mxGetProperty(prhs[0],0,"nRGC"));
  nSC = (int) mxGetScalar(mxGetProperty(prhs[0],0,"nSC"));
  
  /* Activity */
  CAct = mxGetPr(mxGetProperty(prhs[0],0,"CAct"));
  UAct = mxGetProperty(prhs[0],0,"UAct");
  gammaAct = mxGetScalar(mxGetProperty(prhs[0],0,"gammaAct"));
  actScaling = mxGetScalar(mxGetProperty(prhs[0],0,"actScaling"));
  
  /* Competition */
  AComp = mxGetScalar(mxGetProperty(prhs[0],0,"AComp"));
  BComp = mxGetScalar(mxGetProperty(prhs[0],0,"BComp"));
  DComp = mxGetScalar(mxGetProperty(prhs[0],0,"DComp"));
  EComp = mxGetScalar(mxGetProperty(prhs[0],0,"EComp"));
  
  alphaComp = mxGetScalar(mxGetProperty(prhs[0],0,"alphaComp"));
  betaComp = mxGetScalar(mxGetProperty(prhs[0],0,"betaComp"));
  deltaComp = mxGetScalar(mxGetProperty(prhs[0],0,"deltaComp"));
  epsilonComp = mxGetScalar(mxGetProperty(prhs[0],0,"epsilonComp"));
  
  /* Chemical */
  RGCEphA = mxGetPr(mxGetProperty(prhs[0],0,"RGCEphA"));
  RGCEphB = mxGetPr(mxGetProperty(prhs[0],0,"RGCEphB"));
  
  SCephrinA = mxGetPr(mxGetProperty(prhs[0],0,"SCephrinA"));
  SCephrinB = mxGetPr(mxGetProperty(prhs[0],0,"SCephrinB"));
  
  RGCephrinA = mxGetPr(mxGetProperty(prhs[0],0,"RGCephrinA"));
  RGCephrinB = mxGetPr(mxGetProperty(prhs[0],0,"RGCephrinB"));
  SCEphA = mxGetPr(mxGetProperty(prhs[0],0,"SCEphA"));
  SCEphB = mxGetPr(mxGetProperty(prhs[0],0,"SCEphB"));
  
  typeFlag = (int) mxGetScalar(mxGetProperty(prhs[0],0,"typeFlag"));
  
  alphaForwardChem = mxGetScalar(mxGetProperty(prhs[0],0,"alphaForwardChem"));
  betaForwardChem = mxGetScalar(mxGetProperty(prhs[0],0,"betaForwardChem"));
  alphaReverseChem = mxGetScalar(mxGetProperty(prhs[0],0,"alphaReverseChem"));
  betaReverseChem = mxGetScalar(mxGetProperty(prhs[0],0,"betaReverseChem"));
  alphaServoChem = mxGetScalar(mxGetProperty(prhs[0],0,"alphaServoChem"));
  betaServoChem = mxGetScalar(mxGetProperty(prhs[0],0,"betaServoChem"));
  
  servoExp = mxGetScalar(mxGetProperty(prhs[0],0,"servoExp"));
  
  
  /* Count the total number of synapses */
  totalSynapsesInSystem = 0;
  for(i = 0; i < nSC; i++) {  
    totalSynapsesInSystem += totalWeightSC[i];
  }
  
  /* Test functions */
  
  /* printParameterDump(); */
  testAddRemovalSymmetry();
  
  /*  mexPrintf("%d\n", mxGetClassID(prhs[0]));
      mexPrintf("%s\n", mxGetClassName(prhs[0])); */

  step(nSteps);

  mexPrintf("The system has %d synapses total\n", totalSynapsesInSystem);

}


/***************************************************************************/

/* We have to watch out, stored indexes are in matlab convention, but when
   we call functions from MEX, we use C-style indexing. The wrapper function
   stepMex.m takes care of this translation.
*/

int getRandomNewRGC() {

  /* By doing floor we get indexes starting from 0 */
  int RGCidx = (int) floor(nRGC*(rand()/((double) RAND_MAX + 1)));

  return RGCidx;

}

/**************************************************************************/

int getRandomSCWeighted() {

  double r;
  int wSum, SCidx;

  /* We want to weight likelyhood to choose a SC with the number of
     synapses it has. This so we choose any synapse with equal probability */

  r = totalSynapsesInSystem*(rand()/((double) RAND_MAX + 1));
  wSum = 0;
  SCidx = 0;

  while(wSum + totalWeightSC[SCidx] < r) {
    wSum += totalWeightSC[SCidx];
    SCidx += 1;
  }
  
  /* mexPrintf("r = %f, set SCidx = %d\n", r, SCidx); */

  return SCidx;

}

/**************************************************************************/

int getRandomExistingConnectionUnweighted(int SCidx) {

  /* uses global variables presynapticConnections, 
     numerOfConnections and maxConnections */

  /* No range check on SCidx, so watch out ! */

  double nSyn = numPresynapticConnections[SCidx];

  if(nSyn == 0) {
    /* mexPrintf("SC neuron %d does not have any synapses\n", SCidx); */
    return -1;
  }

  /* By doing floor we get indexes starting from 0 */
  int pos = (int) floor(nSyn*(rand()/((double) RAND_MAX + 1)));

  return presynapticConnections[SCidx*maxConnections+pos];

}

/***********************************************************************/

int getRandomExistingConnection(int SCidx) {

  double nSyn = numPresynapticConnections[SCidx];
  int offset;
  double r, wSum;

  if(nSyn == 0) {
    /* mexPrintf("SC neuron %d does not have any synapses\n", SCidx); */
    return -1;
  }

  r = totalWeightSC[SCidx]*(rand()/((double) RAND_MAX + 1));

  offset = SCidx*maxConnections;
  wSum = 0;

  while(wSum + presynapticWeight[offset] < r) {
    wSum += presynapticWeight[offset];
    offset += 1;
  }

  /* mexPrintf("SCidx=%d, r=%f, picked RGCidx=%d\n",
     SCidx,r,presynapticConnections[offset]); */
  
  return presynapticConnections[offset];    

}

/**************************************************************************/

int randDecision(double deltaEnergy) {

  return ((rand()/((double) RAND_MAX + 1)) < 1/(1 + exp(4*deltaEnergy)));

}

/**************************************************************************/

/* Dangerous function without range checks, 
   you better know what you are doing!      */

void addSynapse(int SCidx, int RGCidx) {

  /* uses global variables presynapticConnections, presynapticWeight,
     numPresynapticConnections, totalWeightRGC, totalWeightSC
     maxConnections */

  if(SCidx >= nSC) {  
    mexErrMsgTxt("removeSynapse: SCidx out of range\n");        
  }
  
  /* Check to see if a connection already exists from the SC to the RGC */

  int nCon = (int) numPresynapticConnections[SCidx];
  int foundIt = 0;
  int startOfs = SCidx*maxConnections;
  int i;

  /* mexPrintf("nCon = %d, maxConnections = %d, startOfs = %d\n",
     nCon, maxConnections, startOfs); */

  for(i = 0; i < nCon; i++) {
    if(presynapticConnections[startOfs+i] == RGCidx) { 

      /* mexPrintf("addSynapse: Found connection at i=%d\n", i); */
      foundIt = 1;

      /* A connection already exists, increment weight */
      presynapticWeight[startOfs+i] += 1.0;

      break;
    }
  }
	
  if(!foundIt) {
    /* Synapse does not exist, add it to table */

    /* mexPrintf("Could not find existing synapse\n"); */

    if(nCon >= maxConnections) {
      mexErrMsgTxt("addSynapse: Raise maxConnections, or set it to NaN for auto.");
    }

    /* Append a new synapse at end of table */
    presynapticConnections[startOfs + nCon] = RGCidx; 
    presynapticWeight[startOfs + nCon] = 1.0;

    numPresynapticConnections[SCidx] += 1;
  } 

  totalWeightRGC[RGCidx] += 1;
  totalWeightSC[SCidx] += 1;
  totalSynapsesInSystem += 1;

  /* mexPrintf("Number of %d targets %f\n", RGCidx, totalWeightRGC[RGCidx]); */


}

/**************************************************************************/

/* Dangerous function without range checks, 
   you better know what you are doing!      */

void removeSynapse(int SCidx, int RGCidx) {

  /* uses global variables presynapticConnections, presynapticWeight, 
     numPresynapticConnections, totalWeightRGC, totalWeightSC, maxConnections */

  if(SCidx >= nSC) {  
    mexErrMsgTxt("removeSynapse: SCidx out of range\n");        
  }
    
  int nCon = (int) numPresynapticConnections[SCidx];
  int foundIt = 0;
  int startOfs = SCidx*maxConnections;
  int lastIdx = startOfs+nCon-1;
  int i;

  for(i = 0; i < nCon; i++) {
    if(presynapticConnections[startOfs+i] == RGCidx) {

      /* mexPrintf("removeSynapse: Found connection at i=%d\n", i); */

      foundIt = 1;

      presynapticWeight[startOfs+i] -= 1.0;

      if(presynapticWeight[startOfs+i] <= 0) {
	/* Zero-weight synapse, remove it.
	   We do so by copying the last element to this position, 
	   then clearing last element */

	/* mexPrintf("Removing synapse %d\n", i); */

	presynapticConnections[startOfs+i] = presynapticConnections[lastIdx];
	presynapticWeight[startOfs+i] = presynapticWeight[lastIdx];

	presynapticConnections[lastIdx] = -1;
	presynapticWeight[lastIdx] = 0.0;

	numPresynapticConnections[SCidx] -= 1;

      }

      totalWeightSC[SCidx] -= 1;
      totalWeightRGC[RGCidx] -= 1;
      totalSynapsesInSystem -= 1;

      break;
    }
  }

  if(!foundIt) {
    mexErrMsgTxt("removeSynapse: Trying to remove a non-existing synapse!");
  }

}

/**************************************************************************/

mxArray *setupOutputDouble(const mxArray *inputPtr, const char *name) {

  mxArray *mxInPtr; 
  mxArray *mxOutPtr;
  double *dataPtr;

  mxInPtr = mxGetProperty(inputPtr,0,name);

  if(!mxIsDouble(mxInPtr)) {
    mexErrMsgTxt("setupOutputDouble: variable must be double\n");
  }

  mxOutPtr = mxCreateNumericMatrix(mxGetM(mxInPtr),
				   mxGetN(mxInPtr),
				   mxDOUBLE_CLASS,
				   mxREAL);
  
  if(mxOutPtr == NULL) {
    mexErrMsgTxt("setupOutputDouble: Failed to allocate memory");
  }

  dataPtr = mxGetPr(mxOutPtr);

  memcpy(dataPtr, 
	 mxGetPr(mxInPtr), 
	 mxGetNumberOfElements(mxInPtr)*sizeof(double));

  return mxOutPtr;

}

/****************************************************************************/

mxArray *setupOutputInt(const mxArray *inputPtr, const char *name) {

  mxArray *mxInPtr;
  mxArray *mxOutPtr;
  int *dataPtr;

  mxInPtr = mxGetProperty(inputPtr,0,name);

  if(strcmp(mxGetClassName(mxInPtr),"int32") != 0) {
    mexErrMsgTxt("setupOutputInt: Variables must be int32");
  }

  mxOutPtr = mxCreateNumericMatrix(mxGetM(mxInPtr),
				   mxGetN(mxInPtr),
				   mxINT32_CLASS,
				   mxREAL);
  if(mxOutPtr == NULL) {
    mexErrMsgTxt("setupOutputInt: Failed to allocate memory");
  }
  
  dataPtr = mxGetData(mxOutPtr);

  memcpy(dataPtr, 
	 mxGetData(mxInPtr), 
	 mxGetNumberOfElements(mxInPtr)*sizeof(int));

  return mxOutPtr; 

}

/****************************************************************************/

/* Calling with all these arguments might take a lot of time? */

/**** MAKE THE VARIABLES GLOBAL ***/

double calculateSynapseAdditionEnergy(int SCidx, int RGCidx) {

  double nAxon = totalWeightRGC[RGCidx];
  double nDend = totalWeightSC[SCidx];

  double Ecomp;
  double Echem;
  double Eact;

  double *UActLocal;
  int ofsC;

  mxArray *SCneigh;
  int *SCneighIdx;
  int numSCneighIdx;
  int i, j;
  int nNeighCon, nSelfCon;
  double tmpCW;
  int startOfs;
  int neighCon, selfCon;
  double neighConWeight, selfConWeight;
  int nCon;
  double w;

  /* Competition energy */

  Ecomp = AComp*(pow(nAxon+1,alphaComp)-pow(nAxon,alphaComp))
    + BComp*(pow(nAxon+1,betaComp) -pow(nAxon,betaComp)) 
    + DComp*(pow(nDend+1,deltaComp)-pow(nDend,deltaComp))
    + EComp*(pow(nDend+1,epsilonComp)-pow(nDend,epsilonComp));

  /* Chemical energy */


  switch(typeFlag) {
  case 1: /* Forward signaling only */
    Echem = alphaForwardChem*RGCEphA[RGCidx]*SCephrinA[SCidx]
      - betaForwardChem*RGCEphB[RGCidx]*SCephrinB[SCidx];
    break;
  case 2: /* Forward and reverse signaling */
    Echem = alphaForwardChem*RGCEphA[RGCidx]*SCephrinA[SCidx]
      - betaForwardChem*RGCEphB[RGCidx]*SCephrinB[SCidx]
      + alphaReverseChem*RGCephrinA[RGCidx]*SCEphA[SCidx]
      - betaReverseChem*RGCephrinB[RGCidx]*SCEphB[SCidx];

    break;
  case 3: /* Servomechanism */
    Echem = alphaServoChem*fabs(RGCEphA[RGCidx]*SCephrinA[SCidx]-servoExp)
      + betaServoChem*fabs(RGCEphB[RGCidx]-SCephrinB[SCidx]);

    break;
  default:
    mexErrMsgTxt("stepFast: Unknown type flag!");
    Echem = 0;
    break;

  }

  /* Activity energy */

  Eact = 0.0;

	/* Expensive computation, only do if we have activity component */
  if(gammaAct > 0) {

		/* Loop through all neighbouring SC (excluding self) */
		SCneigh = mxGetCell(neighbourSC,SCidx);
		SCneighIdx = mxGetData(SCneigh);

		UActLocal = mxGetData(mxGetCell(UAct,SCidx));
		ofsC = RGCidx*nRGC;

		for(i = 0; i < nNeighbourSC[SCidx]; i++) {
		
			/* Estimate the neighbours activity by all their incoming RGC synapses */
			nNeighCon = numPresynapticConnections[SCneighIdx[i]];

			tmpCW = 0;
			startOfs = SCneighIdx[i]*maxConnections;


			for(j = 0; j < nNeighCon; j++) {
				/* Sum their contribution, according to weight */
				neighCon = presynapticConnections[startOfs+j];
				neighConWeight = presynapticWeight[startOfs+j];
				tmpCW += CAct[ofsC + neighCon]*neighConWeight;
			}

			/* How far is this SC neighbour from our SC neuron, scale by U func */
			Eact += UActLocal[i] * tmpCW;

		}

		/* Scaling correlations with external cells, since those numbers
			 are affected by increase in neurons */

		Eact *= actScaling;

		nSelfCon = numPresynapticConnections[SCidx]; /* Our SC synapses */

		tmpCW = 0;
		startOfs = SCidx*maxConnections;

		for(j = 0; j < nSelfCon; j++) {
			/* Sum their contribution, according to weight */
			selfCon = presynapticConnections[startOfs+j];
			selfConWeight = presynapticWeight[startOfs+j];
			
			/* U = 1, for synapses on same SC */
			tmpCW += CAct[ofsC + selfCon]*selfConWeight;
		}
		
		/* Since the synapse is not added yet, the correlation with
			 itself was not taken into account, add that term. 
			 Note that there might also be zero synapes between that SC and RGC */
		
		w = 1;
    
		for(j = 0; j < nSelfCon; j++) {
			selfCon = presynapticConnections[startOfs+j];
			
			if(RGCidx == selfCon) {
				/* The synapse we are adding should be included 
					 This is a bit tricky, so we want change as we go from
					 n^2 to (n+1)^2, ie increase of 2*n+1, but we have already
					 included n in the previous look, so just add n+1 */
				w = presynapticWeight[startOfs+j] + 1.0;
				/* mexPrintf("Add found, weight (%d,%d), %f\n", SCidx,RGCidx, w); */
				break;
			}
		}
    
		tmpCW += w * CAct[0]; /* C = 1, for self-comparisons */
    
		Eact += tmpCW;
	
		Eact *= (-gammaAct); /* The 0.5 is removed, not double counting */

	}
	/* Eact was set to 0 if gammaAct is zero or negative */

  /* mexPrintf("Etot = %f, Echem = %f, Eact = %f, Ecomp = %f\n",
     Echem+Ecomp+Eact,Echem,Eact,Ecomp);
  */

  return Echem+Ecomp+Eact;

}


/****************************************************************************/

double calculateSynapseRemovalEnergy(int SCidx, int RGCidx) {

  double nAxon = totalWeightRGC[RGCidx];
  double nDend = totalWeightSC[SCidx];

  double Ecomp;
  double Echem;
  double Eact;

  double *UActLocal;
  int ofsC;

  mxArray *SCneigh;
  int *SCneighIdx;
  int numSCneighIdx;
  int i, j;
  int nNeighCon, nSelfCon;
  double tmpCW;
  int startOfs;
  int neighCon, selfCon;
  double neighConWeight, selfConWeight;
  int nCon;
  double w;

  Ecomp = AComp*(pow(nAxon-1,alphaComp)-pow(nAxon,alphaComp))
    + BComp*(pow(nAxon-1,betaComp) -pow(nAxon,betaComp)) 
    + DComp*(pow(nDend-1,deltaComp)-pow(nDend,deltaComp))
    + EComp*(pow(nDend-1,epsilonComp)-pow(nDend,epsilonComp));

  switch(typeFlag) {
  case 1: /* Forward signaling only */
    Echem = alphaForwardChem*RGCEphA[RGCidx]*SCephrinA[SCidx]
      - betaForwardChem*RGCEphB[RGCidx]*SCephrinB[SCidx];
    break;
  case 2: /* Forward and reverse signaling */
    Echem = alphaForwardChem*RGCEphA[RGCidx]*SCephrinA[SCidx]
      - betaForwardChem*RGCEphB[RGCidx]*SCephrinB[SCidx]
      + alphaReverseChem*RGCephrinA[RGCidx]*SCEphA[SCidx]
      - betaReverseChem*RGCephrinB[RGCidx]*SCEphB[SCidx];

    break;
  case 3: /* Servomechanism */
    Echem = alphaServoChem*fabs(RGCEphA[RGCidx]*SCephrinA[SCidx]-servoExp)
      + betaServoChem*fabs(RGCEphB[RGCidx]-SCephrinB[SCidx]);
    break;
  default:
    mexErrMsgTxt("stepFast: Unknown type flag");
    Echem = 0;
    break;

  }

  /* Activity energy */

  Eact = 0.0;

  /* Loop through all neighbouring SC (excluding self) */
  SCneigh = mxGetCell(neighbourSC,SCidx);
  SCneighIdx = mxGetData(SCneigh);
  UActLocal = mxGetData(mxGetCell(UAct,SCidx));

  ofsC = RGCidx * nRGC;

  for(i = 0; i < nNeighbourSC[SCidx]; i++) {
		
    /* Estimate the neighbours activity by all their incoming RGC synapses */
    nNeighCon = numPresynapticConnections[SCneighIdx[i]];

    tmpCW = 0;
    startOfs = SCneighIdx[i]*maxConnections;

    for(j = 0; j < nNeighCon; j++) {
      /* Sum their contribution, according to weight */
      neighCon = presynapticConnections[startOfs+j];
      neighConWeight = presynapticWeight[startOfs+j];
      tmpCW += CAct[ofsC + neighCon]*neighConWeight;
    }

    /* How far is this SC neighbour from our SC neuron, scale by U func */
    Eact += UActLocal[i] * tmpCW;

  }

  /* Scaling correlations with external cells, since those numbers
     are affected by increase in neurons */

  Eact *= actScaling;

  /* Contribution from own synapses, not to be scaled */
 
  nSelfCon = numPresynapticConnections[SCidx]; /* Our SC synapses */

  tmpCW = 0;
  startOfs = SCidx*maxConnections;

  for(j = 0; j < nSelfCon; j++) {
    /* Sum their contribution, according to weight */
    selfCon = presynapticConnections[startOfs+j];
    selfConWeight = presynapticWeight[startOfs+j];

    if(RGCidx == selfCon) {
      /* We go from n^2 to (n-1)^2, a decrease by -2*n+1 */
      selfConWeight = 2*selfConWeight - 1.0;
      /* mexPrintf("Remove found self (%d,%d), weight %f\n",
	 SCidx,RGCidx,selfConWeight); */
    }

    /* U = 1, for synapses on same SC */
    tmpCW += CAct[ofsC + selfCon]*selfConWeight;
  }
	
  Eact += tmpCW;
    
  Eact *= (-gammaAct);

  /*	mexPrintf("Etot = %f, Echem = %f, Eact = %f, Ecomp = %f\n",
	Ecomp-Echem-Eact,-Echem,-Eact,Ecomp);
  */

  /* Ecomp already took subtraction into account, but not Eact and Echem */
  return (Ecomp-Eact-Echem);

}

/****************************************************************************/

void printParameterDump() {

  mexPrintf("Visual parameter dump\n\n");

  mexPrintf("typeFlag = %d\n", typeFlag);
  mexPrintf("alphaForwardChem = %f\n", alphaForwardChem);
  mexPrintf("typeFlag = %f\n", typeFlag);
  mexPrintf("alphaFor = %f\n", alphaForwardChem);
  mexPrintf("betaForw = %f\n", betaForwardChem);
  mexPrintf("alphaRev = %f\n", alphaReverseChem);
  mexPrintf("betaReve = %f\n", betaReverseChem);
  mexPrintf("alphaSer = %f\n", alphaServoChem);
  mexPrintf("betaServ = %f\n", betaServoChem);
  mexPrintf("RGCEphA; = %f\n", RGCEphA);
  mexPrintf("RGCEphB; = %f\n", RGCEphB);
  mexPrintf("SCephrin = %f\n", SCephrinA);
  mexPrintf("SCephrin = %f\n", SCephrinB);
  mexPrintf("RGCephri = %f\n", RGCephrinA);
  mexPrintf("RGCephri = %f\n", RGCephrinB);
  mexPrintf("SCEphA = %f\n", SCEphA);
  mexPrintf("SCEphB = %f\n", SCEphB);
  mexPrintf("servoExp = %f\n\n", servoExp);

  mexPrintf("AComp = %f\n", AComp);
  mexPrintf("BComp = %f\n", BComp);
  mexPrintf("DComp = %f\n", DComp);
  mexPrintf("EComp = %f\n", EComp);
  mexPrintf("alphaComp = %f\n", alphaComp);
  mexPrintf("betaComp = %f\n", betaComp);
  mexPrintf("deltaComp = %f\n", deltaComp);
  mexPrintf("epsilonComp = %f\n\n", epsilonComp);

  mexPrintf("CAct[0] = %f\n", CAct[0]);
  double *U = mxGetData(mxGetCell(UAct,0));
  mexPrintf("UAct{0}(0) = %f\n\n", U);
  mexPrintf("gammaAct = %f\n", gammaAct);
  mexPrintf("actScaling = %f\n\n", actScaling);

  mexPrintf("presynapticConnections[0] = %f\n", presynapticConnections[0]);
  mexPrintf("presynapticWeight[0] = %f\n", presynapticWeight[0]);
  mexPrintf("numPresynapticConnections[0] = %f\n", numPresynapticConnections[0]);
  mexPrintf("totalWeightRGC[0] = %f\n", totalWeightRGC[0]);
  mexPrintf("totalWeightSC = %f\n\n", totalWeightSC[0]);

  mexPrintf("maxConnections = %f\n", maxConnections);
  mexPrintf("nSC = %d\n", nSC);
  mexPrintf("nRGC = %d\n", nRGC);

}

/**************************************************************************/

void step(int nSteps) {

  double E = 0;
  int i;
  int RGCidx, SCidx;

  mexPrintf("Running %d steps (even sampling of synapses)\n", nSteps);

  /* Step through the iterations */

  for(i = 0; i < nSteps; i++) {

    /* Pick a post synaptic SC at random */
    SCidx = (int) floor(nSC*(rand()/((double) RAND_MAX + 1)));

    RGCidx = getRandomNewRGC();

    if(RGCidx >= 0) {
      E = calculateSynapseAdditionEnergy(SCidx,RGCidx);

      if(randDecision(E)) {
				addSynapse(SCidx,RGCidx);
      }

    } else {
      mexPrintf("No RGC close to SC %d\n", SCidx);
    }

    /* Pick a new post synaptic SC at random, but this time weighted */
    /* SCidx = (int) floor(nSC*(rand()/((double) RAND_MAX + 1))); */
    SCidx = getRandomSCWeighted();

    RGCidx = getRandomExistingConnection(SCidx);

    if(RGCidx >= 0) {
      E = calculateSynapseRemovalEnergy(SCidx,RGCidx);

      if(randDecision(E)) {
				removeSynapse(SCidx,RGCidx);
      }
    }
    else {
      /* mexPrintf("No RGC synapsing on SC %d\n", SCidx); */
    }

		*time += 1.0/totalSynapsesInSystem;

  }

}

/**************************************************************************/

/* Test function to see that the energy to add/remove are symmetric */

void testAddRemovalSymmetry() {

  mexPrintf("Testing add/removal symmetry (global)\n");

  int SCidx = (int) floor(nSC*(rand()/((double) RAND_MAX + 1)));
  int RGCidx = getRandomNewRGC();
  double EA = calculateSynapseAdditionEnergy(SCidx,RGCidx);
  addSynapse(SCidx,RGCidx);
  double ER = calculateSynapseRemovalEnergy(SCidx,RGCidx);
  removeSynapse(SCidx,RGCidx);

  mexPrintf("(SC %d,RGC %d) EA = %f, ER = %f\n", SCidx,RGCidx,EA,ER);

  if(fabs(EA+ER) > 1e-9) {
    mexErrMsgTxt("Too large error!");
  }

  mexPrintf("Testing done.\n");
    
}

/****************************************************************************/
