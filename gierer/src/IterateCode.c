#include <R.h>						/* Needed for Rprintf() function */
#include <math.h>					/* Needed for exp() function */
#include <time.h>


/* iteratedb() and p0db() are the new versions of the code that work on the
 * database of gradients and cellular positions provided by Johannes' framework.
 */

void add_scale_arrays(double *array1, int *array2, int array_length,
		      double scalar1, double scalar2);
void add_scale_arrays2(double *array1, int *array2, int array_length,
		      double scalar1, double scalar2);
void make_new_update_sequence(int *update_seq, int seq_length,
			      int update_method, int *active_cones, int total_active);
double p0(int p0_fn, double alpha, double beta, double sigma, double tau,
	  double x, double y, double u, double v);

void printTime(int it);

double p0db(int p0_fn,
	    double rEA, double reA, double rEB, double reB,
	    double tEA, double teA, double tEB, double teB);

void
iterate(int *axons,		/* Integer of retinal axons */
	double *axon_positions,	/* Matrix of axon co-ordinates (u,v) on retina */
	int *total_cones, /* Integer of number of growth cones */				
	int *cone_axons, /* Vector of axon number for each growth cone */
	int *active_cones,	/* Vector of active growth cones */
	int *total_active,	/* Integer of active growth cones */

	int *tects,		/* Integer of tectum positions */
	double *tect_positions,	/* Matrix of tectum co-ordinates (x,y) */
	int *tect_neighbours_size, /* Integer of tect_neighbours matrix size */
	int *tect_neighbours,	/* Matrix of tectum neighbour indices */
	int *tect_status, /* Vector of active tectal positions 1-active 0-inactive */

	int *connections,	/* Connection matrix */
				
	double *r,		/* Inhibitory term, vector size N */
	double *epsilon,	/* Rate of growth of r */
	double *kappa,		/* Rate of decay of r */

	int *rho,		/* Local density, vector size N */

	int *p0_fn,		/* Vector of values describing p0 to use */
	double *alpha,		/* p0 function parameters */
	double *beta,
	double *sigma,
	double *tau,

	int *iterations,	/* Number of times to iterate */

	int *sequence, /* User provided update sequence. *sequence == 0 if no sequence provided */
	int *s_length, /* Length of user update sequence or length of calculated sequence */

	int *update_method,	/* Integer describing update method */
	int *update_many) /* Boolean, new update sequence each iteration */
{
  /* If there are no iterations, then leave */
  if (*iterations==0)
    return;

  /* Get the R environment, random number generation state */
  GetRNGstate();

  /* Set-up update sequence variables */
  int seq_length = *s_length;
  int update_seq[seq_length];

  /* Construct sequence in C code */
  if (*sequence == 0)
    make_new_update_sequence(&update_seq[0], seq_length, *update_method, active_cones, *total_active);
  /* Use user sequence */
  else
    for (int s = 0; s < seq_length; s++)
      update_seq[s] = *(sequence + s);

  /* Loop for each iteration */
  for (int it = 0; it < *iterations; it++)
    {
      /* Loop for each specified density unit in the update sequence*/		
      for (int unit = 0; unit < seq_length; unit++)
	{
	  /* Get a growth cone */
	  int growth_cone = update_seq[unit];

	  /* Get the axon */
	  int axon = *(cone_axons + growth_cone - 1);

	  /* Get the retinal co-ordinates of the axon */
	  double u = *(axon_positions + axon - 1);
	  double v = *(axon_positions + axon + *axons - 1);

	  /* Get the tectal location */
	  int location = *(connections + growth_cone - 1);

	  /* Get the tectal co-ordinates of the tectal position */
	  double x = *(tect_positions + location - 1);
	  double y = *(tect_positions + location + *tects - 1);

	  /* Calculate the potential value at the current location */
	  double p = p0(*(p0_fn + location - 1), *alpha, *beta, *sigma, *tau, x, y, u, v) + *(r + location - 1);

	  /* Initially set the minimums as relating to the original location */
	  double min_p = p;
	  int min_loc = location;

	  /* Loop through each neighbouring tectal location */
	  for (int n = 0; n < *tect_neighbours_size; n++)
	    {
	      /* The neighbours matrix is changed to a vector by column, when passed to this C code */
	      int neighbour = *(tect_neighbours + location + n*(*tects) - 1);
	      if (neighbour == 0)
		break;

	      /* Check if neighbouring location is active */
	      int status = *(tect_status + neighbour - 1);

	      /* if active */
	      if (status == 1)
		{
		  /* Get x, y and p for each neighbour */
		  double neighb_x = *(tect_positions + neighbour - 1);
		  double neighb_y = *(tect_positions + neighbour + *tects - 1);
		  double neighb_p = p0(*(p0_fn + neighbour - 1), *alpha, *beta, *sigma, *tau, neighb_x, neighb_y, u, v) + *(r + neighbour - 1);

		  /* If a new value of p is lower, then update values */
		  if (neighb_p < min_p)
		    {
		      min_p = neighb_p;
		      min_loc = neighbour;
		    }
		}
	    }

	  /* Update rho and connections */
	  (*(rho + location - 1))--;
	  (*(rho + min_loc - 1))++;
	  *(connections + growth_cone - 1) = min_loc;

	  /* Update r */
	  add_scale_arrays(r, rho, *tects, 1 - *kappa, *epsilon / *total_cones);
	}

      /* Make new update sequence for next iteration if wanted*/
      if ((*update_many==1)&&(*sequence==0))
	{		
	  make_new_update_sequence(&update_seq[0], seq_length, *update_method, active_cones, *total_active);
	}
    }

  PutRNGstate();

  return;
}


void
iteratedb(int *axons,		/* Integer of retinal axons */
	double *axon_positions,	/* Matrix of axon co-ordinates (u,v) on retina + gradients*/
	int *total_cones, /* Integer of number of growth cones */				
	int *cone_axons, /* Vector of axon number for each growth cone */
	int *active_cones,	/* Vector of active growth cones */
	int *total_active,	/* Integer of active growth cones */

	int *tects,		/* Integer of tectum positions */
	double *tect_positions,	/* Matrix of tectum co-ordinates (x,y)  + gradients*/
	int *tect_neighbours_size, /* Integer of tect_neighbours matrix size */
	int *tect_neighbours,	/* Matrix of tectum neighbour indices */
	int *tect_status, /* Vector of active tectal positions 1-active 0-inactive */

	int *connections,	/* Connection matrix */
				
	double *r,		/* Inhibitory term, vector size N */
	double *epsilon,	/* Rate of growth of r */
	double *kappa,		/* Rate of decay of r */

	int *rho,		/* Local density, vector size N */

	int *p0_fn,		/* Vector of values describing p0 to use */
	double *alpha,		/* p0 function parameters */
	double *beta,
	double *sigma,
	double *tau,

	int *iterations,	/* Number of times to iterate */

	int *sequence, /* User provided update sequence. *sequence == 0 if no sequence provided */
	int *s_length, /* Length of user update sequence or length of calculated sequence */

	int *update_method,	/* Integer describing update method */
	int *update_many) /* Boolean, new update sequence each iteration */
{
  double rEA, reA, rEB, reB, tEA, teA, tEB, teB;
  double x, y, u, v, p, neighb_p, *ptr;
  /* If there are no iterations, then leave */
  if (*iterations==0)
    return;

  /* Get the R environment, random number generation state */
  GetRNGstate();

  /* Set-up update sequence variables */
  int seq_length = *s_length;
  int update_seq[seq_length];
  int nneighs = *tect_neighbours_size;

  /* SJE mods */
  double dt = 1.0 / (*total_cones);
  double k1 = 1.0 - (*kappa * dt);
  double k2 = *epsilon * dt;
  
  
  /* Construct sequence in C code */
  if (*sequence == 0)
    make_new_update_sequence(&update_seq[0], seq_length, *update_method,
			     active_cones, *total_active);
  /* Use user sequence */
  else
    for (int s = 0; s < seq_length; s++)
      update_seq[s] = *(sequence + s);

  /* Loop for each iteration */
  for (int it = 0; it < *iterations; it++) {

    /* Loop for each specified density unit in the update sequence*/
    if ( (it % 50) == 0) {
      printTime(it);
    }
    for (int unit = 0; unit < seq_length; unit++) {
      /* Get a growth cone */
      int growth_cone = update_seq[unit];

      /* Get the axon */
      int axon = *(cone_axons + growth_cone - 1);
      
      /* Get the retinal co-ordinates of the axon */
      ptr = (axon_positions + axon - 1);
      u = *ptr; ptr += *axons;
      v = *ptr; ptr += *axons;
      rEA = *ptr; ptr += *axons;
      rEB = *ptr; ptr += *axons;
      reA = *ptr; ptr += *axons;
      reB = *ptr;

      /* Get the tectal location */
      int location = *(connections + growth_cone - 1);
      
      /* Get the co-ordinates of the tectal position */
      ptr = (tect_positions + location - 1);
      x = *ptr; ptr += *tects;
      y = *ptr; ptr += *tects;
      teA = *ptr; ptr += *tects;
      teB = *ptr; ptr += *tects;
      tEA = *ptr; ptr += *tects;
      tEB = *ptr;
      
      /* Calculate the potential value at the current location */
      /* double p = p0(*(p0_fn + location - 1), *alpha, *beta, *sigma, *tau, x, y, u, v) + *(r + location - 1); */
      p = p0db(*(p0_fn + location - 1),
	       rEA, reA, rEB, reB,
	       tEA, teA, tEB, teB) + *(r + location - 1);
      
      /* Initially set the minimums as relating to the original location */
      double min_p = p;
      int min_loc = location;
      
      /* Loop through each neighbouring tectal location */
      for (int n = 0; n < nneighs; n++) {
	/* The neighbours matrix is changed to a vector by column,
	   when passed to this C code */
	int neighbour = *(tect_neighbours + location + n*(*tects) - 1);
	if (neighbour == 0)
	  break;

	/* Check if neighbouring location is active */
	int status = *(tect_status + neighbour - 1);
	
	/* if active */
	if (status == 1) {
	  /* Get x, y and p for each neighbour */
	  ptr = (tect_positions + neighbour - 1);
	  x = *ptr; ptr += *tects;
	  y = *ptr; ptr += *tects;
	  teA = *ptr; ptr += *tects;
	  teB = *ptr; ptr += *tects;
	  tEA = *ptr; ptr += *tects;
	  tEB = *ptr;
	  
	  neighb_p = p0db(*(p0_fn + neighbour - 1),
			  rEA, reA, rEB, reB,
			  tEA, teA, tEB, teB) + *(r + neighbour - 1);
	  
	  /* neighb_p = p0(*(p0_fn + neighbour - 1), *alpha, *beta,
	   *sigma, *tau, neighb_x, neighb_y, u, v) + *(r + neighbour - 1); */
	  /* If a new value of p is lower, then update values */
	  if (neighb_p < min_p) {
	    min_p = neighb_p;
	    min_loc = neighbour;
	  }
	}
      }

      /* Update rho and connections */
      (*(rho + location - 1))--;
      (*(rho + min_loc - 1))++;
      *(connections + growth_cone - 1) = min_loc;
      
      /* Update r */
      add_scale_arrays2(r, rho, *tects, k1, k2);
    }
    
    /* Make new update sequence for next iteration if wanted*/
    if ((*update_many==1)&&(*sequence==0)) {
      make_new_update_sequence(&update_seq[0], seq_length,
			       *update_method, active_cones, *total_active);
    }
    /* End of one epoch. */
  }
  
  PutRNGstate();

  return;
}

/* Adds arrays which have been multiplied by scalars */
/* Replaces array1 with new values */
void add_scale_arrays(double *array1, int *array2,
		      int array_length, double scalar1, double scalar2)
{
  for (int i = 0; i < array_length; i++)
    {
      *(array1 + i) = (*(array1 + i) * scalar1) + (*(array2 + i) * scalar2);
    }
  return;
}


void add_scale_arrays2(double *array1, int *array2,
		       int array_length, double scalar1, double scalar2) {
  /* Faster implementation of add_scale_arrays */
  int i;
  for (i = array_length; i > 0; i--) {
    *array1 = (*array1 * scalar1) + (*array2 * scalar2);
    array1++; array2++;
  }
}

/* Calculate the appropriate p0 */
double p0(int p0_fn, double alpha, double beta, double sigma,
	  double tau, double x, double y, double u, double v)
{
  double result;

  double random1 = unif_rand();
  double random2 = unif_rand();
  double random3 = unif_rand();
  double random4 = unif_rand();
	
  switch(p0_fn)
    {
      /* p.gierer2, Gierer's reciprocal inhibition model for p */
    case 2:;
      result = 	(exp(-alpha*u) / exp(-alpha*x)) + (exp(-alpha*x) / exp(-alpha*u)) +
	(exp(-alpha*x) * sigma * random1) +
	(exp(-beta*v) / exp(-beta*y)) + (exp(-beta*y) / exp(-beta*v)) +
	(exp(-beta*y) * sigma * random2);
      break;
      /* p.gierer3, Gierer's single gradient model for p */
    case 3:;
      result =	(exp(-2*alpha*u) + random1*sigma) * (exp(alpha*x) + random2*tau) +
	(exp(-2*beta*v) + random3*sigma) * (exp(beta*y) + random4*tau);
      break;
    }
  return result;
}

double p0db(int p0_fn,
	    double rEA, double reA, double rEB, double reB,
	    double tEA, double teA, double tEB, double teB)
{
  double result;
  double invert = -1.0;
  switch(p0_fn) {
    case 1:			/* provide no gradient signal. */
      result = 1.0;
      break;
      /* p.gierer2, Gierer's reciprocal inhibition model for p */
    case 2:
      result =     (rEA*teA) + (reA*tEA) +  /* forward + reverse A system*/
	(invert * ((rEB*teB) + (reB*tEB)));  /* forward + reverse B system*/
      break;
    case 3:
      result =     (rEA*teA) + (reA*tEA) +  /* forward + reverse A system*/
	fabs(rEB-teB) + fabs(reB-tEB);  /* forward + reverse B system*/
      break;
    case 4:
      result =     (rEA*teA) - (rEB*teB); /* forward system only, both type I */
      break;
    default:
      Rprintf("Error: no such switch in p0db (%d).\n", p0_fn);
    }
  return result;
}

/* Create a new update sequence, using every active growth cone */
void make_new_update_sequence(int *update_seq, int seq_length,
			      int update_method, int *active_cones, int total_active)
{
  switch(update_method)
    {
      /* sample with replacement */
    case 1:
      for (int i = 0; i < seq_length; i++)
	{
	  double d_random = floor(unif_rand() * total_active);	/* Random from 0 to total_active - 1 */
	  int i_random = d_random;	/* Cast double as int */
	  *(update_seq + i) =	*(active_cones + i_random);
	}
      break;

      /* sample without replacement */
      /* Fisher-Yates Shuffle */
      /* http://en.wikipedia.org/wiki/Fisher-Yates_shuffle */
    case 2:
      *update_seq = *active_cones;
      for (int i = 1; i < seq_length; i++)
	{
	  double d_j = floor(unif_rand() * (i+1));
	  int j = d_j;
	  *(update_seq + i) =	*(update_seq + j);
	  *(update_seq + j) = *(active_cones + i);
	}
      break;
    }
  return;
}


void printTime(int it) {
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Rprintf ( "C Epoch %d %s", it, asctime (timeinfo) );
}
