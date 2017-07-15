/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_linear.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"
#include <time.h>

#include "memory.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

SolveLinear::SolveLinear(SPPARKS *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg != 1) error->all("Illegal solve command");

  random = new RandomPark(ranmaster->uniform());
  prob = NULL;

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
  srand( time(NULL) );
  prob_intval = NULL;
  prob_addr = NULL;
  sum_intval = NULL;
/*********************************/
}

/* ---------------------------------------------------------------------- */

SolveLinear::~SolveLinear()
{
  delete [] prob;
  delete random;
  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
  memory->destroy_1d_T_array(prob_intval, 0);
//  memory->destroy_1d_T_array(prob_addr, 0);
  memory->sfree(prob_addr);
  memory->destroy_2d_double_array(sum_intval);
/*********************************/
}

/* ---------------------------------------------------------------------- */

SolveLinear *SolveLinear::clone()
{
  int narg = 1;
  char *arg[1];
  arg[0] = style;

  SolveLinear *ptr = new SolveLinear(spk,narg,arg);

  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveLinear::init(int n, double *propensity)
{
  delete [] prob;
  nevents = n;
  prob = new double[n];

  sum = 0.0;
  num_active = 0;


  for (int i = 0; i < n; i++) {
    if (propensity[i] > 0.0) num_active++;
    prob[i] = propensity[i];
    sum += propensity[i];
  }
//printf("solve->init():line91:num_active=%d\n",num_active);
  
}

/* ---------------------------------------------------------------------- */
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
void SolveLinear::init_intval(int n, IntervalItem *propensity_intval)
{
	nreactions = n;
	
	// allocate memory of interval
  memory->destroy_1d_T_array(prob_intval, 0);
  memory->create_1d_T_array(prob_intval, 0, n, "chemistry:prob_intval");
	// allocate memory of interval address
//  memory->destroy_1d_T_array(prob_addr, 0);
//  memory->create_1d_T_array(prob_addr, 0, n, "chemistry:prob_addr");
  prob_addr=(IntervalItem**)memory->smalloc((n+1)*sizeof(prob_intval),
				"chemistry:prob_addr");
	
  memory->destroy_2d_double_array(sum_intval);
  sum_intval=memory->create_2d_double_array(n+1, 2, "chemistry:sum_intval");

  for (int i = 0; i < n; i++) {
    prob_intval[i].id = i;
	if(propensity_intval[i].inf <= propensity_intval[i].sup) {
	    prob_intval[i].inf = propensity_intval[i].inf;
	    prob_intval[i].sup = propensity_intval[i].sup;
	}
	else {
	    prob_intval[i].inf = propensity_intval[i].sup;
	    prob_intval[i].sup = propensity_intval[i].inf;
	}
    prob_intval[i].wid = prob_intval[i].sup-prob_intval[i].inf;
	
//	prob[i] = (prob_intval[i].inf+prob_intval[i].sup)/2.0;	//take the average approach to find the single-valued prob.
	
	prob_addr[i] = &prob_intval[i];
  }
  
	update_sum_intval();
/*	
for(int i = 0; i < n; i++){
printf("sum_intval[%d]=[%f, %f] ", i,sum_intval[i][L_BOUND],sum_intval[i][U_BOUND]);  
printf("prob_intval[%d]=[%f, %f] ", i,prob_intval[i].inf,prob_intval[i].sup);
printf("prob_addr[%d]=[%f, %f] \n", i,prob_addr[i]->inf,prob_addr[i]->sup);
	}  
*/
}
	


void SolveLinear::update_sum_intval()
{
    // sort interval valued prob_intval based on their widths
	qsort(prob_addr, nreactions, sizeof *prob_addr, compare_interval_by_wid);

	// calculate the sum of interval probabilities out of the sorted 
	sum_intval[0][L_BOUND]=prob_addr[0]->inf;
	sum_intval[0][U_BOUND]=prob_addr[0]->sup;
	for (int i = 1; i < nreactions; i++) {	
			sum_intval[i][L_BOUND]=sum_intval[i-1][U_BOUND]+prob_addr[i]->inf;
			sum_intval[i][U_BOUND]=sum_intval[i-1][L_BOUND]+prob_addr[i]->sup;
	}
	
	// introducing a Null Event so that the summation is up to a precise number
	// here, we simply choose the upper bound of the previous sum_intval
	sum_intval[nreactions][L_BOUND]=sum_intval[nreactions-1][U_BOUND];
	sum_intval[nreactions][U_BOUND]=sum_intval[nreactions-1][U_BOUND];
}


void SolveLinear::update_intval(int n, int *indices, IntervalItem *propensity_intval)
{
  int m;
  for (int i = 0; i < n; i++) {
    m = indices[i];
	prob_intval[m].inf=propensity_intval[m].inf;
	prob_intval[m].sup=propensity_intval[m].sup;
	prob[m] = (prob_intval[m].inf+prob_intval[m].sup)/2.0;	//take the average approach to find the single-valued prob.
  }
  update_sum_intval();
}

/****
Identify the random set of events that will be fired, return the (int)  number of events
****/
int SolveLinear::event_intval(int *ireactionList, IntervalItem *pdt)
{
	update_sum_intval();
	if (num_active == 0) 
		return 0;
	
/***	double fraction = sum_intval[nreactions][U_BOUND] * random->uniform();	***/
	
/***/
	double fraction = sum_intval[nreactions][U_BOUND] 
						* ((double)rand())/RAND_MAX;
/***/
	int index_L = 0, index_U = nreactions;
	int i, numEvents;
//printf("\n fraction=%f ",fraction);

//for(i=0;i<nreactions+1;i++) 
//printf(" sum_intval[%d]=[%f, %f] ", i,sum_intval[i][L_BOUND],sum_intval[i][U_BOUND]); 		
	
	// find the lowest index of possible events
	for(i=0; i<nreactions; i++) {
		if (fraction <= sum_intval[i][U_BOUND]) { index_L = i; break;}
	}

	// find the highest index of possible events
	for(i=nreactions; i>0; i--) {
		if (fraction >= sum_intval[i-1][L_BOUND]) { index_U = i; break;}
	}
	
	numEvents = index_U-index_L+1;

	pdt->inf=0.0; pdt->sup=0.0;
/***	double u = random->uniform();	***/
/***/	double u = ((double)rand())/RAND_MAX;	/***/ //A more random number generator
	IntervalItem sum; 
	sum.inf=0.0; sum.sup=0.0;
	for(i=0; i<nreactions; i++) {
		sum.inf += prob_intval[i].inf; 
		sum.sup += prob_intval[i].sup;
	}
	
//printf("\n numberEvents=%d: %d~%d \n",numEvents,index_L,index_U); 		
	for(i=0; i<numEvents; i++) {
		if(index_L+i != nreactions) // make sure it is not a null event - the last event is a null event
		{
			ireactionList[i]= prob_addr[index_L+i]->id;
		//	*pdt += -1.0/sum_intval[nreactions][U_BOUND] * log(u); //sum of time for "numEvents" number of events fired at the same time
			if(pdt->inf<=0.0) pdt->inf = -1.0/sum.sup * log(u); //sum of time for "numEvents" number of events fired at the same time
			pdt->sup += 1.0/sum.inf;
			sum.inf -= prob_addr[index_L+i]->inf;
		//	sum.sup -= prob_addr[index_L+i]->sup;
		}
		else
			ireactionList[i]= nreactions; // the last event is a null event
	}	
	pdt->sup = pdt->sup * (-log(u));
	//time increment
	//*********************************************
	// the minimum time of any event occurs: this returns only one event time
//	*pdt = -1.0/sum_intval[nreactions][U_BOUND] * log(random->uniform());	
	
	//**********************************************
	// the average time of the random set of fired events
/*	double kL=0.0, kU=0.0; 
	int m=0;
	for(i=0; i<numEvents; i++) {
		if(index_L+i != nreactions) // make sure it is not a null event - the last event is a null event
		{
			kL += -1.0/prob_addr[index_L+i]->sup;
			kU += -1.0/prob_addr[index_L+i]->inf; 
			m++;
		}
		else
			ireactionList[i]= nreactions; // the last event is a null event
	}	
	*pdt = (kL+kU)/(2*m) * log(random->uniform());	
*/
	//*********************************************
//printf(" dt=%f  ",*pdt);
	  
	return numEvents;
}

/*********************************/
/* ---------------------------------------------------------------------- */

void SolveLinear::update(int n, int *indices, double *propensity)
{
  int m;
  for (int i = 0; i < n; i++) {
    m = indices[i];
    if (prob[m] > 0.0) num_active--;
    if (propensity[m] > 0.0) num_active++;
    sum -= prob[m];
    prob[m] = propensity[m];
    sum += propensity[m];
  }
}
/* ---------------------------------------------------------------------- */

void SolveLinear::update(int n, double *propensity)
{
  if (prob[n] > 0.0) num_active--;
  if (propensity[n] > 0.0) num_active++;
  sum -= prob[n];
  prob[n] = propensity[n];
  sum += propensity[n];
}
/* ---------------------------------------------------------------------- */


void SolveLinear::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}
/* ---------------------------------------------------------------------- */

int SolveLinear::event(double *pdt)
{
  int m;

  if (num_active == 0) {
    sum = 0.0;
    return -1;
  }

/********** commented by Yan Wang and replaced by the following ***  
  double fraction = sum * random->uniform(); 
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
 double fraction = sum * ((double)rand())/RAND_MAX; //A more random number generator
/****************************/
  double partial = 0.0;

  for (m = 0; m < nevents; m++) {
    partial += prob[m];
    if (partial > fraction) break;
  }

/********** commented by Yan Wang and replaced by the following ***  
  *pdt = -1.0/sum * log(random->uniform());
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
  *pdt = -1.0/sum * log(((double)rand())/RAND_MAX);
/****************************/

 return m;
}

