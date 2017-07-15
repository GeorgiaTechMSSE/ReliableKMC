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

#ifdef SOLVE_CLASS
SolveStyle(linear,SolveLinear)

#else

#ifndef SPK_SOLVE_LINEAR_H
#define SPK_SOLVE_LINEAR_H

#include "solve.h"

namespace SPPARKS_NS {

class SolveLinear : public Solve {
 public:
  SolveLinear(class SPPARKS *, int, char **);
  ~SolveLinear();
  SolveLinear *clone();

  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
	void init_intval(int, IntervalItem *);
	void update_intval(int, int *, IntervalItem *);
	void update_sum_intval();
	int event_intval(int *, IntervalItem *);
/*********************************/

 private:
  class RandomPark *random;
  int nevents;
  double *prob;
  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
	int nreactions;					//the total number of reactions/events
	IntervalItem *prob_intval;		//the interval propensities
	IntervalItem **prob_addr;		//keep track of the sorted interval propensities based on the widths
	double **sum_intval;			//keep sums of sorted interval propensities
/*********************************/

};

}

#endif
#endif
