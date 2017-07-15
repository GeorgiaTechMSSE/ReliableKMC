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

#include "string.h"
#include "solve.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Solve::Solve(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);
  sum = 0.0;
}

/* ---------------------------------------------------------------------- */

Solve::~Solve()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */

double Solve::get_total_propensity()
{
  return sum;
}

/* ---------------------------------------------------------------------- */

int Solve::get_num_active()
{
  return num_active;
}



/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/

int compare_interval_by_id(const void *a, const void *b)  
{
   const struct IntervalItem *const *x = (struct IntervalItem **)a;
   const struct IntervalItem *const *y = (struct IntervalItem **)b;
   if ( (*x)->id > (*y)->id ) return 1;
   if ( (*x)->id < (*y)->id ) return -1;
   return 0;
}

int compare_interval_by_wid(const void *a, const void *b)  
{
   const struct IntervalItem *const *x = (struct IntervalItem **)a;
   const struct IntervalItem *const *y = (struct IntervalItem **)b;
   if ( (*x)->wid > (*y)->wid ) return 1;
   if ( (*x)->wid < (*y)->wid ) return -1;
   return 0;
} 

int compare_interval_by_inf(const void *a, const void *b)  
{
   const struct IntervalItem *const *x = (struct IntervalItem **)a;
   const struct IntervalItem *const *y = (struct IntervalItem **)b;
   if ( (*x)->inf > (*y)->inf ) return 1;
   if ( (*x)->inf < (*y)->inf ) return -1;
   return 0;
}

int compare_interval_by_sup(const void *a, const void *b)  
{
   const struct IntervalItem *const *x = (struct IntervalItem **)a;
   const struct IntervalItem *const *y = (struct IntervalItem **)b;
   if ( (*x)->sup > (*y)->sup ) return 1;
   if ( (*x)->sup < (*y)->sup ) return -1;
   return 0;
}  
/*********************************/