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
#include "string.h"
#include "stdlib.h"
#include "app_diffusion.h"
#include "solve.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED,TOP};
enum{NO_ENERGY,LINEAR,NONLINEAR};
enum{DEPOSITION,NNHOP,SCHWOEBEL};
enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};  // from app_lattice.cpp

#define DELTAEVENT 100000

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppDiffusion::AppDiffusion(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  // these can be changed by model choice, see below

  ninteger = 1;
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  numrandom = 1;

  // parse arguments

  if (narg < 3) error->all("Illegal app_style command");
  if (strcmp(arg[1],"off") == 0) engstyle = NO_ENERGY;
  else if (strcmp(arg[1],"linear") == 0) engstyle = LINEAR;
  else if (strcmp(arg[1],"nonlinear") == 0) engstyle = NONLINEAR;
  else error->all("Illegal app_style command");

  if (strcmp(arg[2],"hop") == 0) {
    if (narg != 3) error->all("Illegal app_style command");
    hopstyle = NNHOP;
  } else if (strcmp(arg[2],"schwoebel") == 0) {
    if (narg != 5) error->all("Illegal app_style command");
    hopstyle = SCHWOEBEL;
    nsmax = atoi(arg[3]);
    nsmin = atoi(arg[4]);
  } else error->all("Illegal app_style command");

  // increment delpropensity by 1 for nonlinear energy
  // increment delpropensity and delevent by 1 for Schwoebel hops
  // change allow_rejection to 1 for linear energy and non-Schwoebel hops

  if (engstyle == NONLINEAR) delpropensity++;
  if (hopstyle == SCHWOEBEL) delpropensity++;
  if (hopstyle == SCHWOEBEL) delevent++;
  if (engstyle == LINEAR && hopstyle == NNHOP) allow_rejection = 1;

  create_arrays();

  esites = psites = NULL;
  echeck = pcheck = NULL;
  maxevent = 0;
  events = NULL;
  firstevent = NULL;

  hbarrier = sbarrier = NULL;
  /****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	hbarrier_intval = sbarrier_intval = NULL;
/*****************************/

  ecoord = NULL;
  hopsite = NULL;
  marklist = NULL;
  mark = NULL;

  allocated = 0;

  // default settings for app-specific commands

  depflag = 0;
  barrierflag = 0;

  // statistics

  ndeposit = ndeposit_failed = 0;
  nfirst = nsecond = 0;
}

/* ---------------------------------------------------------------------- */

AppDiffusion::~AppDiffusion()
{
  delete [] esites;
  delete [] psites;
  delete [] echeck;
  delete [] pcheck;
  memory->sfree(events);
  memory->sfree(firstevent);

  memory->destroy_2d_double_array(hbarrier);
  memory->destroy_2d_double_array(sbarrier);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	memory->destroy_2d_T_array(hbarrier_intval);
	memory->destroy_2d_T_array(sbarrier_intval);

/*****************************/	
  
  delete [] ecoord;

  delete [] hopsite;
  delete [] marklist;
  memory->sfree(mark);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppDiffusion::input_app(char *command, int narg, char **arg)
{
    double **barrier;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	IntervalItem **barrier_intval;
/*****************************/	

  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

  if (strcmp(command,"ecoord") == 0) {
    if (engstyle != NONLINEAR)
      error->all("Can only use ecoord command with "
		 "app_style diffusion nonlinear");
    if (narg != 2) error->all("Illegal ecoord command");

    int lo,hi;
    bounds(arg[0],maxneigh,lo,hi);
    double value = atof(arg[1]);

    for (int i = lo; i <= hi; i++) ecoord[i] = value;

  } 
  else if (strcmp(command,"deposition") == 0) 
  {
    if (narg < 1) error->all("Illegal deposition command");
    if (strcmp(arg[0],"off") == 0) {
      if (narg != 1 ) error->all("Illegal deposition command");
      depflag = 0;
      return;
    }
//printf("maxneigh=%d\n",maxneigh);

/******************************* commented out by Yan Wang (June 2010) and replaced by the following
    if (narg != 7) error->all("Illegal deposition command");
    depflag = 1;
    deprate = atof(arg[0]);
    dir[0] = atof(arg[1]);
    dir[1] = atof(arg[2]);
    dir[2] = atof(arg[3]);
    d0 = atof(arg[4]);
    coordlo = atoi(arg[5]);
    coordhi = atoi(arg[6]);
    if (deprate < 0.0) error->all("Illegal deposition command: deposite rate should not be negative");
    if (domain->dimension == 2 && (dir[1] >= 0.0 || dir[2] != 0.0))
      error->all("Illegal deposition command: deposite direction does not match dimension");
    if (domain->dimension == 3 && dir[2] >= 0.0)
      error->all("Illegal deposition command: deposite direction does not match dimension");
    if (d0 < 0.0) error->all("Illegal deposition command: capture distance should not be negative");
    if (coordlo < 0 || coordhi > maxneigh || coordlo > coordhi)
      error->all("Illegal deposition command");
******************************************/

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
    if (narg != 8) error->all("Illegal deposition command");
    depflag = 1;
    deprate_intval.inf = atof(arg[0]);
	deprate_intval.sup = atof(arg[1]);
	deprate = (deprate_intval.inf+deprate_intval.sup)/2.0;
    dir[0] = atof(arg[2]);
    dir[1] = atof(arg[3]);
    dir[2] = atof(arg[4]);
    d0 = atof(arg[5]);
    coordlo = atoi(arg[6]);
    coordhi = atoi(arg[7]);
    if (deprate_intval.inf < 0.0 || deprate_intval.sup < 0.0)
		error->all("Illegal deposition command: deposite rate should not be negative");
    if (domain->dimension == 2 && (dir[1] >= 0.0 || dir[2] != 0.0))
      error->all("Illegal deposition command: deposite direction does not match dimension");
    if (domain->dimension == 3 && dir[2] >= 0.0)
      error->all("Illegal deposition command: deposite direction does not match dimension");
    if (d0 < 0.0) error->all("Illegal deposition command: capture distance should not be negative");
    if (coordlo < 0 || coordhi > maxneigh || coordlo > coordhi)
      error->all("Illegal deposition command");
/********************************/	  
	  
    double len = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0] /= len;
    dir[1] /= len;
    dir[2] /= len;

  } 
  else if (strcmp(command,"barrier") == 0) 
  {
    if (narg < 1) error->all("Illegal barrier command");
    barrierflag = 1;

    if (strcmp(arg[0],"none") == 0) {
      if (narg != 1) error->all("Illegal barrier command");
      barrierflag = 0;
      return;
    } else if (strcmp(arg[0],"hop") == 0) {
      barrier = hbarrier;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
		barrier_intval = hbarrier_intval;
/*****************************/	

    } else if (strcmp(arg[0],"schwoebel") == 0) {
      barrier = sbarrier;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
		barrier_intval = sbarrier_intval;
/*****************************/	
	  
    } else error->all("Illegal barrier command");
	
    if (barrier == sbarrier && hopstyle != SCHWOEBEL)
      error->all("Cannot define Schwoebel barrier without Schwoebel model");


/**********************   commented out by Yan Wang, June 2010, replaced by the following	  
    if (narg < 2 || narg > 4) error->all("Illegal barrier command");
    if (narg == 2) {
      double q = atof(arg[1]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
		for (j = 0; j <= maxneigh; j++)
			barrier[i][j] = q;

    } else if (narg == 3) {
      int delta = atoi(arg[1]);
      double q = atof(arg[2]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
		for (j = 0; j <= maxneigh; j++)
			if (j-i == delta) barrier[i][j] = q;

    } else {
      int ilo,ihi,jlo,jhi;
      bounds(arg[1],maxneigh,ilo,ihi);
      bounds(arg[2],maxneigh,jlo,jhi);
      double q = atof(arg[3]);

      for (int i = ilo; i <= ihi; i++)
		for (int j = jlo; j <= jhi; j++)
			barrier[i][j] = q;
    }
********************/


/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  	  
    if (narg < 3 || narg > 5) error->all("Illegal barrier command");
    if (narg == 3) {
      double qL = atof(arg[1]);
	  double qU = atof(arg[2]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
		for (j = 0; j <= maxneigh; j++) {
			barrier[i][j] = (qL+qU)/2.0; 
			barrier_intval[i][j].inf = qL; 
			barrier_intval[i][j].sup = qU;
		}
    } else if (narg == 4) {
      int delta = atoi(arg[1]);
      double qL = atof(arg[2]);
      double qU = atof(arg[3]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
		for (j = 0; j <= maxneigh; j++)
			if (j-i == delta) {
				barrier[i][j] = (qL+qU)/2.0;
				barrier_intval[i][j].inf = qL;
				barrier_intval[i][j].sup = qU;
			}

    } else {
      int ilo,ihi,jlo,jhi;
      bounds(arg[1],maxneigh,ilo,ihi);
      bounds(arg[2],maxneigh,jlo,jhi);
      double qL = atof(arg[3]);
      double qU = atof(arg[4]);
	
      for (int i = ilo; i <= ihi; i++)
		for (int j = jlo; j <= jhi; j++) {
			barrier[i][j] = (qL+qU)/2.0;
			barrier_intval[i][j].inf = qL;
			barrier_intval[i][j].sup = qU;
		}
    }
/******************************/

	
  } else error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppDiffusion::grow_app()
{
  lattice = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppDiffusion::init_app()
{
  if (depflag && nprocs > 1)
    error->all("Cannot perform deposition in parallel");
  if (depflag && nsector > 1)
    error->all("Cannot perform deposition with multiple sectors");

  if (!allocated) allocate_data();
  allocated = 1;

  dimension = domain->dimension;

  // sweeping timestep

  dt_sweep = 1.0/maxneigh;

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (lattice[i] < VACANT || lattice[i] > TOP) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all("One or more sites have invalid values");

}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppDiffusion::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = pcheck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;
  
 
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppDiffusion::site_energy(int i)
{
  // energy only non-zero for OCCUPIED sites when energy included in model

  if (lattice[i] != OCCUPIED || engstyle == NO_ENERGY) return 0.0;

  // energy is a non-linear function of coordination number
  // calculate from user-specified tabulated values

  if (engstyle == NONLINEAR) {
    int n = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) n++;
    return ecoord[n];
  }

  // energy is a linear function of coordination number, just count bonds

  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   null bin extends to size maxneigh
------------------------------------------------------------------------- */

void AppDiffusion::site_event_rejection(int i, RandomPark *random)
{
  double einitial,edelta;

  // OCCUPIED site exchanges with random neighbor if VACANT

  if (lattice[i] != OCCUPIED) return;
  int iran = (int) (maxneigh*random->uniform());
  if (iran > maxneigh) iran = maxneigh-1;
  int j = neighbor[i][iran];
  if (lattice[j] != VACANT) return;

  // accept or reject via energy and barrier model
  // factor of 2 in edelta accounts for energy change of neighbors of I,J

  int hop = 0;
  if (engstyle != NO_ENERGY) einitial = site_energy(i);

  lattice[i] = VACANT;
  lattice[j] = OCCUPIED;

  if (engstyle == NO_ENERGY) {
    if (!barrierflag) hop = 1;
    else if (temperature > 0.0) {
      if (random->uniform() < exp(-hbarrier[ncoord(i)-1][ncoord(j)]*t_inverse))
	  hop = 1;
    }

  } else {
    edelta = site_energy(j) - einitial;

    if (!barrierflag) {
      if (edelta <= 0.0) hop = 1;
      else if (temperature > 0.0) {
	if (random->uniform() < exp(-2.0*edelta*t_inverse)) hop = 1;
      }
    } else if (temperature > 0.0) {
      if (edelta <= 0.0) {
	if (random->uniform() < 
	    exp(-hbarrier[ncoord(i)-1][ncoord(j)]*t_inverse)) hop = 1;
      } else {
	if (random->uniform() < 
	    exp((-2.0*edelta-hbarrier[ncoord(i)-1][ncoord(j)]) * t_inverse))
	  hop = 1;
      }
    }
  }

  if (hop) {
    naccept++;
    nfirst++;
  } else {
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusion::site_propensity(int i)
{
  if (engstyle == NO_ENERGY) return site_propensity_no_energy(i);
  else if (engstyle == LINEAR) return site_propensity_linear(i);
  else return site_propensity_nonlinear(i);
}
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
IntervalItem AppDiffusion::site_propensity_intval(int i)
{
	return site_propensity_intval_linear(i);
}
/*****************************/ 	

/* ---------------------------------------------------------------------- */

double AppDiffusion::site_propensity_no_energy(int i)
{
  int j,ihop,nhop1,nhop2,eflag;
  double einitial,edelta,probone,proball;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engstyle, compute edelta between initial and final state
  // factor of 2 in edelta accounts for energy change of neighbors of I,J
  // if barrierflag, compute coordination of sites I,J
  // for 1st neigh hop, set delta = 1 to remove hopping atom itself
  // propensity is function of edelta, barrier, up/down hill, temperature

  clear_events(i);

  if (lattice[i] != OCCUPIED) return 0.0;

  // nhop1 = 1st neigh hops, nhop2 = 2nd neigh hops
  // hopsite = all possible hop sites

  nhop1 = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) hopsite[nhop1++] = neighbor[i][j];
  nhop2 = nhop1;
  if (hopstyle == SCHWOEBEL) 
    nhop2 += schwoebel_enumerate(i,&hopsite[nhop1]);

  // loop over all possible hops

  proball = 0.0;
  double **barrier = hbarrier;
  int delta = 1;

  for (ihop = 0; ihop < nhop2; ihop++) {
    j = hopsite[ihop];
    if (ihop == nhop1) {
      barrier = sbarrier;
      delta = 0;
    }

    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
    probone = 0.0;

    if (!barrierflag) probone = 1.0;
    else if (temperature > 0.0)
      probone = exp(-barrier[ncoord(i)-delta][ncoord(j)]*t_inverse);
    
    if (probone > 0.0) {
      eflag = (ihop < nhop1) ? NNHOP : SCHWOEBEL;
      add_event(i,j,probone,eflag);
      proball += probone;
    }
    
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate,DEPOSITION);
    proball += deprate;
  }

  return proball;
}

/* ---------------------------------------------------------------------- */

double AppDiffusion::site_propensity_linear(int i)
{
  int j,ihop,nhop1,nhop2,eflag;
  double einitial,edelta,probone,proball;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engstyle, compute edelta between initial and final state
  // factor of 2 in edelta accounts for energy change of neighbors of I,J
  // if barrierflag, compute coordination of sites I,J
  // for 1st neigh hop, set delta = 1 to remove hopping atom itself
  // propensity is function of edelta, barrier, up/down hill, temperature

  clear_events(i);

  if (lattice[i] != OCCUPIED) return 0.0;

  // nhop1 = 1st neigh hops, nhop2 = 2nd neigh hops
  // hopsite = all possible hop sites

  nhop1 = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) hopsite[nhop1++] = neighbor[i][j];
  nhop2 = nhop1;
  if (hopstyle == SCHWOEBEL) 
    nhop2 += schwoebel_enumerate(i,&hopsite[nhop1]);

  // loop over all possible hops

  einitial = site_energy(i);

  proball = 0.0;
  double **barrier = hbarrier;
  int delta = 1;

  for (ihop = 0; ihop < nhop2; ihop++) {
    j = hopsite[ihop];
    if (ihop == nhop1) {
      barrier = sbarrier;
      delta = 0;
    }

    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
    probone = 0.0;

    edelta = site_energy(j) - einitial;
      
    if (!barrierflag) {
      if (edelta <= 0.0) probone = 1.0;
      else if (temperature > 0.0) 
	probone = exp(-2.0*edelta*t_inverse);
    } 
	else if (temperature > 0.0) {
      if (edelta <= 0.0)
	probone = exp(-barrier[ncoord(i)-delta][ncoord(j)]*t_inverse);
      else
	probone = exp((-2.0*edelta-barrier[ncoord(i)-delta][ncoord(j)]) * 
		      t_inverse);
    }
    
    if (probone > 0.0) {
      eflag = (ihop < nhop1) ? NNHOP : SCHWOEBEL;
      add_event(i,j,probone,eflag);
      proball += probone;

    }
    
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate,DEPOSITION);
    proball += deprate;
  }

  return proball;
}


/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
IntervalItem AppDiffusion::site_propensity_intval_linear(int i)
{
  int j,ihop,nhop1,nhop2,eflag;
  double einitial,edelta; 
  double tmp;
  IntervalItem probone,proball;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engstyle, compute edelta between initial and final state
  // factor of 2 in edelta accounts for energy change of neighbors of I,J
  // if barrierflag, compute coordination of sites I,J
  // for 1st neigh hop, set delta = 1 to remove hopping atom itself
  // propensity is function of edelta, barrier, up/down hill, temperature

  clear_events(i);

  if (lattice[i] != OCCUPIED) {
    proball.inf = 0.0;
	proball.sup = 0.0;
	proball.wid = 0.0;
	return proball;
  }

  // nhop1 = 1st neigh hops, nhop2 = 2nd neigh hops
  // hopsite = all possible hop sites

  nhop1 = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) hopsite[nhop1++] = neighbor[i][j];
  nhop2 = nhop1;
  if (hopstyle == SCHWOEBEL) 
    nhop2 += schwoebel_enumerate(i,&hopsite[nhop1]);

  // loop over all possible hops

  einitial = site_energy(i);

    proball.inf = 0.0;
	proball.sup = 0.0;
	proball.wid = 0.0;
  double **barrier = hbarrier;
  IntervalItem **barrier_intval = hbarrier_intval;
  
  int delta = 1;

  for (ihop = 0; ihop < nhop2; ihop++) 
  {
    j = hopsite[ihop];
    if (ihop == nhop1) {
      barrier = sbarrier;
	  barrier_intval = sbarrier_intval;
      delta = 0;
    }

    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
    
    probone.inf = 0.0;
	probone.sup = 0.0;
	probone.wid = 0.0;
    edelta = site_energy(j) - einitial;
      
    if (!barrierflag) {
      if (edelta <= 0.0) {
		 probone.inf = 1.0;
		 probone.sup = 1.0;
		 probone.wid = 0.0;  
	  }
      else if (temperature > 0.0) {
		 double t = exp(-2.0*edelta*t_inverse);
		 probone.inf = t;
		 probone.sup = t;
		 probone.wid = 0.0;
	  }
    } 
	else if (temperature > 0.0) {
      if (edelta <= 0.0) {
		 probone.inf = exp(-barrier_intval[ncoord(i)-delta][ncoord(j)].sup*t_inverse);
		 probone.sup = exp(-barrier_intval[ncoord(i)-delta][ncoord(j)].inf*t_inverse);
		 probone.wid = probone.sup - probone.inf;
			if (probone.wid < 0){
				probone.wid = -probone.wid;
				tmp = probone.inf; probone.inf = probone.sup; probone.sup = tmp;
			}
		}
      else {
		 probone.inf = exp((-2.0*edelta-barrier_intval[ncoord(i)-delta][ncoord(j)].sup)*t_inverse);
		 probone.sup = exp((-2.0*edelta-barrier_intval[ncoord(i)-delta][ncoord(j)].inf)*t_inverse);
		 probone.wid = probone.sup - probone.inf;
			if (probone.wid < 0){
				probone.wid = -probone.wid;
				tmp = probone.inf; probone.inf = probone.sup; probone.sup = tmp;
			}
		}
    }
    
    if ( probone.inf>0.0 || probone.sup>0.0 ) {
      eflag = (ihop < nhop1) ? NNHOP : SCHWOEBEL;
      add_event(i,j,probone,(probone.inf+probone.sup)/2.0,eflag);
      proball.inf += probone.inf;
      proball.sup += probone.sup;
      proball.wid = proball.sup-proball.inf;
			if (proball.wid < 0){
				proball.wid = -proball.wid;
				tmp = proball.inf; proball.inf = proball.sup; proball.sup = tmp;
			}

    }
    
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate_intval,(deprate_intval.inf+deprate_intval.sup)/2.0,DEPOSITION);
      proball.inf += deprate_intval.inf;
      proball.sup += deprate_intval.sup;
      proball.wid = proball.sup-proball.inf;
			if (proball.wid < 0){
				proball.wid = -proball.wid;
				tmp = proball.inf; proball.inf = proball.sup; proball.sup = tmp;
			}
  }

  return proball;
}

/*****************************/ 	

/* ---------------------------------------------------------------------- */

double AppDiffusion::site_propensity_nonlinear(int i)
{
  int j,k,m,nsites,ihop,nhop1,nhop2,eflag;
  double einitial,efinal,edelta,probone,proball;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engstyle, compute edelta between initial and final state
  // since eng is nonlinear, this must include eng of neighbor sites of I,J
  // if barrierflag, compute coordination of sites I,J
  // for 1st neigh hop, set delta = 1 to remove hopping atom itself
  // propensity is function of edelta, barrier, up/down hill, temperature

  clear_events(i);

  if (lattice[i] != OCCUPIED) return 0.0;

  // nhop1 = 1st neigh hops, nhop2 = 2nd neigh hops
  // hopsite = all possible hop sites

  nhop1 = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) hopsite[nhop1++] = neighbor[i][j];
  nhop2 = nhop1;
  if (hopstyle == SCHWOEBEL) 
    nhop2 += schwoebel_enumerate(i,&hopsite[nhop1]);

  // loop over all possible hops

  proball = 0.0;
  double **barrier = hbarrier;
  int delta = 1;

  for (ihop = 0; ihop < nhop2; ihop++) {
    j = hopsite[ihop];
    if (ihop == nhop1) {
      barrier = sbarrier;
      delta = 0;
    }

    probone = 0.0;

    // einitial = i,j and their neighbors
    // use pcheck[] to avoid recomputing energy of same site

    einitial = site_energy(i) + site_energy(j);
    nsites = 0;
    psites[nsites++] = i;
    psites[nsites++] = j;
    pcheck[i] = pcheck[j] = 1;
      
    for (k = 0; k < numneigh[i]; k++) {
      m = neighbor[i][k];
      if (pcheck[m]) continue;
      einitial += site_energy(m);
      psites[nsites++] = m;
      pcheck[m] = 1;
    }
    for (k = 0; k < numneigh[j]; k++) {
      m = neighbor[j][k];
      if (pcheck[m]) continue;
      einitial += site_energy(m);
      psites[nsites++] = m;
      pcheck[m] = 1;
    }
    
    for (m = 0; m < nsites; m++) pcheck[psites[m]] = 0;
    
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
      
    // efinal = i,j and their neighbors
    // use pcheck[] to avoid recomputing energy of same site
    
    efinal = site_energy(i) + site_energy(j);
    nsites = 0;
    psites[nsites++] = i;
    psites[nsites++] = j;
    pcheck[i] = pcheck[j] = 1;
    
    for (k = 0; k < numneigh[i]; k++) {
      m = neighbor[i][k];
      if (pcheck[m]) continue;
      efinal += site_energy(m);
      psites[nsites++] = m;
      pcheck[m] = 1;
    }
    for (k = 0; k < numneigh[j]; k++) {
      m = neighbor[j][k];
      if (pcheck[m]) continue;
      efinal += site_energy(m);
      psites[nsites++] = m;
      pcheck[m] = 1;
    }
    
    for (m = 0; m < nsites; m++) pcheck[psites[m]] = 0;
    
    edelta = efinal - einitial;

    if (!barrierflag) {
      if (edelta <= 0.0) probone = 1.0;
      else if (temperature > 0.0) 
	probone = exp(-edelta*t_inverse);
    } else if (temperature > 0.0) {
      if (edelta <= 0.0)
	probone = exp(-barrier[ncoord(i)-delta][ncoord(j)]*t_inverse);
      else
	probone = exp((-edelta-barrier[ncoord(i)-delta][ncoord(j)]) * 
		      t_inverse);
    }
    
    if (probone > 0.0) {
      eflag = (ihop < nhop1) ? NNHOP : SCHWOEBEL;
      add_event(i,j,probone,eflag);
      proball += probone;
    }
    
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate,DEPOSITION);
    proball += deprate;
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusion::site_event(int i, class RandomPark *random)
{
  if (engstyle == NO_ENERGY || engstyle == LINEAR)
/***************** commented out by Yan Wang (June 2010) and replaced by the following
	return site_event_linear(i,random);
************************************/
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
	return site_event_intval_linear(i,random);
/******************************/
  else return site_event_nonlinear(i,random);
}

/* ---------------------------------------------------------------------- */

void AppDiffusion::site_event_linear(int i, class RandomPark *random)
{
  int j,k,kk,kkk,m,mm,mmm,isite;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one("Did not reach event propensity threshhold");
  }

  // deposition or hop event
  // for deposition event, find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  if (events[ievent].style == DEPOSITION) {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = OCCUPIED;
    i = j = m;
  } else {
    j = events[ievent].destination;
    if (events[ievent].style == NNHOP) nfirst++;
    else nsecond++;
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
  }

  // compute propensity changes for self and swap site and their neighs
  // 1,2 neighs for NNHOP and 1,2,3 neighs for SCHWOEBEL
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);

  esites[nsites++] = isite;
  echeck[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);

    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  if (hopstyle == NNHOP) {
    nsites += neighbor2(i,&esites[nsites]);
    nsites += neighbor2(j,&esites[nsites]);
  } else {
    nsites += neighbor3(i,&esites[nsites]);
    nsites += neighbor3(j,&esites[nsites]);
  }

  solve->update(nsites,esites,propensity);
  // sanity check on all propensity values

  /*
  printf("EVENT %d %d\n",i,j);
  for (m = 0; m < nlocal; m++) {
    if (fabs(propensity[m]-site_propensity(m)) > 1.0e-6) {
      printf("BAD PROP = %d %d %d %g %g\n",
	     id[i],id[j],id[m],propensity[m],site_propensity(m));
      error->one("BAD DONE");
    }
  }
  */

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;


}


/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
void AppDiffusion::site_event_intval_linear(int i, class RandomPark *random)
{
  int j,k,kk,kkk,m,mm,mmm,isite;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event
/********************************************************/
/*printf("\n firstevent[%d]=%d i2site=%d \n",i,firstevent[i], i2site[i]);

for (int m = 0; m < nevents; m++)
		{ printf("events[%d](styl %d, dest.%d, prop. %f, next %d,) ", 
			m, events[m].style, events[m].destination, events[m].propensity,  events[m].next); 
		}   
*/		
/******************************************************/		
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one("Did not reach event propensity threshhold");
  }

  // deposition or hop event
  // for deposition event, find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  if (events[ievent].style == DEPOSITION) {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = OCCUPIED;
    i = j = m;
  } else {
    j = events[ievent].destination;
    if (events[ievent].style == NNHOP) nfirst++;
    else nsecond++;
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
  }

  // compute propensity changes for self and swap site and their neighs
  // 1,2 neighs for NNHOP and 1,2,3 neighs for SCHWOEBEL
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
//  propensity[isite] = site_propensity(i);

	propensity_intval[isite] = site_propensity_intval(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
//    propensity[isite] = site_propensity(j);
	propensity_intval[isite] = site_propensity_intval(j);

    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  if (hopstyle == NNHOP) {
    nsites += neighbor2(i,&esites[nsites]);
    nsites += neighbor2(j,&esites[nsites]);
  } else {
    nsites += neighbor3(i,&esites[nsites]);
    nsites += neighbor3(j,&esites[nsites]);
  }

  solve->update(nsites,esites,propensity);
 
////	solve->update_intval(nsites,esites,propensity_intval);	//the update of only a few related ones is wrong because of even the number of events has changed!
	solve->init_intval(nlocal,propensity_intval);
/******************
for(int i=0;i<nsites;i++) printf("esites[%d]=%d ",i,esites[i]);
printf("\n updated \n");
for(int m=0; m<set[0].nlocal; m++) 
	printf("propensity_intval[%d]=(%d,%f,%f,%f) ",
			m, propensity_intval[m].id, propensity_intval[m].inf,
			   propensity_intval[m].sup,propensity_intval[m].wid);
*****************/
  // sanity check on all propensity values

  /*
  printf("EVENT %d %d\n",i,j);
  for (m = 0; m < nlocal; m++) {
    if (fabs(propensity[m]-site_propensity(m)) > 1.0e-6) {
      printf("BAD PROP = %d %d %d %g %g\n",
	     id[i],id[j],id[m],propensity[m],site_propensity(m));
      error->one("BAD DONE");
    }
  }
  */

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;

/****************************************/
/*printf("\n event[%d] (%d->%d) fired \n",ievent,i,j);
for (int m = 0; m < nevents; m++)
		{ printf("events[%d](styl %d, dest.%d, prop. %f, next %d,) ", 
			m, events[m].style, events[m].destination, events[m].propensity,  events[m].next); 
		}   
*/		
/***************************************/
}
/*****************************/ 	


/* ---------------------------------------------------------------------- */

void AppDiffusion::site_event_nonlinear(int i, class RandomPark *random)
{
  int j,k,kk,kkk,m,mm,mmm,isite;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one("Did not reach event propensity threshhold");
  }

  // deposition or hop event
  // for deposition event, find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  if (events[ievent].style == DEPOSITION) {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = OCCUPIED;
    i = j = m;
  } else {
    j = events[ievent].destination;
    if (events[ievent].style == NNHOP) nfirst++;
    else nsecond++;
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
  }

  // compute propensity changes for self and swap site and their neighs
  // 1,2,3 neighs for NNHOP and 1,2,3,4 neighs for SCHWOEBEL
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(i);
/*****************************/ 	
  esites[nsites++] = isite;
  echeck[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(j);
/*****************************/ 	
    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  if (hopstyle == NNHOP) {
    nsites += neighbor3(i,&esites[nsites]);
    nsites += neighbor3(j,&esites[nsites]);
  } else {
    nsites += neighbor4(i,&esites[nsites]);
    nsites += neighbor4(j,&esites[nsites]);
  }

  solve->update(nsites,esites,propensity);

  // sanity check on all propensity values

  /*
  printf("EVENT %d %d\n",i,j);
  for (m = 0; m < nlocal; m++) {
    if (fabs(propensity[m]-site_propensity(m)) > 1.0e-6) {
      printf("BAD PROP = %d %d %d %g %g\n",
	     id[i],id[j],id[m],propensity[m],site_propensity(m));
      error->one("BAD DONE");
    }
  }
  */

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 2nd neighbors of site I
------------------------------------------------------------------------- */

int AppDiffusion::neighbor2(int i, int *sites)
{
  int k,kk,m,mm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(m);
/*****************************/ 	
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(mm);
/*****************************/ 	
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  return nsites;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 3rd neighbors of site I
------------------------------------------------------------------------- */

int AppDiffusion::neighbor3(int i, int *sites)
{
  int k,kk,kkk,m,mm,mmm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(m);
/*****************************/ 	
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(mm);
/*****************************/ 	
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(mmm);
/*****************************/ 	
	  sites[nsites++] = isite;
	  echeck[isite] = 1;
	}
      }
    }
  }

  return nsites;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 4th neighbors of site I
------------------------------------------------------------------------- */

int AppDiffusion::neighbor4(int i, int *sites)
{
  int k,kk,kkk,kkkk,m,mm,mmm,mmmm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(m);
/*****************************/ 	
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(mm);
/*****************************/ 	
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(mmm);
/*****************************/ 	
	  sites[nsites++] = isite;
	  echeck[isite] = 1;
	}
	for (kkkk = 0; kkkk < numneigh[mmm]; kkkk++) {
	  mmmm = neighbor[mmm][kkkk];
	  isite = i2site[mmmm];
	  if (isite >= 0 && echeck[isite] == 0) {
	    propensity[isite] = site_propensity(mmmm);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	propensity_intval[isite] = site_propensity_intval(mmmm);
/*****************************/ 	
	    sites[nsites++] = isite;
	    echeck[isite] = 1;
	  }
	}
      }
    }
  }

  return nsites;
}

/* ---------------------------------------------------------------------- */

int AppDiffusion::ncoord(int i)
{
  int count = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) count++;
  return count;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusion::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppDiffusion::add_event(int i, int destination, 
			      double propensity, int eventflag)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;
  events[freeevent].propensity = propensity;
  events[freeevent].destination = destination;
  events[freeevent].style = eventflag;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
void AppDiffusion::add_event(int i, int destination, 
			      IntervalItem propensity_intval, double propensity, int eventflag)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;
  events[freeevent].propensity_intval = propensity_intval;
  events[freeevent].propensity = propensity;
  events[freeevent].destination = destination;
  events[freeevent].style = eventflag;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}
/*****************************/

/* ----------------------------------------------------------------------
   enumerate Schwoebel hop events centered around OCCUPIED site I
   assume mark array is currently cleared, use it, clear it when done
------------------------------------------------------------------------- */

int AppDiffusion::schwoebel_enumerate(int i, int *site)
{
  int j,k,m,jneigh,kneigh,count;

  int nhop = 0;

  // if coord(I) > Nmax, no hops possible

  count = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) count++;
  if (count > nsmax) return nhop;

  // mark first neighbors of site I as vacant = 1, occupied = 2

  for (jneigh = 0; jneigh < numneigh[i]; jneigh++) {
    j = neighbor[i][jneigh];
    if (lattice[j] == VACANT) mark[j] = 1;
    else if (lattice[j] == OCCUPIED) mark[j] = 2;
  }

  // loop over 1st and 2nd neighbors of site I to find possible hops to K
  // if K not vacant, no hop
  // if mark(K) = 1 or 2, K is a 1st neigh, not a 2nd neigh
  // if mark(K) = 30, already seen as a possible hop
  // mark(K) = 10 or 20, it has a 1st neigh that is vacant or occupied
  // mark(K) = 30, has both a vacant/occupied 1st neigh, so consider it
  // if coord(K) < Nmin, no hop possible
  // if all criteria met, then it is a candidate hop, add to site[]
  
  int nlist = 0;
  for (jneigh = 0; jneigh < numneigh[i]; jneigh++) {
    j = neighbor[i][jneigh];
    for (kneigh = 0; kneigh < numneigh[j]; kneigh++) {
      k = neighbor[j][kneigh];
      if (lattice[k] != VACANT) continue;
      if (mark[k] == 1 || lattice[k] == 2) continue;
      if (mark[k] == 30) continue;
      if (mark[k] == 10*mark[j]) continue;
      if (mark[k] == 0) marklist[nlist++] = k;
      mark[k] += 10*mark[j];
      if (mark[k] != 30) continue;

      count = 0;
      for (m = 0; m < numneigh[k]; m++)
	if (lattice[neighbor[k][m]] == OCCUPIED) count++;
      if (count < nsmin) continue;

      site[nhop++] = k;
    }
  }

  // clear marked sites, 1st and 2nd neighbors

  for (j = 0; j < numneigh[i]; j++) mark[neighbor[i][j]] = 0;
  for (k = 0; k < nlist; k++) mark[marklist[k]] = 0;

  return nhop;
}

/* ----------------------------------------------------------------------
   identify a VACANT site to deposit an atom
   return -1 if could not find a suitable site
------------------------------------------------------------------------- */

int AppDiffusion::find_deposition_site(RandomPark *random)
{
  // pick a random position at top of box

  double start[3];
  start[0] = domain->boxxlo + domain->xprd*random->uniform();
  if (dimension == 2) {
    start[1] = domain->boxyhi;
    start[2] = 0.0;
  } else {
    start[1] = domain->boxylo + domain->yprd*random->uniform();
    start[2] = domain->boxzhi;
  }

  // for each vacant site:
  // discard site if neighbor count not between coordlo and coordhi
  // find site whose projected distance is closest to start point

  int i,j,ncount;
  double dist2start;

  int closesite = -1;
  double closedist = 1.0e20;

  for (i = 0; i < nlocal; i++) {
    if (lattice[i] != VACANT) continue;
    ncount = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) ncount++;
    if (ncount < coordlo || ncount > coordhi) continue;

    if (exceed_limit(i,start,dist2start)) continue;
    if (dist2start < closedist) {
      closedist = dist2start;
      closesite = i;
    }
  }

  if (closesite < 0) ndeposit_failed++;
  else ndeposit++;

  return closesite;
}

/* ----------------------------------------------------------------------
   test if site M is within normal distance d0 from incident line
   if so, return 0 and dist2start, else return 1
   site M really becomes a periodic image in XY of M, adjusted via iprd/jprd
   dist2start = dist from site M to starting point of incident line
   dist2start is dist along incident line from start point to
     normal projection point of M
------------------------------------------------------------------------- */

int AppDiffusion::exceed_limit(int m, double *start, double &dist2start)
{
  int increment,iprd,jprd;

  iprd = jprd = 0;
  double d0sq = d0*d0;

  double distsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  double newdistsq = distsq_to_line(m,start,iprd-1,jprd,dist2start);
  if (newdistsq < distsq) increment = -1;
  else increment = 1;

  iprd += increment;
  newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  while (newdistsq < distsq) {
    distsq = newdistsq;
    iprd += increment;
    newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  }
  iprd -= increment;

  if (dimension == 3) {
    newdistsq = distsq_to_line(m,start,iprd,jprd-1,dist2start);
    if (newdistsq < distsq) increment = -1;
    else increment = 1;

    jprd += increment;
    newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
    while (newdistsq < distsq) {
      distsq = newdistsq;
      jprd += increment;
      newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
    }
  }
  jprd -= increment;

  if (distsq > d0sq) return 1;
  distsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  return 0;
}

/* ----------------------------------------------------------------------
   compute normal distsq from site M to incident line of deposition
   site M really becomes a periodic image in XY of M, adjusted via iprd/jprd
   also compute and return dist2start
   dist2start = dist from site M to starting point of incident line
   dist2start is dist along incident line from start point to
     normal projection point of M
------------------------------------------------------------------------- */

double AppDiffusion::distsq_to_line(int m, double *start,
				    int iprd, int jprd, double &dist2start)
{
  double dot,distsq;
  double delta[3],projection[3],offset[3];

  delta[0] = xyz[m][0] + iprd*domain->xprd - start[0];
  delta[1] = xyz[m][1] + jprd*domain->yprd - start[1];
  delta[2] = xyz[m][2] - start[2];
    
  dist2start = dir[0]*delta[0] + dir[1]*delta[1] + dir[2]*delta[2];
  projection[0] = dist2start*dir[0];
  projection[1] = dist2start*dir[1];
  projection[2] = dist2start*dir[2];
  
  offset[0] = delta[0] - projection[0];
  offset[1] = delta[1] - projection[1];
  offset[2] = delta[2] - projection[2];
  return offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2];
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 0 to nmax,
     (3) i* = 0 to nmax, (4) *j = 0 to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void AppDiffusion::bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),0);
    nhi = MIN(atoi(str),nmax);
  } else if (strlen(str) == 1) {
    nlo = 0;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 0;
    nhi = MIN(atoi(ptr+1),nmax);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),0);
    nhi = nmax;
  } else {
    nlo = MAX(atoi(str),0);
    nhi = MIN(atoi(ptr+1),nmax);
  }
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppDiffusion::allocate_data()
{
  // for no_energy or linear:
  //   make esites large enough for 2 sites and their 1,2 neighbors
  //   do not need psites
  // for nonlinear:
  //   make esites large enough for 2 sites and their 1,2,3 neighbors
  //   make psites large enough for 2 sites and their 1st neighbors
  // Schwoebel hops add one level of neighbor dependence to esites

  if ((engstyle == NO_ENERGY || engstyle == LINEAR) && hopstyle == NNHOP) {
    int emax = 1 + maxneigh + maxneigh*maxneigh;
    esites = new int[2*emax];
    psites = NULL;
  } else if ((engstyle == NO_ENERGY || engstyle == LINEAR) && 
	     hopstyle == SCHWOEBEL) {
    int emax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
    esites = new int[2*emax];
    psites = NULL;
  } else if (engstyle == NONLINEAR && hopstyle == NNHOP) {
    int emax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
    int pmax = 1 + maxneigh;
    esites = new int[2*emax];
    psites = new int[2*pmax];
  } else if (engstyle == NONLINEAR && hopstyle == SCHWOEBEL) {
    int emax = 1 + maxneigh + maxneigh*maxneigh + 
      maxneigh*maxneigh*maxneigh + maxneigh*maxneigh*maxneigh*maxneigh;
    int pmax = 1 + maxneigh;
    esites = new int[2*emax];
    psites = new int[2*pmax];
  }

  echeck = new int[nlocal+nghost];
  pcheck = new int[nlocal+nghost];

  firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");

  ecoord = new double[maxneigh+1];
  for (int i = 0; i <= maxneigh; i++) ecoord[i] = 0.0;

  hbarrier = 
    memory->create_2d_double_array(maxneigh+1,maxneigh+1,"app:hbarrier");
  sbarrier = 
    memory->create_2d_double_array(maxneigh+1,maxneigh+1,"app:sbarrier");

  for (int i = 0; i <= maxneigh; i++)
    for (int j = 0; j <= maxneigh; j++)
      hbarrier[i][j] = sbarrier[i][j] = 0.0;

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	memory->create_2d_T_array(hbarrier_intval,maxneigh+1,maxneigh+1,
			"diffusion:hbarrier_intval");
	memory->create_2d_T_array(sbarrier_intval,maxneigh+1,maxneigh+1,
			"diffusion:sbarrier_intval");
	for (int i = 0; i <= maxneigh; i++)
		for (int j = 0; j <= maxneigh; j++)
		{
			hbarrier_intval[i][j].inf
				= hbarrier_intval[i][j].sup
				= hbarrier_intval[i][j].wid = 0.0;
			sbarrier_intval[i][j].inf
				= sbarrier_intval[i][j].sup
				= sbarrier_intval[i][j].wid = 0.0;
		}
/*****************************/	

  hopsite = new int[maxneigh*maxneigh + maxneigh];
  marklist = new int[maxneigh*maxneigh];

  mark = NULL;
  if (hopstyle == SCHWOEBEL)
    mark = (int *) memory->smalloc((nlocal+nghost)*sizeof(int),"app:mark");
  if (mark)
    for (int i = 0; i < nlocal+nghost; i++) mark[i] = 0;
}
