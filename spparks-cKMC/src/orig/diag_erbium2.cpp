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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_erbium.h"
#include "app.h"
#include "app_erbium.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{VACANCY,ERBIUM,HYDROGEN,HELIUM,ELEMENT4};   // same as AppErbium
enum{VAC,ELEM1,ELEM2,ELEM3,ELEM4,EVENTS,ONE,TWO,THREE,FOUR,CTRL};

/* ---------------------------------------------------------------------- */

DiagErbium::DiagErbium(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"erbium") != 0)
    error->all("Diag_style erbium requires app_style erbium");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all("Illegal diag_style erbium command");
  }

  if (nlist == 0) error->all("Illegal diag_style erbium command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagErbium::~DiagErbium()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagErbium::init()
{
  apperbium = (AppErbium *) app;

  int none = apperbium->none;
  int ntwo = apperbium->ntwo;
  int nthree = apperbium->nthree;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
  int nfour = apperbium->nfour;
  int nctrl = apperbium->numControlSpecies;
/*****************************/
  for (int i = 0; i < nlist; i++) {
  
/*********commented by Yan Wang and replaced by the following *******
    if (strcmp(list[i],"er") == 0) which[i] = ER;
    else if (strcmp(list[i],"h") == 0) which[i] = H;
    else if (strcmp(list[i],"he") == 0) which[i] = HE;
    else if (strcmp(list[i],"vac") == 0) which[i] = VAC;
    else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
*****************************************************************/	
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
    int ispecies = apperbium->find_species(list[i]);
	if (ispecies != -1) which[i] = ispecies;
	else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
/******************************/	
    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all("Invalid value setting in diag_style erbium: s");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all("Invalid value setting in diag_style erbium: d");
      index[i] = n - 1;
    } else if (list[i][0] == 't') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all("Invalid value setting in diag_style erbium: t");
      index[i] = n - 1;
    }
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
	else if (list[i][0] == 'q') {
      which[i] = FOUR;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nfour) 
	error->all("Invalid value setting in diag_style erbium: q");
      index[i] = n - 1;
    }
	else if (list[i][0] == 'c') {
      which[i] = CTRL;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nctrl) 
	error->all("Invalid value setting in diag_style erbium: c");
      index[i] = n - 1;
    }
/****************************/	
	else error->all("Invalid value setting in diag_style erbium");
  }

  siteflag = 0;
  for (int i = 0; i < nlist; i++)
/****commented by Yan Wang  and replaced by the following ********
    if (which[i] == ELEM1 || which[i] == ELEM2 || which[i] == ELEM3 || which[i] == VAC)
**************************************************************/    
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
    if (which[i] == ELEM1 || which[i] == ELEM2 || which[i] == ELEM3 ||
		which[i] == ELEM4 || which[i] == VAC || which[i] == CTRL)
/****************************/		
      siteflag = 1;

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagErbium::compute()
{
  int sites[6],ivalue;

  if (siteflag) {
/****commented by Yan Wang  and replaced by the following ********
    sites[ERBIUM] = sites[HYDROGEN] = sites[HELIUM] = sites[VACANCY] = 0;
**************************************************************/    
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
     sites[ERBIUM] = sites[HYDROGEN] = sites[HELIUM] = sites[VACANCY] = sites[ELEMENT4] = 0;
/****************************/	
    int *element = apperbium->element;
    int nlocal = apperbium->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == ELEM1) ivalue = sites[ERBIUM];
    else if (which[i] == ELEM2) ivalue = sites[HYDROGEN];
    else if (which[i] == ELEM3) ivalue = sites[HELIUM];
    else if (which[i] == ELEM4) ivalue = sites[ELEMENT4];
    else if (which[i] == VAC) ivalue = sites[VACANCY];
    else if (which[i] == EVENTS) ivalue = apperbium->nevents;
    else if (which[i] == ONE) ivalue = apperbium->scount[index[i]];
    else if (which[i] == TWO) ivalue = apperbium->dcount[index[i]];
    else if (which[i] == THREE) ivalue = apperbium->tcount[index[i]];
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
    else if (which[i] == FOUR) ivalue = apperbium->qcount[index[i]];
    else if (which[i] == CTRL) ivalue = apperbium->currCtrlSpecSite[index[i]];
/********************************/   
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErbium::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErbium::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
