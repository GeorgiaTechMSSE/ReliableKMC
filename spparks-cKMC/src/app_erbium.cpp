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
#include "mpi.h"
#include "app_erbium.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <cstring>
#include <stdlib.h>

using namespace SPPARKS_NS;

enum{NOOP,FCC,OCTA,TETRA};
enum{VACANCY,ERBIUM,HYDROGEN,HELIUM};      // same as DiagErbium

#define DELTAEVENT 100000
#define PI 3.1415926535897932
#define EPSILON	0.000001

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
#define MEM_DELTA	4

int compare_site(const void *a, const void *b)  
{
   const struct csite *const *x = (struct csite **)a;
   const struct csite *const *y = (struct csite **)b;
   if ( (*x)->keyvalue > (*y)->keyvalue ) return 1;
   if ( (*x)->keyvalue < (*y)->keyvalue ) return -1;
   return 0;
}

/*********************************/

/* ---------------------------------------------------------------------- */

AppErbium::AppErbium(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  create_arrays();

  if (narg != 1) error->all("Illegal app_style command");

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = 0;
  srate = drate = trate = NULL;
  spropensity = dpropensity = tpropensity = NULL;
  stype = sinput = soutput = NULL;
  dtype = dinput = doutput = NULL;
  ttype = tinput = toutput = NULL;
  scount = dcount = tcount = NULL;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/
  //parameters for events with 4 (quadruple) neighbors involved
  nfour = 0;
  qrate = NULL;
  qpropensity = NULL;
  qtype = qinput = qoutput = NULL;
  qcount = NULL; 
  
  //1D arrays storing the indices of control events corresponding to "controlEvent_list"
  // single, double, triple, quadruple events, store -1 if not controlled event
  sctrlIndex = NULL;
  dctrlIndex = NULL; 
  tctrlIndex = NULL;
  qctrlIndex = NULL;

  numControlSpecies = 0;			//# of control species
  controlSpec_list = NULL;		//list of the control species
  controlSpec_rate = NULL;		//control species reaction rate
  controlSpec_dir = NULL;		//control species reaction direction vector
  controlSpec_coord = NULL;		//lower and upper bounds of neighbors that control species would react 

  numControlEvents = 0;			//# of controlled events
  controlEvent_list = NULL;		//list of the indices for controlled events 
								// numControlEventsX2 array, 1st column: type of event; 2nd column: index
  controlEvent_dir = NULL;		// numControlEventsX7 array, controlled reaction direction vector
								// first 3 columns: mean values of x,y,z directions, 
								// next 3 columns: standard deviations of direction vectors 
								// last column: the allowed angular deviation from the above vector
  controlEvent_coord = NULL;		////numControlEventsX2 array, lower and upper bounds of neighbors that controlled event would react 
  
  ctrlStartTime = NULL;		//the time when the control specie reactions start, 0 by default
  currCtrlSpecSite = NULL;		//keep track of current index of control sites as the time advances
  numCtrlSpecSites = NULL;		// # of sites where control species reside
  ctrlSpecSites = NULL;		//the list of sites where control species reside initially
  csSiteArray = NULL;
  
  ctrlSpecSitesList = NULL;	//the actual 2D list of control species sites 
  maxNumCtrlSites = 0;		//keep track of the maxium number of control species sites
/*****************************/  
}

/* ---------------------------------------------------------------------- */

AppErbium::~AppErbium()
{

  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);

  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(trate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(tpropensity);
  memory->sfree(stype);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->destroy_2d_int_array(dtype);
  memory->destroy_2d_int_array(dinput);
  memory->destroy_2d_int_array(doutput);
  memory->destroy_2d_int_array(ttype);
  memory->destroy_2d_int_array(tinput);
  memory->destroy_2d_int_array(toutput);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(tcount);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/ 
  memory->sfree(qrate);
  memory->sfree(qpropensity);
  memory->destroy_2d_int_array(qtype);
  memory->destroy_2d_int_array(qinput);
  memory->destroy_2d_int_array(qoutput);
  memory->sfree(qcount);
  
  memory->sfree(sctrlIndex);
  memory->sfree(dctrlIndex);
  memory->sfree(tctrlIndex);
  memory->sfree(qctrlIndex);
  
  memory->destroy_2d_int_array(controlSpec_list);
  memory->destroy_2d_double_array(controlSpec_dir);
  memory->destroy_2d_int_array(controlSpec_coord);
  memory->destroy_1d_T_array(ctrlSpecSites, 0);
  memory->sfree(ctrlSpecSites_addr);
  
  memory->destroy_2d_int_array(controlEvent_list);
  memory->destroy_2d_double_array(controlEvent_dir);
  memory->destroy_2d_int_array(controlEvent_coord);

  memory->destroy_2d_double_array(csSiteArray);

  memory->sfree(ctrlStartTime);	
  memory->sfree(currCtrlSpecSite);	
  memory->sfree(numCtrlSpecSites);

  memory->destroy_2d_T_array(ctrlSpecSitesList);
/*******************************/  
}

/* ---------------------------------------------------------------------- */
/**************************************************************
 * Commented by Yan Wang (July  2010) and replaced by the following
 ************************************************************  
void AppErbium::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all("Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 5) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) stype[none] = FCC;
      else if (strcmp(arg[1],"oct") == 0) stype[none] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) stype[none] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"er") == 0) sinput[none] = ERBIUM;
      else if (strcmp(arg[2],"h") == 0) sinput[none] = HYDROGEN;
      else if (strcmp(arg[2],"he") == 0) sinput[none] = HELIUM;
      else if (strcmp(arg[2],"vac") == 0) sinput[none] = VACANCY;
      else error->all("Illegal event command");

      srate[none] = atof(arg[3]);

      if (strcmp(arg[4],"er") == 0) soutput[none] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) soutput[none] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) soutput[none] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) soutput[none] = VACANCY;
      else error->all("Illegal event command");

      none++;
      
    } else if (rstyle == 2) {
      if (narg != 8) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) dtype[ntwo][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) dtype[ntwo][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) dtype[ntwo][0] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) dtype[ntwo][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) dtype[ntwo][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) dtype[ntwo][1] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[3],"er") == 0) dinput[ntwo][0] = ERBIUM;
      else if (strcmp(arg[3],"h") == 0) dinput[ntwo][0] = HYDROGEN;
      else if (strcmp(arg[3],"he") == 0) dinput[ntwo][0] = HELIUM;
      else if (strcmp(arg[3],"vac") == 0) dinput[ntwo][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[4],"er") == 0) dinput[ntwo][1] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) dinput[ntwo][1] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) dinput[ntwo][1] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) dinput[ntwo][1] = VACANCY;
      else error->all("Illegal event command");

      drate[ntwo] = atof(arg[5]);

      if (strcmp(arg[6],"er") == 0) doutput[ntwo][0] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) doutput[ntwo][0] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) doutput[ntwo][0] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) doutput[ntwo][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[7],"er") == 0) doutput[ntwo][1] = ERBIUM;
      else if (strcmp(arg[7],"h") == 0) doutput[ntwo][1] = HYDROGEN;
      else if (strcmp(arg[7],"he") == 0) doutput[ntwo][1] = HELIUM;
      else if (strcmp(arg[7],"vac") == 0) doutput[ntwo][1] = VACANCY;
      else error->all("Illegal event command");

      ntwo++;

    } else if (rstyle == 3) {
      if (narg != 11) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) ttype[nthree][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) ttype[nthree][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) ttype[nthree][0] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) ttype[nthree][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) ttype[nthree][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) ttype[nthree][1] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[3],"fcc") == 0) ttype[nthree][2] = FCC;
      else if (strcmp(arg[3],"oct") == 0) ttype[nthree][2] = OCTA;
      else if (strcmp(arg[3],"tet") == 0) ttype[nthree][2] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[4],"er") == 0) tinput[nthree][0] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) tinput[nthree][0] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) tinput[nthree][0] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) tinput[nthree][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[5],"er") == 0) tinput[nthree][1] = ERBIUM;
      else if (strcmp(arg[5],"h") == 0) tinput[nthree][1] = HYDROGEN;
      else if (strcmp(arg[5],"he") == 0) tinput[nthree][1] = HELIUM;
      else if (strcmp(arg[5],"vac") == 0) tinput[nthree][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[6],"er") == 0) tinput[nthree][2] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) tinput[nthree][2] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) tinput[nthree][2] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) tinput[nthree][2] = VACANCY;
      else error->all("Illegal event command");

      trate[nthree] = atof(arg[7]);

      if (strcmp(arg[8],"er") == 0) toutput[nthree][0] = ERBIUM;
      else if (strcmp(arg[8],"h") == 0) toutput[nthree][0] = HYDROGEN;
      else if (strcmp(arg[8],"he") == 0) toutput[nthree][0] = HELIUM;
      else if (strcmp(arg[8],"vac") == 0) toutput[nthree][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[9],"er") == 0) toutput[nthree][1] = ERBIUM;
      else if (strcmp(arg[9],"h") == 0) toutput[nthree][1] = HYDROGEN;
      else if (strcmp(arg[9],"he") == 0) toutput[nthree][1] = HELIUM;
      else if (strcmp(arg[9],"vac") == 0) toutput[nthree][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[10],"er") == 0) toutput[nthree][2] = ERBIUM;
      else if (strcmp(arg[10],"h") == 0) toutput[nthree][2] = HYDROGEN;
      else if (strcmp(arg[10],"he") == 0) toutput[nthree][2] = HELIUM;
      else if (strcmp(arg[10],"vac") == 0) toutput[nthree][2] = VACANCY;
      else error->all("Illegal event command");

      nthree++;

    }
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************  
	else if (rstyle == 4) {
      if (narg != 14) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) qtype[nfour][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) qtype[nfour][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) qtype[nfour][0] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) qtype[nfour][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) qtype[nfour][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) qtype[nfour][1] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[3],"fcc") == 0) qtype[nfour][2] = FCC;
      else if (strcmp(arg[3],"oct") == 0) qtype[nfour][2] = OCTA;
      else if (strcmp(arg[3],"tet") == 0) qtype[nfour][2] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[4],"fcc") == 0) qtype[nfour][3] = FCC;
      else if (strcmp(arg[4],"oct") == 0) qtype[nfour][3] = OCTA;
      else if (strcmp(arg[4],"tet") == 0) qtype[nfour][3] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[5],"er") == 0) qinput[nfour][0] = ERBIUM;
      else if (strcmp(arg[5],"h") == 0) qinput[nfour][0] = HYDROGEN;
      else if (strcmp(arg[5],"he") == 0) qinput[nfour][0] = HELIUM;
      else if (strcmp(arg[5],"vac") == 0) qinput[nfour][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[6],"er") == 0) qinput[nfour][1] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) qinput[nfour][1] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) qinput[nfour][1] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) qinput[nfour][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[7],"er") == 0) qinput[nfour][2] = ERBIUM;
      else if (strcmp(arg[7],"h") == 0) qinput[nfour][2] = HYDROGEN;
      else if (strcmp(arg[7],"he") == 0) qinput[nfour][2] = HELIUM;
      else if (strcmp(arg[7],"vac") == 0) qinput[nfour][2] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[8],"er") == 0) qinput[nfour][3] = ERBIUM;
      else if (strcmp(arg[8],"h") == 0) qinput[nfour][3] = HYDROGEN;
      else if (strcmp(arg[8],"he") == 0) qinput[nfour][3] = HELIUM;
      else if (strcmp(arg[8],"vac") == 0) qinput[nfour][3] = VACANCY;
      else error->all("Illegal event command");

      qrate[nfour] = atof(arg[9]);

      if (strcmp(arg[10],"er") == 0) qoutput[nfour][0] = ERBIUM;
      else if (strcmp(arg[10],"h") == 0) qoutput[nfour][0] = HYDROGEN;
      else if (strcmp(arg[10],"he") == 0) qoutput[nfour][0] = HELIUM;
      else if (strcmp(arg[10],"vac") == 0) qoutput[nfour][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[11],"er") == 0) qoutput[nfour][1] = ERBIUM;
      else if (strcmp(arg[11],"h") == 0) qoutput[nfour][1] = HYDROGEN;
      else if (strcmp(arg[11],"he") == 0) qoutput[nfour][1] = HELIUM;
      else if (strcmp(arg[11],"vac") == 0) qoutput[nfour][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[12],"er") == 0) qoutput[nfour][2] = ERBIUM;
      else if (strcmp(arg[12],"h") == 0) qoutput[nfour][2] = HYDROGEN;
      else if (strcmp(arg[12],"he") == 0) qoutput[nfour][2] = HELIUM;
      else if (strcmp(arg[12],"vac") == 0) qoutput[nfour][2] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[13],"er") == 0) qoutput[nfour][3] = ERBIUM;
      else if (strcmp(arg[13],"h") == 0) qoutput[nfour][3] = HYDROGEN;
      else if (strcmp(arg[13],"he") == 0) qoutput[nfour][3] = HELIUM;
      else if (strcmp(arg[13],"vac") == 0) qoutput[nfour][3] = VACANCY;
      else error->all("Illegal event command");

      nfour++;

    }
/***************************************
	else error->all("Illegal event command");
  } else error->all("Unrecognized command");
}
*****************************************************/


/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
void AppErbium::input_app(char *command, int narg, char **arg)
{

  if (strcmp(command,"add_species") == 0) add_species(narg,arg);
  else if (strcmp(command,"control_species") == 0) add_control_species(narg,arg);
  else if (strcmp(command,"event") == 0) parse_event(narg,arg);
  else if (strcmp(command,"control_event") == 0) add_control_event(narg,arg);
  else error->all("Unrecognized command");
/*******
for(int j=0; j<none; j++)	
	printf("sinput[%d]=%d soutput[%d]=%d\n",j,sinput[j],j,soutput[j]);
for(int j=0; j<ntwo; j++)
	for(int k=0; k<2; k++)
		printf("dinput[%d][%d]=%d doutput[%d][%d]=%d\n",j,k,dinput[j][k],j,k,doutput[j][k]);	
for(int j=0; j<nthree; j++)	
	for(int k=0; k<3; k++)
		printf("tinput[%d][%d]=%d toutput[%d][%d]=%d\n",j,k,tinput[j][k],j,k,toutput[j][k]);	
for(int j=0; j<nfour; j++)	
	for(int k=0; k<4; k++)
		printf("qinput[%d][%d]=%d qoutput[%d][%d]=%d\n",j,k,qinput[j][k],j,k,qoutput[j][k]);	
*******/
}

void AppErbium::parse_event(int narg, char **arg)
{
	int ispecies;
	
    if (narg < 1) error->all("Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 5) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) stype[none] = FCC;
      else if (strcmp(arg[1],"oct") == 0) stype[none] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) stype[none] = TETRA;
      else error->all("Illegal event command");
/*      if (strcmp(arg[2],"er") == 0) sinput[none] = ERBIUM;
      else if (strcmp(arg[2],"h") == 0) sinput[none] = HYDROGEN;
      else if (strcmp(arg[2],"he") == 0) sinput[none] = HELIUM;
      else if (strcmp(arg[2],"vac") == 0) sinput[none] = VACANCY;
      else error->all("Illegal event command");
*/
	    ispecies = find_species(arg[2]);
	    if (ispecies == -1) error->all("Unknown species in event command");
		sinput[none] = ispecies;

	srate[none] = atof(arg[3]);

/*	
      if (strcmp(arg[4],"er") == 0) soutput[none] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) soutput[none] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) soutput[none] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) soutput[none] = VACANCY;
      else error->all("Illegal event command");
*/
    ispecies = find_species(arg[4]);
    if (ispecies == -1) error->all("Unknown species in event command");
	soutput[none] = ispecies;

      none++;
      
    } 
	
	else if (rstyle == 2) {
      if (narg != 8) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) dtype[ntwo][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) dtype[ntwo][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) dtype[ntwo][0] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) dtype[ntwo][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) dtype[ntwo][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) dtype[ntwo][1] = TETRA;
      else error->all("Illegal event command");
	  
/*      if (strcmp(arg[3],"er") == 0) dinput[ntwo][0] = ERBIUM;
      else if (strcmp(arg[3],"h") == 0) dinput[ntwo][0] = HYDROGEN;
      else if (strcmp(arg[3],"he") == 0) dinput[ntwo][0] = HELIUM;
      else if (strcmp(arg[3],"vac") == 0) dinput[ntwo][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[4],"er") == 0) dinput[ntwo][1] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) dinput[ntwo][1] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) dinput[ntwo][1] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) dinput[ntwo][1] = VACANCY;
      else error->all("Illegal event command");
*/
		for(int i=0; i<rstyle; i++)
		{
		    ispecies = find_species(arg[1+rstyle+i]);
		    if (ispecies == -1) error->all("Unknown species in event command");
			dinput[ntwo][i] = ispecies;
		}
   
      drate[ntwo] = atof(arg[5]);

/*      if (strcmp(arg[6],"er") == 0) doutput[ntwo][0] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) doutput[ntwo][0] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) doutput[ntwo][0] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) doutput[ntwo][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[7],"er") == 0) doutput[ntwo][1] = ERBIUM;
      else if (strcmp(arg[7],"h") == 0) doutput[ntwo][1] = HYDROGEN;
      else if (strcmp(arg[7],"he") == 0) doutput[ntwo][1] = HELIUM;
      else if (strcmp(arg[7],"vac") == 0) doutput[ntwo][1] = VACANCY;
      else error->all("Illegal event command");
*/
		for(int i=0; i<rstyle; i++)
		{
		    ispecies = find_species(arg[1+rstyle+rstyle+1+i]);
		    if (ispecies == -1) error->all("Unknown species in event command");
			doutput[ntwo][i] = ispecies;
		}

      ntwo++;

    }

	else if (rstyle == 3) {
      if (narg != 11) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) ttype[nthree][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) ttype[nthree][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) ttype[nthree][0] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) ttype[nthree][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) ttype[nthree][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) ttype[nthree][1] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[3],"fcc") == 0) ttype[nthree][2] = FCC;
      else if (strcmp(arg[3],"oct") == 0) ttype[nthree][2] = OCTA;
      else if (strcmp(arg[3],"tet") == 0) ttype[nthree][2] = TETRA;
      else error->all("Illegal event command");
/*      if (strcmp(arg[4],"er") == 0) tinput[nthree][0] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) tinput[nthree][0] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) tinput[nthree][0] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) tinput[nthree][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[5],"er") == 0) tinput[nthree][1] = ERBIUM;
      else if (strcmp(arg[5],"h") == 0) tinput[nthree][1] = HYDROGEN;
      else if (strcmp(arg[5],"he") == 0) tinput[nthree][1] = HELIUM;
      else if (strcmp(arg[5],"vac") == 0) tinput[nthree][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[6],"er") == 0) tinput[nthree][2] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) tinput[nthree][2] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) tinput[nthree][2] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) tinput[nthree][2] = VACANCY;
      else error->all("Illegal event command");
*/
		for(int i=0; i<rstyle; i++)
		{
		    ispecies = find_species(arg[1+rstyle+i]);
		    if (ispecies == -1) error->all("Unknown species in event command");
			tinput[nthree][i] = ispecies;
		}

      trate[nthree] = atof(arg[7]);

/*      if (strcmp(arg[8],"er") == 0) toutput[nthree][0] = ERBIUM;
      else if (strcmp(arg[8],"h") == 0) toutput[nthree][0] = HYDROGEN;
      else if (strcmp(arg[8],"he") == 0) toutput[nthree][0] = HELIUM;
      else if (strcmp(arg[8],"vac") == 0) toutput[nthree][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[9],"er") == 0) toutput[nthree][1] = ERBIUM;
      else if (strcmp(arg[9],"h") == 0) toutput[nthree][1] = HYDROGEN;
      else if (strcmp(arg[9],"he") == 0) toutput[nthree][1] = HELIUM;
      else if (strcmp(arg[9],"vac") == 0) toutput[nthree][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[10],"er") == 0) toutput[nthree][2] = ERBIUM;
      else if (strcmp(arg[10],"h") == 0) toutput[nthree][2] = HYDROGEN;
      else if (strcmp(arg[10],"he") == 0) toutput[nthree][2] = HELIUM;
      else if (strcmp(arg[10],"vac") == 0) toutput[nthree][2] = VACANCY;
      else error->all("Illegal event command");
*/
		for(int i=0; i<rstyle; i++)
		{
		    ispecies = find_species(arg[1+rstyle+rstyle+1+i]);
		    if (ispecies == -1) error->all("Unknown species in event command");
			toutput[nthree][i] = ispecies;
		}
      nthree++;

    }
 
	else if (rstyle == 4) {
      if (narg != 14) error->all("Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) qtype[nfour][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) qtype[nfour][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) qtype[nfour][0] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) qtype[nfour][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) qtype[nfour][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) qtype[nfour][1] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[3],"fcc") == 0) qtype[nfour][2] = FCC;
      else if (strcmp(arg[3],"oct") == 0) qtype[nfour][2] = OCTA;
      else if (strcmp(arg[3],"tet") == 0) qtype[nfour][2] = TETRA;
      else error->all("Illegal event command");
      if (strcmp(arg[4],"fcc") == 0) qtype[nfour][3] = FCC;
      else if (strcmp(arg[4],"oct") == 0) qtype[nfour][3] = OCTA;
      else if (strcmp(arg[4],"tet") == 0) qtype[nfour][3] = TETRA;
      else error->all("Illegal event command");
/*      if (strcmp(arg[5],"er") == 0) qinput[nfour][0] = ERBIUM;
      else if (strcmp(arg[5],"h") == 0) qinput[nfour][0] = HYDROGEN;
      else if (strcmp(arg[5],"he") == 0) qinput[nfour][0] = HELIUM;
      else if (strcmp(arg[5],"vac") == 0) qinput[nfour][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[6],"er") == 0) qinput[nfour][1] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) qinput[nfour][1] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) qinput[nfour][1] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) qinput[nfour][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[7],"er") == 0) qinput[nfour][2] = ERBIUM;
      else if (strcmp(arg[7],"h") == 0) qinput[nfour][2] = HYDROGEN;
      else if (strcmp(arg[7],"he") == 0) qinput[nfour][2] = HELIUM;
      else if (strcmp(arg[7],"vac") == 0) qinput[nfour][2] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[8],"er") == 0) qinput[nfour][3] = ERBIUM;
      else if (strcmp(arg[8],"h") == 0) qinput[nfour][3] = HYDROGEN;
      else if (strcmp(arg[8],"he") == 0) qinput[nfour][3] = HELIUM;
      else if (strcmp(arg[8],"vac") == 0) qinput[nfour][3] = VACANCY;
      else error->all("Illegal event command");
*/
		for(int i=0; i<rstyle; i++)
		{
		    ispecies = find_species(arg[1+rstyle+i]);
		    if (ispecies == -1) error->all("Unknown species in event command");
			qinput[nfour][i] = ispecies;
		}

      qrate[nfour] = atof(arg[9]);

/*      if (strcmp(arg[10],"er") == 0) qoutput[nfour][0] = ERBIUM;
      else if (strcmp(arg[10],"h") == 0) qoutput[nfour][0] = HYDROGEN;
      else if (strcmp(arg[10],"he") == 0) qoutput[nfour][0] = HELIUM;
      else if (strcmp(arg[10],"vac") == 0) qoutput[nfour][0] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[11],"er") == 0) qoutput[nfour][1] = ERBIUM;
      else if (strcmp(arg[11],"h") == 0) qoutput[nfour][1] = HYDROGEN;
      else if (strcmp(arg[11],"he") == 0) qoutput[nfour][1] = HELIUM;
      else if (strcmp(arg[11],"vac") == 0) qoutput[nfour][1] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[12],"er") == 0) qoutput[nfour][2] = ERBIUM;
      else if (strcmp(arg[12],"h") == 0) qoutput[nfour][2] = HYDROGEN;
      else if (strcmp(arg[12],"he") == 0) qoutput[nfour][2] = HELIUM;
      else if (strcmp(arg[12],"vac") == 0) qoutput[nfour][2] = VACANCY;
      else error->all("Illegal event command");
      if (strcmp(arg[13],"er") == 0) qoutput[nfour][3] = ERBIUM;
      else if (strcmp(arg[13],"h") == 0) qoutput[nfour][3] = HYDROGEN;
      else if (strcmp(arg[13],"he") == 0) qoutput[nfour][3] = HELIUM;
      else if (strcmp(arg[13],"vac") == 0) qoutput[nfour][3] = VACANCY;
      else error->all("Illegal event command");
*/
		for(int i=0; i<rstyle; i++)
		{
		    ispecies = find_species(arg[1+rstyle+rstyle+1+i]);
		    if (ispecies == -1) error->all("Unknown species in event command");
			qoutput[nfour][i] = ispecies;
		}

      nfour++;

    }

	else error->all("Illegal event command");


}


void AppErbium::add_species(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal species command");

  // grow species arrays

  int n = numspecies + narg;
  sname = (char **) memory->srealloc(sname,n*sizeof(char *),
					  "chemistry:sname");
  pcount = (int *) memory->srealloc(pcount,n*sizeof(int),
					 "chemistry:pcount");

  for (int iarg = 0; iarg < narg; iarg++) {
    if (find_species(arg[iarg]) >= 0) {
      char *str = new char[128];
      sprintf(str,"Species ID %s already exists",arg[iarg]);
      error->all(str);
    }
    int nlen = strlen(arg[iarg]) + 1;
    sname[numspecies+iarg] = new char[nlen];
    strcpy(sname[numspecies+iarg],arg[iarg]);
    pcount[numspecies+iarg] = 0;
//printf("sname[%d]=%s\n",numspecies+iarg,sname[numspecies+iarg]);  
  }
  numspecies += narg;

}

void AppErbium::add_control_species(int narg, char **arg)
{
    if (narg != 9) error->all("Illegal control_species command: wrong number of arguments");
	int ispecie = find_species(arg[0]);
    if (ispecie == -1) {
      char *str = new char[128];
      sprintf(str,"Source control specie %s does not exist",arg[0]);
      error->all(str);
    }
	else {
	
	int n = numControlSpecies + 1;
	controlSpec_list = memory->grow_2d_int_array(controlSpec_list,n,2,
					 "canm:controlSpec_list");
	controlSpec_rate = (double *) memory->srealloc(controlSpec_rate,n*sizeof(double),
					  "canm:controlSpec_rate");
	controlSpec_dir = memory->grow_2d_double_array(controlSpec_dir,n,3,
					  "canm:controlSpec_dir");
	controlSpec_coord = memory->grow_2d_int_array(controlSpec_coord,n,2,
					  "canm:controlSpec_coord");
	
	ctrlStartTime = (double *) memory->srealloc(ctrlStartTime,n*sizeof(double),
						"canm:ctrlStartTime");			//the time when the controll reaction start
	numCtrlSpecSites = (int *)  memory->srealloc(numCtrlSpecSites,n*sizeof(int),
						"canm:numCtrlSpecSites");			//# of sites where control species reside
	currCtrlSpecSite = (int *)  memory->srealloc(currCtrlSpecSite,n*sizeof(int),
						"canm:currCtrlSpecSite");		//keep track of current index of control sites as the time advances

						
	//// add initial default values 
	ctrlStartTime[numControlSpecies] = 0.0;
	numCtrlSpecSites[numControlSpecies] = 0;
	currCtrlSpecSite[numControlSpecies] = 0;
	
	////// add source specie ID
	controlSpec_list[numControlSpecies][0] = ispecie;
	
//    deprate_intval.inf = atof(arg[0]);
//	deprate_intval.sup = atof(arg[1]);
//	deprate = (deprate_intval.inf+deprate_intval.sup)/2.0;

	////// add control specie reaction reate
	controlSpec_rate[numControlSpecies] = atof(arg[1]);		//control species reaction rate
    if (controlSpec_rate[numControlSpecies] < 0.0)
		error->all("Illegal control_species command: reaction rate should not be negative");

	////// add target specie ID	
	int ispecie2 = find_species(arg[2]);
    if (ispecie2 == -1) {
      char *str = new char[128];
      sprintf(str,"Target control specie %s does not exist",arg[2]);
      error->all(str);
    }
	
	controlSpec_list[numControlSpecies][1] = ispecie2;
	
	///// add the direction of control reaction
    controlSpec_dir[numControlSpecies][0] = atof(arg[3]);		//control species reaction direction vector - x
    controlSpec_dir[numControlSpecies][1] = atof(arg[4]);		//control species reaction direction vector - y
    controlSpec_dir[numControlSpecies][2] = atof(arg[5]);		//control species reaction direction vector - z
/*    if (domain->dimension == 2 && 
		(controlSpec_dir[numControlSpecies][1] >= 0.0 || 
		 controlSpec_dir[numControlSpecies][2] != 0.0) )
      error->all("Illegal control_species command: reaction direction does not match dimension");
    if (domain->dimension == 3 && controlSpec_dir[numControlSpecies][2] >= 0.0)
      error->all("Illegal control_species command: reaction direction does not match dimension");
*/
	  
	//control species react when the number of neighbors is greater than this lower bound
    controlSpec_coord[numControlSpecies][0] = atoi(arg[6]);	
	//control species react when the number of neighbors is less than this upper bound	
    controlSpec_coord[numControlSpecies][1] = atoi(arg[7]);			
    if (controlSpec_coord[numControlSpecies][0] < 0 || 
		controlSpec_coord[numControlSpecies][1] > maxneigh || 
		controlSpec_coord[numControlSpecies][0] > controlSpec_coord[numControlSpecies][1])
      error->all("Illegal control_species command: the number of neighbors exceeds maximum limit");
	

	///// add the start time of the control reaction
    ctrlStartTime[numControlSpecies] = atof(arg[8]);		//the start time of the controlled reaction
	if (numControlSpecies > 1 && ctrlStartTime[numControlSpecies]<ctrlStartTime[numControlSpecies-1])
      error->all("Illegal control_species command: the control species start time should be non-decreasing");
	
	numControlSpecies++;
/*	
printf("control_species: %d rate=%f dir=[%f %f %f] cood=(%d %d)",
	controlSpec_list[numControlSpecies],controlSpec_rate[numControlSpecies],
	controlSpec_dir[numControlSpecies][0],controlSpec_dir[numControlSpecies][1],controlSpec_dir[numControlSpecies][2],
	controlSpec_coord[numControlSpecies][0],controlSpec_coord[numControlSpecies][1] );	
*/
	} 
	

	
}

void AppErbium::add_control_event(int narg, char **arg)
{
    if (narg != 11) error->all("Illegal control_event command: wrong number of arguments");
	int rtype = atoi(arg[0]);
	int iEvent = atoi(arg[1]);
	
    if (rtype < 1 || rtype > 4) {
      char *str = new char[128];
      sprintf(str,"Illegal control_event command: event type %s does not exist",arg[0]);
      error->all(str);
    }
    if ( iEvent < 0 || 
		 (rtype == 1 && iEvent >= none) ||
		 (rtype == 2 && iEvent >= ntwo) ||
 		 (rtype == 3 && iEvent >= nthree) ||
		 (rtype == 4 && iEvent >= nfour) 
		) {
      char *str = new char[128];
      sprintf(str,"Illegal control_event command: event No.%s of type %s does not exist",arg[1],arg[0]);
      error->all(str);
    }
	else {
	
	int n = numControlEvents + 1;
	controlEvent_list = memory->grow_2d_int_array(controlEvent_list,n,2,
					 "canm:controlEvent_list");
	controlEvent_dir = memory->grow_2d_double_array(controlEvent_dir,n,7,
					  "canm:controlEvent_dir");
	controlEvent_coord = memory->grow_2d_int_array(controlEvent_coord,n,2,
					  "canm:controlEvent_coord");
		
	////// add event type
	controlEvent_list[numControlEvents][0] = rtype;
	////// add event id
	controlEvent_list[numControlEvents][1] = iEvent;
		
	///// add the direction of controlled reaction - the orientation vector
    controlEvent_dir[numControlEvents][0] = atof(arg[2]);	//control reaction direction vector - x
    controlEvent_dir[numControlEvents][1] = atof(arg[3]);	//control reaction direction vector - y
    controlEvent_dir[numControlEvents][2] = atof(arg[4]);	//control reaction direction vector - z

	///// add the direction of controlled reaction - the pivot/focus coordinates
    controlEvent_dir[numControlEvents][3] = atof(arg[5]);	//pivot/focus coordinate - x
    controlEvent_dir[numControlEvents][4] = atof(arg[6]);	//pivot/focus coordinate - y
    controlEvent_dir[numControlEvents][5] = atof(arg[7]);	//pivot/focus coordinate - z

    controlEvent_dir[numControlEvents][6] = atof(arg[8]);	//allowed angular deviation +/-0~180 degree from direction vector
    if (controlEvent_dir[numControlEvents][6] < -180 || 
		controlEvent_dir[numControlEvents][6] > 180 )
      error->all("Illegal control_event command: the angular allowance should be +/-0~180 degree");

	//controlled react when the number of neighbors is greater than this lower bound
    controlEvent_coord[numControlEvents][0] = atoi(arg[9]);	
	//control species react when the number of neighbors is less than this upper bound	
    controlEvent_coord[numControlEvents][1] = atoi(arg[10]);			
    if (controlEvent_coord[numControlEvents][0] < 0 || 
		controlEvent_coord[numControlEvents][1] > maxneigh || 
		controlEvent_coord[numControlEvents][0] > controlEvent_coord[numControlEvents][1])
      error->all("Illegal control_event command: the number of neighbors exceeds maximum limit");
/*
printf("\nbefore add_control_event: ");	  
for(int i=0;i<none;i++)
	printf("sctrlIndex[%d]=%d ",i,sctrlIndex[i]);	  
for(int i=0;i<ntwo;i++)
	printf("dctrlIndex[%d]=%d ",i,dctrlIndex[i]);	  
for(int i=0;i<nthree;i++)
	printf("tctrlIndex[%d]=%d ",i,tctrlIndex[i]);	  
for(int i=0;i<nfour;i++)
	printf("qctrlIndex[%d]=%d ",i,qctrlIndex[i]);	  
*/	
	////update the control event index associated with the event
	switch (rtype)
	{
	  case 1: 
		sctrlIndex[iEvent] = numControlEvents;
		break;
	  case 2: 
		dctrlIndex[iEvent] = numControlEvents;
		break;
	  case 3: 
		tctrlIndex[iEvent] = numControlEvents;
		break;
	  case 4: 
		qctrlIndex[iEvent] = numControlEvents;
		break;
	}
/*	
printf("\nafter add_control_event: ");	  
for(int i=0;i<none;i++)
	printf("sctrlIndex[%d]=%d ",i,sctrlIndex[i]);	  
for(int i=0;i<ntwo;i++)
	printf("dctrlIndex[%d]=%d ",i,dctrlIndex[i]);	  
for(int i=0;i<nthree;i++)
	printf("tctrlIndex[%d]=%d ",i,tctrlIndex[i]);	  
for(int i=0;i<nfour;i++)
	printf("qctrlIndex[%d]=%d ",i,qctrlIndex[i]);
*/	
/*	
printf("control_event: type=%d no.%d dir=[%f %f %f] cood=(%d %d)\n",
	controlEvent_list[numControlEvents][0],controlEvent_list[numControlEvents][1],
	controlEvent_dir[numControlEvents][0],controlEvent_dir[numControlEvents][1],controlEvent_dir[numControlEvents][2],
	controlEvent_coord[numControlEvents][0],controlEvent_coord[numControlEvents][1] );	
*/
	numControlEvents++;

	}
}


void AppErbium::generate_control_species_sites(int ispecie)
{
	int n = 0;
	double min_distance = 1e20;
	int mSpec;
//for(int i=0;i<nlocal;i++) printf("site[%d]:iarray=(%d,%d)  ",i,iarray[0][i],iarray[1][i]);	

//printf("generate_control_species_sites():line898:numCtrlSpecSites[%d]=%d \n",ispecie,numCtrlSpecSites[ispecie]);	
	csSiteArray=memory->grow_2d_double_array(csSiteArray,0,2,
					  "canm:csSiteArray");
//printf("generate_control_species_sites():line899:csSiteArray=%d \n",csSiteArray);	
	
	for(int i=0; i<nlocal; i++)
	{
		if(numCtrlSpecSites[ispecie] >= n)
		{	
			n += MEM_DELTA;
			csSiteArray=memory->grow_2d_double_array(csSiteArray,n,2,
					  "canm:csSiteArray");	
//printf("generate_control_species_sites():line910:csSiteArray=%d \n",csSiteArray);	
		}
		if(iarray[1][i]==controlSpec_list[ispecie][0])		//find one site which is the control species
		{
			csSiteArray[numCtrlSpecSites[ispecie]][0]=i;			//store the index of site
			
			//// calculate the keyvalue: here it is the time that the control reaction should occur at the controlled site, 
			//// which is calculated by the dot product between the xyz coodinates and the controlled reaction direction
			csSiteArray[numCtrlSpecSites[ispecie]][1] = 
						xyz[i][0]*controlSpec_dir[ispecie][0]+		
						xyz[i][1]*controlSpec_dir[ispecie][1]+		
						xyz[i][2]*controlSpec_dir[ispecie][2];	
			// find the minimum distance			
			if( min_distance > csSiteArray[numCtrlSpecSites[ispecie]][1])
				min_distance = csSiteArray[numCtrlSpecSites[ispecie]][1];
			numCtrlSpecSites[ispecie]++;

		}
	}

////printf("numCtrlSpecSites[%d]=%d \n",ispecie,numCtrlSpecSites[ispecie]);	
	
//for(int jj=0;jj<numCtrlSpecSites[ispecie];jj++) 
//	printf("csSiteArray[%d]=(%f %f) \n",jj,csSiteArray[jj][0],csSiteArray[jj][1]);

	//update the maximum number of control species sites among all species
	maxNumCtrlSites= (maxNumCtrlSites > numCtrlSpecSites[ispecie] ? maxNumCtrlSites : numCtrlSpecSites[ispecie]);
//printf("generate_control_species_sites():line936:maxNumCtrlSites=%d \n",maxNumCtrlSites);
	

	if(ctrlSpecSites) 
		memory->destroy_1d_T_array(ctrlSpecSites,0);
		
	memory->create_1d_T_array(ctrlSpecSites, 0, numCtrlSpecSites[ispecie]-1,
						     "canm:ctrlSpecSites");

	if(ctrlSpecSites_addr) 
////		memory->destroy_1d_T_array(ctrlSpecSites_addr,0);
////	memory->create_1d_T_array(ctrlSpecSites_addr, 0, numCtrlSpecSites-1,
////						     "canm:ctrlSpecSites_addr");
		memory->sfree(ctrlSpecSites_addr);			
		
	ctrlSpecSites_addr = (csite**) memory->smalloc(numCtrlSpecSites[ispecie]*sizeof(ctrlSpecSites),
					"canm:ctrlSpecSites_addr");						 

	for(int j=0; j<numCtrlSpecSites[ispecie]; j++)
	{
		ctrlSpecSites[j].index = (int)csSiteArray[j][0];
		ctrlSpecSites[j].keyvalue = 
			(csSiteArray[j][1]-min_distance)/controlSpec_rate[ispecie]+ctrlStartTime[ispecie];
		ctrlSpecSites_addr[j] = &ctrlSpecSites[j];	
//printf("ctrlSpecSites_addr[%d].keyvalue=%f\n",j,ctrlSpecSites_addr[j]->keyvalue);
	}
//	memory->destroy_2d_double_array(csSiteArray); //don't need the 2d array any more

	
	qsort(ctrlSpecSites_addr, numCtrlSpecSites[ispecie], sizeof *ctrlSpecSites_addr, compare_site);
////printf("\nafter sorting\n");
////for(int j=0; j<numCtrlSpecSites[ispecie]; j++) printf("ctrlSpecSites_addr[%d].keyvalue=%f\n",j,ctrlSpecSites_addr[j]->keyvalue);

/*

*/
}

void AppErbium::setup_control_species_sites(int ispecies)
{
//printf("setup_control_species_sites %d \n",ispecies);
	for(int jj=0;jj<numCtrlSpecSites[ispecies];jj++)
	{
//printf("ctrlSpecSites_addr[%d]->index=%d ",jj,ctrlSpecSites_addr[jj]->index);
		ctrlSpecSitesList[jj][ispecies].index = ctrlSpecSites_addr[jj]->index;
		ctrlSpecSitesList[jj][ispecies].keyvalue = ctrlSpecSites_addr[jj]->keyvalue;
//printf("ctrlSpecSitesList[%d][%d].keyvalue=%f ",jj,ispecies,ctrlSpecSitesList[jj][ispecies].keyvalue);
	}	

}

int AppErbium::find_species(char *str)
{
  for (int i = 0; i < numspecies; i++)
    if (strcmp(str,sname[i]) == 0) return i;
  return -1;
}

/**************************************************/



/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppErbium::grow_app()
{
  type = iarray[0];
  element = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppErbium::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = new int[3 + 3*maxneigh];
  }

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] < FCC || type[i] > TETRA) flag = 1;
/*************commented by Yan Wang and replaced by the following ********	
    if (element[i] < ERBIUM || element[i] > VACANCY) flag = 1;
***********************************************************************/
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
    if (element[i] < 0 || element[i] >= numspecies) flag = 1;	
/******************************/ 
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all("One or more sites have invalid values");
  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/

////	update_control_species_sites();
 
	ptrCtrlSpec = 0;
	
/******************************/ 
  
}


/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/
 void AppErbium::update_control_species_sites()
 {
	// initialize
	maxNumCtrlSites = 0;
		
	if (numControlSpecies > 0)
	{
	

		for(int ispec=0; ispec<numControlSpecies; ispec++)
		{
			// initialize
			numCtrlSpecSites[ispec] = 0;
			
			generate_control_species_sites(ispec);
//printf("numCtrlSpecSites[%d]=%d\n",ispec,numCtrlSpecSites[ispec]);
//printf("maxNumCtrlSites=%d\n",maxNumCtrlSites);
			memory->grow_2d_T_array(ctrlSpecSitesList, maxNumCtrlSites, numControlSpecies,
							    "canm:ctrlSpecSitesList");	
			setup_control_species_sites(ispec);					 
		}

 	} 
 }
 
/******************************/ 
 
/* ---------------------------------------------------------------------- */

void AppErbium::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates

  if (temperature == 0.0)
    error->all("Temperature cannot be 0.0 for app erbium");

  for (int m = 0; m < none; m++) {
    spropensity[m] = srate[m];
    scount[m] = 0;
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = exp(-drate[m]/temperature);
    dcount[m] = 0;
  }
  for (int m = 0; m < nthree; m++) {
    tpropensity[m] = exp(-trate[m]/temperature);
    tcount[m] = 0;
  }
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  

 for (int m = 0; m < nfour; m++) {
    qpropensity[m] = exp(-qrate[m]/temperature);
    qcount[m] = 0;
  }
//printf("setup_app():line1065 \n");    
/*****************************/ 
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppErbium::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppErbium::site_propensity(int i)
{
  int j,k,m;

  // valid single, double, triple events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;

  // currently fcc sites have no events
  // can remove this later
/******* commented by Yan Wang, July 2010
  if (type[i] == FCC) return proball = 0.0;
*******/
  // single-site events
  for (m = 0; m < none; m++) {
    if (type[i] != stype[m] || element[i] != sinput[m]) continue;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/
/*****************************/ 
	add_event(i,1,m,spropensity[m],-1,-1);
    proball += spropensity[m];
  }

  // double-site events
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < ntwo; m++) {
      if (type[i] != dtype[m][0] || element[i] != dinput[m][0]) continue;
      if (type[j] != dtype[m][1] || element[j] != dinput[m][1]) continue;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/
	if( dctrlIndex[m]!=-1 )
	{
	
	
		//// the direction vector, right now, is determinstric as the specified mean values. Could be changed to random later on
		double ax=controlEvent_dir[dctrlIndex[m]][0];
		double ay=controlEvent_dir[dctrlIndex[m]][1];
		double az=controlEvent_dir[dctrlIndex[m]][2];
				
		if(ax*ax+ay*ay+az*az < EPSILON)	// if the directional vector is zero, should follow direction to the pivot/focus point
		{
			//// the specified direction should be from the current site i to the pivot/focus point -- "converge"
			ax=controlEvent_dir[dctrlIndex[m]][3] - xyz[i][0];
			ay=controlEvent_dir[dctrlIndex[m]][4] - xyz[i][1];
			az=controlEvent_dir[dctrlIndex[m]][5] - xyz[i][2];
			
			//// the specified direction should be from the pivot/focus point to the current site i instead -- "diverge"
			if(controlEvent_dir[dctrlIndex[m]][6]<0)
			{
				ax=-ax;
				ay=-ay;
				az=-az;
			}
		}

		//// the projected path from the current site i to the first neighboring site j
		double bx=xyz[j][0] - xyz[i][0];
		double by=xyz[j][1] - xyz[i][1];
		double bz=xyz[j][2] - xyz[i][2];
		
		double d1=ay*bz - az*by;
		double d2=az*bx - ax*bz;
		double d3=ax*by - ay*bx;
			
		double sinval = sqrt( (d1*d1+d2*d2+d3*d3)/((bx*bx+by*by+bz*bz)*(ax*ax+ay*ay+az*az)) );
		double dotprod = ax*bx+ay*by+az*bz;		//dot product: positive if angle<90; negative if angle>90
//printf("\n[%f %f %f][%f %f %f] sin=%f-%f dotp=%f", ax,ay,az, bx,by,bz, 
//	sinval,sin(controlEvent_dir[dctrlIndex[m]][6]*PI/180),dotprod); 

		if( fabs(controlEvent_dir[dctrlIndex[m]][6])<90 && 
			dotprod <= 0 											|| //specified<90 but calculated>90
			fabs(controlEvent_dir[dctrlIndex[m]][6])<90 && dotprod>=0 &&
			sinval-sin(fabs(controlEvent_dir[dctrlIndex[m]][6])*PI/180) > EPSILON || //specified<90,calculated<90,but>specified
			fabs(controlEvent_dir[dctrlIndex[m]][6])>=90 && dotprod<0 &&
			sin(fabs(controlEvent_dir[dctrlIndex[m]][6])*PI/180)-sinval >	EPSILON  //specified>90,calculated>90,but>specified
		  )
		{ 
//printf("-skipped ");
			continue;			
		}
	}
/*****************************/ 
      add_event(i,2,m,dpropensity[m],j,-1);
      proball += dpropensity[m];
//printf("site_propensity():line 1155: site %d (type %d) neighb %d event added proball=%f\n",i,m,j,proball);		
    }
  }

  // triple-site events
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (int kk = 0; kk < numneigh[i]; kk++) {	
      if (jj == kk) continue;
      k = neighbor[i][kk];
      for (m = 0; m < nthree; m++) {
	if (type[i] != ttype[m][0] || element[i] != tinput[m][0]) continue;
	if (type[j] != ttype[m][1] || element[j] != tinput[m][1]) continue;
	if (type[k] != ttype[m][2] || element[k] != tinput[m][2]) continue;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/
	if( tctrlIndex[m]!=-1 )
	{
		//// the direction vector, right now, is determinstric as the specified mean values. Could be changed to random later on
		double ax=controlEvent_dir[tctrlIndex[m]][0];
		double ay=controlEvent_dir[tctrlIndex[m]][1];
		double az=controlEvent_dir[tctrlIndex[m]][2];

		if(ax*ax+ay*ay+az*az < EPSILON)	// if the directional vector is zero, should follow direction to the pivot/focus point
		{
			//// the specified direction should be from the current site i to the pivot/focus point instead
			ax=controlEvent_dir[dctrlIndex[m]][3] - xyz[i][0];
			ay=controlEvent_dir[dctrlIndex[m]][4] - xyz[i][1];
			az=controlEvent_dir[dctrlIndex[m]][5] - xyz[i][2];
			
			//// the specified direction should be from the pivot/focus point to the current site i instead -- "diverge"
			if(controlEvent_dir[dctrlIndex[m]][6]<0)
			{
				ax=-ax;
				ay=-ay;
				az=-az;
			}
		}
		
		//// the projected path from the current site i to the first neighboring site j
		double bx=xyz[j][0] - xyz[i][0];
		double by=xyz[j][1] - xyz[i][1];
		double bz=xyz[j][2] - xyz[i][2];
		
		double d1=ay*bz - az*by;
		double d2=az*bx - ax*bz;
		double d3=ax*by - ay*bx;
		
		double sinval = sqrt( (d1*d1+d2*d2+d3*d3)/((bx*bx+by*by+bz*bz)*(ax*ax+ay*ay+az*az)) );
		double dotprod = ax*bx+ay*by+az*bz;		//dot product: positive if angle<90; negative if angle>90
//printf("\n[%f %f %f][%f %f %f] sin=%f-%f dotp=%f", ax,ay,az, bx,by,bz, 
//	sinval,sin(controlEvent_dir[dctrlIndex[m]][6]*PI/180),dotprod); 

		if( fabs(controlEvent_dir[dctrlIndex[m]][6])<90 && 
			dotprod <= 0 											|| //specified<90 but calculated>90
			fabs(controlEvent_dir[dctrlIndex[m]][6])<90 && dotprod>=0 &&
			sinval-sin(fabs(controlEvent_dir[dctrlIndex[m]][6])*PI/180) > EPSILON || //specified<90,calculated<90,but>specified
			fabs(controlEvent_dir[dctrlIndex[m]][6])>=90 && dotprod<0 &&
			sin(fabs(controlEvent_dir[dctrlIndex[m]][6])*PI/180)-sinval >	EPSILON  //specified>90,calculated>90,but>specified
		  )
		{
			continue;			
		}
	}
/*****************************/ 
	add_event(i,3,m,tpropensity[m],j,k);
	proball += tpropensity[m];
//printf("site_propensity():line 1207: site %d neighbor %d event added proball=%f\n",i,m,proball);		
      }
    }
  }

 
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/  
  int k3;
  // quadruple-site events
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (int kk = 0; kk < numneigh[i]; kk++) {	
      if (jj == kk) continue;
	  for (int kk3 = 0; kk3 < numneigh[i]; kk3++) {	
		if (jj == kk3 || kk == kk3) continue;
		    
			k3 = neighbor[i][kk3];
		    for (m = 0; m < nfour; m++) {
			if (type[i] != qtype[m][0] || element[i] != qinput[m][0]) continue;
			if (type[j] != qtype[m][1] || element[j] != qinput[m][1]) continue;
			if (type[k] != qtype[m][2] || element[k] != qinput[m][2]) continue;
			if (type[k3] != qtype[m][3] || element[k3] != qinput[m][3]) continue;

	if( qctrlIndex[m]!=-1 )
	{
		//// the direction vector, right now, is determinstric as the specified mean values. Could be changed to random later on
		double ax=controlEvent_dir[qctrlIndex[m]][0];
		double ay=controlEvent_dir[qctrlIndex[m]][1];
		double az=controlEvent_dir[qctrlIndex[m]][2];

		if(ax*ax+ay*ay+az*az < EPSILON)	// if the directional vector is zero, should follow direction to the pivot/focus point
		{
			//// the specified direction should be from the current site i to the pivot/focus point instead
			ax=controlEvent_dir[dctrlIndex[m]][3] - xyz[i][0];
			ay=controlEvent_dir[dctrlIndex[m]][4] - xyz[i][1];
			az=controlEvent_dir[dctrlIndex[m]][5] - xyz[i][2];
			
			//// the specified direction should be from the pivot/focus point to the current site i instead -- "diverge"
			if(controlEvent_dir[dctrlIndex[m]][6]<0)
			{
				ax=-ax;
				ay=-ay;
				az=-az;
			}			
		}
		
		//// the projected path from the current site i to the first neighboring site j
		double bx=xyz[j][0] - xyz[i][0];
		double by=xyz[j][1] - xyz[i][1];
		double bz=xyz[j][2] - xyz[i][2];
		
		double d1=ay*bz - az*by;
		double d2=az*bx - ax*bz;
		double d3=ax*by - ay*bx;

		double sinval = sqrt( (d1*d1+d2*d2+d3*d3)/((bx*bx+by*by+bz*bz)*(ax*ax+ay*ay+az*az)) );
		double dotprod = ax*bx+ay*by+az*bz;		//dot product: positive if angle<90; negative if angle>90
//printf("\n[%f %f %f][%f %f %f] sin=%f-%f dotp=%f", ax,ay,az, bx,by,bz, 
//	sinval,sin(controlEvent_dir[dctrlIndex[m]][6]*PI/180),dotprod); 

		if( fabs(controlEvent_dir[dctrlIndex[m]][6])<90 && 
			dotprod <= 0 											|| //specified<90 but calculated>90
			fabs(controlEvent_dir[dctrlIndex[m]][6])<90 && dotprod>=0 &&
			sinval-sin(fabs(controlEvent_dir[dctrlIndex[m]][6])*PI/180) > EPSILON || //specified<90,calculated<90,but>specified
			fabs(controlEvent_dir[dctrlIndex[m]][6])>=90 && dotprod<0 &&
			sin(fabs(controlEvent_dir[dctrlIndex[m]][6])*PI/180)-sinval >	EPSILON  //specified>90,calculated>90,but>specified
		  )
		{
			continue;			
		}

	}		
			add_event(i,4,m,qpropensity[m],j,k,k3);
			proball += qpropensity[m];
			}
	  
	  }
    }
  }

/*******************************/	
 
  return proball;
}

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
IntervalItem AppErbium::site_propensity_intval(int i)
{
	IntervalItem a = {0,0.0,0.0,0.0};
	return a;
}
/*****************************/
/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppErbium::site_event(int i, class RandomPark *random)
{
  int j,k,m,n;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform single, double, or triple event

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;
  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
 int k3 = events[ievent].k3partner;
/*********************************/

  int ibef = 0;
  for (int jj = 0; jj < nlocal; jj++)
    if (element[jj] == HYDROGEN) ibef++;
  int iel1 = element[i];
  int iel2 = element[j];

  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
  } else if (rstyle == 2) {
    element[i] = doutput[which][0];
    element[j] = doutput[which][1];
    dcount[which]++;
  } else if (rstyle == 3) {
    element[i] = toutput[which][0];
    element[j] = toutput[which][1];
    element[k] = toutput[which][2];
    tcount[which]++;
  }
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
 else if (rstyle == 4) {
    element[i] = qoutput[which][0];
    element[j] = qoutput[which][1];
    element[k] = qoutput[which][2];
    element[k3] = qoutput[which][3];
    qcount[which]++;
	}
/*******************************/	
  // compute propensity changes for participating sites and first neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;
  int isite = i2site[i];
//printf("site_event():line 1359\n");  
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if (rstyle >= 2) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(m);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  if (rstyle >= 3) {
    for (n = 0; n < numneigh[k]; n++) {
      m = neighbor[k][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(m);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }
  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
  if (rstyle >= 4) {
    for (n = 0; n < numneigh[k3]; n++) {
      m = neighbor[k3][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(m);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

/*****************************/ 
  solve->update(nsites,esites,propensity);

  // clear echeck array
  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;

  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (July  2010)
 ******************************/ 
	
///// check if the ctrlSpecSites need to be updated	
	int updateNeeded = 0;
	while (ptrCtrlSpec < numControlSpecies && time > ctrlStartTime[ptrCtrlSpec])
	{
		ptrCtrlSpec++;
//printf("ptrCtrlSpec=%d\n",ptrCtrlSpec);
		updateNeeded = 1;
	}
	if (updateNeeded)
	{
//printf("update when ptrCtrlSpec=%d\n",ptrCtrlSpec);
		update_control_species_sites();
	}	
///// If there is any event at the site of control species need to be fired (which is
///// tracked by the control event clock time), fire site_event again at those sites

	int c;
	for(int iSpec=0; iSpec<numControlSpecies; iSpec++)
	{
		while(	(currCtrlSpecSite[iSpec] < numCtrlSpecSites[iSpec]) &&
				(ctrlSpecSitesList[currCtrlSpecSite[iSpec]][iSpec].keyvalue < time) )
		{
//printf("\ncurrCtrlSpecSite[%d]=%d ",iSpec,ctrlSpecSites_addr[currCtrlSpecSite[iSpec]]->index);
			///// Find the site where control reaction occurs
			c = ctrlSpecSitesList[currCtrlSpecSite[iSpec]][iSpec].index;
			///// Update its associated element to the target specie
			
//printf(" befor:element[%d]=%d ",c,element[c]);
			element[c] = controlSpec_list[iSpec][1];
			
//printf(" after:element[%d]=%d \n",c,element[c]);


/*			///// Update the propensity of the site as well as its neighbors
			  isite = i2site[c];
			//  isite = c;
			  propensity[isite] = site_propensity(c);
//			  esites[nsites++] = isite;
			  echeck[isite] = 1;
//printf("isite=%d \n",isite);
			  for (n = 0; n < numneigh[c]; n++) {
			    m = neighbor[c][n];
			    isite = i2site[m];
//printf("isite=%d ",isite);
			    if (isite >= 0 && echeck[isite] == 0) {
			      propensity[isite] = site_propensity(m);
//			      esites[nsites++] = isite;
			      echeck[isite] = 1;
			    }
			  }
*/	
			///// Update the pointer to the current control reaction site 
			currCtrlSpecSite[iSpec]++;
		}
	}
/*******************************/ 
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppErbium::clear_events(int i)
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

void AppErbium::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner)
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

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
  events[freeevent].k3partner = -1;
/******************************/
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}



/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  

void AppErbium::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner, int k3partner)
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

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].k3partner = k3partner;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/********************************/

/* ----------------------------------------------------------------------
   grow list of stored reactions for single, double, or triple
------------------------------------------------------------------------- */

void AppErbium::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/erbium:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/erbium:spropensity");
    stype = (int *) 
      memory->srealloc(stype,n*sizeof(int),"app/erbium:stype");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/erbium:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/erbium:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/erbium:scount");
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/
    sctrlIndex = (int *) 
      memory->srealloc(sctrlIndex,n*sizeof(int),"app/canm:sctrlIndex"); 
	sctrlIndex[none] = -1;	//by default not a controlled event  
/*****************************/ 
  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/erbium:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/erbium:dpropensity");
    dtype = memory->grow_2d_int_array(dtype,n,2,"app/erbium:dtype");
    dinput = memory->grow_2d_int_array(dinput,n,2,"app/erbium:dinput");
    doutput = memory->grow_2d_int_array(doutput,n,2,"app/erbium:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/erbium:dcount");
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
    dctrlIndex = (int *) 
      memory->srealloc(dctrlIndex,n*sizeof(int),"app/canm:dctrlIndex"); 
	dctrlIndex[ntwo] = -1;	//by default not a controlled event    
/*****************************/ 
  } else if (rstyle == 3) {
    int n = nthree + 1;
    trate = (double *) 
      memory->srealloc(trate,n*sizeof(double),"app/erbium:trate");
    tpropensity = (double *) 
      memory->srealloc(tpropensity,n*sizeof(double),"app/erbium:tpropensity");
    ttype = memory->grow_2d_int_array(ttype,n,3,"app/erbium:ttype");
    tinput = memory->grow_2d_int_array(tinput,n,3,"app/erbium:tinput");
    toutput = memory->grow_2d_int_array(toutput,n,3,"app/erbium:toutput");
    tcount = (int *) 
      memory->srealloc(tcount,n*sizeof(int),"app/erbium:tcount");
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
    tctrlIndex = (int *) 
      memory->srealloc(tctrlIndex,n*sizeof(int),"app/canm:tctrlIndex"); 
	tctrlIndex[nthree] = -1;	//by default not a controlled event  
/*****************************/ 
  }
  
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
	else if (rstyle == 4) {
    int n = nfour + 1;
    qrate = (double *) 
      memory->srealloc(qrate,n*sizeof(double),"app/erbium:qrate");
    qpropensity = (double *) 
      memory->srealloc(qpropensity,n*sizeof(double),"app/erbium:qpropensity");
    qtype = memory->grow_2d_int_array(qtype,n,4,"app/erbium:qtype");
    qinput = memory->grow_2d_int_array(qinput,n,4,"app/erbium:qinput");
    qoutput = memory->grow_2d_int_array(qoutput,n,4,"app/erbium:qoutput");
    qcount = (int *) 
      memory->srealloc(qcount,n*sizeof(int),"app/erbium:qcount");
    qctrlIndex = (int *) 
      memory->srealloc(qctrlIndex,n*sizeof(int),"app/canm:qctrlIndex"); 
	qctrlIndex[nfour] = -1;	//by default not a controlled event  
  }
  
  
 /********************************/
}
