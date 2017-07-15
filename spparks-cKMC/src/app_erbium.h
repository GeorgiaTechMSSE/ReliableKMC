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

#ifdef APP_CLASS
AppStyle(erbium,AppErbium)

#else

#ifndef SPK_APP_ERBIUM_H
#define SPK_APP_ERBIUM_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppErbium : public AppLattice {
  friend class DiagErbium;

 public:
  AppErbium(class SPPARKS *, int, char **);
  ~AppErbium();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/  
  IntervalItem site_propensity_intval(int);
/*****************************/
 private:
  int engstyle;
  int *type,*element;      // variables on each lattice site
  int firsttime;

  int *esites;
  int *echeck;

  int none,ntwo,nthree;
  double *srate,*drate,*trate;
  double *spropensity,*dpropensity,*tpropensity;
  int *stype,**dtype,**ttype;
  int *sinput,**dinput,**tinput;
  int *soutput,**doutput,**toutput;
  int *scount,*dcount,*tcount;
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
  int nfour;			//parameters for events with 4 (quadruple) neighbors involved
  double *qrate;
  double *qpropensity;
  int **qtype;
  int **qinput;
  int **qoutput;
  int *qcount;
  
  //1D array storing the index of control events corresponding to "controlEvent_list"
  // store -1 if not controlled event
  int *sctrlIndex, *dctrlIndex, *tctrlIndex, *qctrlIndex;

  int numspecies;                  // # of unique species
  char **sname;                  // ID of each species
  int *pcount;               // counts for each species

  
  int numControlSpecies;			//# of control species
  int** controlSpec_list;		//list of the control species
  double* controlSpec_rate;		//control species reaction rate
  double** controlSpec_dir;		//control species reaction direction vector
  int** controlSpec_coord;		////lower and upper bounds of neighbors that control species would react 

  int numControlEvents;			//# of controlled events
  int** controlEvent_list;		//list of the indices for controlled events 
								// numControlEventsX2 array, 1st column: type of event; 2nd column: index
  double** controlEvent_dir;	// numControlEventsX7 array, controlled reaction direction vector
								// first 3 columns: mean values of x,y,z directions, 
								// next 3 columns: standard deviations of direction vectors 
								// last column: the allowed angular deviation from the above vector
  int** controlEvent_coord;		////numControlEventsX2 array, lower and upper bounds of neighbors that controlled event would react 

  double *ctrlStartTime;			//the time when the controll reaction start
////  double *checkPointCtrlSpec;	//the list of time to update ctrlSpecSites
  int	ptrCtrlSpec;			//the index pointing to the index of *ctrlStartTime when should update ctrlSpecSites
  int *numCtrlSpecSites;			//# of sites where control species reside
  int *currCtrlSpecSite;		//keep track of current index of control sites as the time advances
  csite *ctrlSpecSites;			//list of control species sites used in generating the list
  csite **ctrlSpecSites_addr;	//the address of the sorted control species sites used in generating the list
  csite **ctrlSpecSitesList;	//the actual 2D list of control species sites 
  int  maxNumCtrlSites;			//keep track of the maxium number of control species sites
  double **csSiteArray;		//a 2d double array only internally used in generating the list of control specie sites

  
/******************************/

  struct Event {           // one event for an owned site
    int style;             // reaction style = SINGLE,DOUBLE,TRIPLE
    int which;             // which reaction of this type
    int jpartner,kpartner; // which J,K neighbors of I are part of event

/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
	int k3partner;			//the third neighbor
/*******************************/	
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, int, double, int, int);
/****************************
 * Inserted for interval based KMC 
 *     - By Yan Wang (June  2010)
 ******************************/ 
  void add_event(int, int, int, double, int, int, int);	//overloaded function
  void add_species(int, char **);		//add species of elements
  void add_control_species(int, char **);		//add control species that reations are controllable
  void add_control_event(int, char **);		//add events that are controllable
  int find_species(char *);
  void parse_event(int, char **);
  void generate_control_species_sites(int);
  void setup_control_species_sites(int);
  void update_control_species_sites();
/********************************/  
  void grow_reactions(int);
};

}

#endif
#endif
