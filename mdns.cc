#include "options.h"
#include "dns.h"
#include "MultiIntegrator.h"

const double ProblemVersion=1.0;

#ifndef DEPEND
#if !COMPLEX
#error Navier requires COMPLEX=1 in options.h
#endif
#endif

const char *problem="Multispectral Direct Numerical Simulation of Turbulence";

const char *method="MDNS";
const char *integrator="RK5";
const char *ic="Equipartition";
//const char *linearity="Power";
const char *forcing="WhiteNoiseBanded";

// Vocabulary
Real nuH=0.0, nuL=0.0;
int pH=1;
int pL=0;
unsigned Nx=15;
unsigned Ny=15;
Real eta=0.0;
Complex force=0.0;
Real kforce=1.0;
Real deltaf=1.0;
unsigned movie=0;
unsigned rezero=0;
unsigned spectrum=1;
Real icalpha=1.0;
Real icbeta=1.0;

// global variables for dns.h
int xpad=1;
int ypad=1;

// DNS setup routines
DNS::DNS()
{
  cout << "DNS::DNS()" << endl;
  //  exit(1);
}

DNS::~DNS()
{
}

void DNS::InitialConditions()
{
  cout << "DNS::InitialConditions()" << endl;
  exit(1);
}


InitialConditionBase *InitialCondition;
ForcingBase *Forcing;

//Global variables for MultiIntegrator.h
MultiProblem *MProblem;
unsigned Ngrids=1;
const char *subintegrator; 

class MDNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "MDNS";}
  const char *Abbrev() {return "MDNS";}
  MDNSVocabulary();

  Table<InitialConditionBase> *InitialConditionTable;
  InitialConditionBase *NewInitialCondition(const char *& key) {
    return InitialConditionTable->Locate(key);
  }

  Table<ForcingBase> *ForcingTable;
  ForcingBase *NewForcing(const char *& key) {
    return ForcingTable->Locate(key);
  }
};

MDNSVocabulary MDNS_Vocabulary;

// define the Grid class:
class Grid {
  // FIXME: flesh this out
};

//class MDNS : public DNS, public MultiProblem {
//class MDNS : public DNS {
class MDNS : public MultiProblem {
private:
  unsigned mx, my; // size of data arrays
  array1<Grid *> G;
public:
  MDNS();
  ~MDNS();
  //  enum Field {OMEGA,TRANSFER,EK,Nfields};
  enum Field {OMEGA,Nfields};

  Table<InitialConditionBase> *InitialConditionTable;
  //void Initialize();
  //void InitialConditions(unsigned Ngrids0);
  void InitialConditions();

  void Project(unsigned g);
  void Prolong(unsigned g);
  
  void Output(int it);
  virtual void Source(const vector2&, const vector2&, double);


};

// Global problem.
MDNS *MDNSProblem;


/****** Vocabulary *****/

MDNSVocabulary::MDNSVocabulary()
{
  cout << "MDNSVocabulary::MDNSVocabulary()"<<endl;
  //exit(1);
  Vocabulary=this;

  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  //  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  // VOCAB(spectrum,0,1,"Spectrum flag (0=off, 1=on)");
  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");

  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  VOCAB(icalpha,0.0,0.0,"initial condition parameter");
  VOCAB(icbeta,0.0,0.0,"initial condition parameter");
  INITIALCONDITION(Zero);
  //  INITIALCONDITION(Constant);
  //  INITIALCONDITION(Equipartition); 

  VOCAB(nuH,0.0,REAL_MAX,"High-wavenumber viscosity");
  VOCAB(nuL,0.0,REAL_MAX,"Low-wavenumber viscosity");
  VOCAB(pH,0,0,"Power of Laplacian for high-wavenumber viscosity");
  VOCAB(pL,0,0,"Power of Laplacian for molecular viscosity");

  VOCAB_NOLIMIT(forcing,"Forcing type");
  ForcingTable=new Table<ForcingBase>("forcing");

  VOCAB(eta,0.0,REAL_MAX,"vorticity injection rate");
  VOCAB(force,(Complex) 0.0, (Complex) 0.0,"constant external force");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");
  //FORCING(None);
  //FORCING(WhiteNoiseBanded);

  METHOD(MDNS); 

  subintegrator="rk5";
  VOCAB_NOLIMIT(subintegrator,"subintegrator for multi-integrator");
  INTEGRATOR(MultiIntegrator);
  VOCAB(Ngrids,1,INT_MAX,"Number of multispectral grids");
  cout << "... done MDNSVocabulary::MDNSVocabulary()"<<endl;
}


MDNS::MDNS() 
{
  MDNSProblem=this;
  MProblem=this;
  check_compatibility(DEBUG);
}

MDNS::~MDNS() 
{
}

//void MDNS::InitialConditions(unsigned Ngrids0) {
void MDNS::InitialConditions()
{
  cout << "MDNS::InitialConditions()" << endl;
  nfields=Nfields;
  // Vocabulary
  Ngrids=::Ngrids;
  MultiProblem::InitialConditions(Ngrids);

  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  for(unsigned i=0; i < Ngrids; ++i) {
    NY[Nfields*i+OMEGA]=Nx*my;
    //NY[Nfields*i+TRANSFER]=
    //NY[Nfields*i+MOMENT]=
  }
  
  Allocator(); // FIXME: what is this?

  G.Allocate(Ngrids);
  for(unsigned int i=0; i < Ngrids; ++i) {
    vector y,m,T;
    Dimension(y,Y[Nfields*i+OMEGA]);
    //Dimension(T,Y[Nfields*i+TRANSFER]);
    //Dimension(m,Y[Nfields*i+MOMENT]);
    G[i]=new Grid; // FIXME: must pass data
  }
  
  //for(unsigned int i=0; i < Ngrids; ++i) G[i]->initialspectrum();
  
}

void MDNS::Project(unsigned g) 
{
  // FIXME
}
void MDNS::Prolong(unsigned g) 
{
  // FIXME
}


void MDNS::Output(int it)
{
  // FIXME
}
void DNS::Output(int it)
{
  // FIXME
  // This is just to get some random table thing working.
}

void MDNS::Source(const vector2&, const vector2&, double)
{
  // FIXME
}


