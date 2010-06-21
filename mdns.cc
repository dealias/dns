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

//***** Initial Conditions *****//
InitialConditionBase *InitialCondition;
class Constant : public InitialConditionBase {
public:
  const char *Name() {return "Constant";}
  void Set(Complex *w, unsigned n) {
    for(unsigned i=0; i < n; i++) {
      w[i]=Complex(icalpha,icbeta);
    }
  }
};

//***** Forcing *****//
ForcingBase *Forcing;

//***** Global variables for MultiIntegrator.h *****//
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

// The Grid class contains all the information for calculating source
// terms on a subgrid.
class Grid : public DNS {
private:
  //MDNS *parent;
  //protected:
public:
  Grid();
  //  Grid() {}
  ~Grid();
  void InitialConditions();

  void  NonLinearSource(const vector2& Src, const vector2& Y, double t);
  
  // FIXME: flesh this out
};


//class MDNS : public DNS, public MultiProblem {
//class MDNS : public DNS {
class MDNS : public MultiProblem {
private:
  //  array2<Complex> w; // Vorticity field
  //array2<Complex> s;
  //vector2 w,s;
  unsigned mx, my; // size of data arrays
  array1<Grid *> G;
public:
  friend class Grid;
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

  // Source functions
  void Source(const vector2&, const vector2&, double);
  void ConservativeSource(const vector2& Src, const vector2& Y, double t);
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void ExponentialSource(const vector2& Src, const vector2& Y, double t);
};

// Global problem.
MDNS *MDNSProblem;

Grid::Grid()
{
  cout << MDNSProblem->nfields << endl;
  cout << MDNSProblem->mx;
  Nx=2;
  cout << "grid " << endl;
}

Grid::~Grid()
{
  cout << "no grid!"<<endl;
}

void Grid::InitialConditions()
{
    // load vocabulary from global variables
  xpad=::xpad;
  ypad=::ypad;
  nuH=::nuH;
  nuL=::nuL;
  pH=::pH;
  pL=::pL;
  Nx=::Nx;
  Ny=::Ny;
  //  spectrum=::spectrum;
  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");

  k0=1.0;
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  xorigin=mx-1;

  NY[OMEGA]=Nx*my;
  NY[TRANSFER]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;

  cout << "\nALLOCATING FFT BUFFERS" << endl;
  align=sizeof(Complex);

  Allocator(align);  // FIXME: refers to kernel.h, but should go somewhere else?

  //Dimension(T,nshells);

  // TODO: block could be used as work arrays for more than one grid

  
  unsigned int Nxmy=Nx*my;
  unsigned int nbuf=3*Nxmy;
  unsigned int Nx0=Nx+xpad;
  unsigned int Ny0=Ny+ypad;
  int my0=Ny0/2+1;
  xorigin=mx-1;

  w.Dimension(Nx,my);
  f0.Dimension(Nx,my);
  block=fftwpp::ComplexAlign(nbuf);
  f1.Dimension(Nx,my,block);
  g0.Dimension(Nx,my,block+Nxmy);
  g1.Dimension(Nx,my,block+2*Nxmy);
  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,2);
  F[1]=f1;
  G[0]=g0;
  G[1]=g1;



  // FIXME: ic=<whatever> at command line hangs program, but fine if in p file.
  InitialCondition=MDNS_Vocabulary.NewInitialCondition(ic);
  w.Set(Y[OMEGA]);
  InitialCondition->Set(w,NY[OMEGA]);
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);

  for(unsigned i=0; i < nshells; i++)
    Y[EK][i]=0.0;

  // TODO: more initialization, etc!
}

void Grid::NonLinearSource(const vector2& Src, const vector2& Y, double t)
{
  DNS::NonLinearSource(Src,Y,t);
}

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
  INITIALCONDITION(Constant);
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
  align=sizeof(Complex);
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
  
  Allocator(align); // FIXME: what is this?

  G.Allocate(Ngrids);
  for(unsigned int i=0; i < Ngrids; ++i) {
    vector y,m,T;
    Dimension(y,Y[Nfields*i+OMEGA]);
    //Dimension(T,Y[Nfields*i+TRANSFER]);
    //Dimension(m,Y[Nfields*i+MOMENT]);
    G[i]=new Grid; // FIXME: must pass data as well
  }
  for(unsigned int i=0; i < Ngrids; ++i) 
    G[i]->InitialConditions();
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
  // This is just to get some random table thing working.
}

void MDNS::Source(const vector2& Src, const vector2& Y, double t)
{
  cout << "Grid " << grid <<endl;
  ConservativeSource(Src,Y,t);
  //NonConservativeSource(Src,Y,t);
  // FIXME
  // Calculate source routines for all of the grids.
}

void MDNS::ConservativeSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  //Transfer(Src,Y);
  //LinearSource(Src,Y,t);
}

void MDNS::NonConservativeSource(const vector2& Src, const vector2& Y, double t)
{
  //Moments(Src,Y,t);
}

void MDNS::LinearSource(const vector2& Src, const vector2& Y, double t) 
{
  //vector y,source;
  //Dimension(y,Y[OMEGA]);
  //Dimension(source,Src[OMEGA]);
  //G[grid]->LinearSource(source,y,t);
}

void MDNS::NonLinearSource(const vector2& Src, const vector2& Y, double) 
{
  G[grid]->NonLinearSource(Src,Y,t);
}

void MDNS::ExponentialSource(const vector2& Src, const vector2& Y, double t) 
{
  //NonLinearSource(Src,Y,t);
  //Transfer(Src,Y);
  //NonConservativeSource(Src,Y,t);
}

//void MDNS::Stochastic(const vector2& Y, double t, double dt) 
//{
  /*
  vector y,T;
  Dimension(y,Y[OMEGA]);
  Dimension(T,Y[TRANSFER]);
  G[grid]->Stochastic(T,y,t,dt);
  Sigma();
  */
//}

//void MDNS::Moments(const vector2& Src, const vector2& Y, double t) 
//{
  /*
  vector y;
  array2<Var> moment;
  Dimension(y,Y[OMEGA]);
  moment.Dimension(Nmoments,Nshells,Src[MOMENT]);
  G[grid]->Moments(moment,y,t);
  */
//}
  
//void MDNS::Transfer(const vector2& Src, const vector2& Y) 
//{
  /*
  vector y,source,T;
  Dimension(y,Y[OMEGA]);
  Dimension(source,Src[OMEGA]);
  Dimension(T,Src[TRANSFER]);
  G[grid]->Transfer(T,source,y);
  */
//}


