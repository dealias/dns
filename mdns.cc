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
unsigned radix=4;
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
}

DNS::~DNS()
{
}

void DNS::InitialConditions()
{
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

// ***** Vocabulary ***** //
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


//***** Base class for functionoids *****//
class FunctRR {
 public:
  virtual void f(const Real w2, const Real k2, int n,va_list args) = 0;
 };
typedef FunctRR* FunctRRPtr;

unsigned gN(unsigned N, unsigned g) {
  return radix == 1 ? pow(2,(int) g)*(N+1)-1 : N;
};
unsigned gm(unsigned m, unsigned g) {
  return radix == 1 ? pow(2,(int) g)*m : m;
};


//***** Grid class based on DNS *****//
class Grid : public DNS {
private:
  unsigned Invisible;
  int Invisible2;
  unsigned myg;

public:
  Grid();
  ~Grid();

  void InitialConditions(unsigned g);
  Real gk(Real k, unsigned g) {return k*pow(sqrt((double) radix),g);};

  void NonLinearSource(const vector2& Src, const vector2& Y, double t);

  void ComputeInvariants(Real& E, Real& Z, Real& P);

  void loopwF(const FunctRRPtr,int n,...);
  class Invariants : public FunctRR {
    void f(const Real w2, const Real k02, int n,va_list args) {
      double * Z=va_arg(args,double *);
      double * E=va_arg(args,double *);
      double * P=va_arg(args,double *);
      *Z += w2;
      *E += w2/k02;
      *P += w2*k02;
    }
  };
};

void Grid::loopwF(const FunctRRPtr F,int n,...)
{
  va_list args;
  va_start(args,n);
  w.Set(Y[OMEGA]);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j > Invisible || I2 > Invisible2) {
	Real w2=abs2(wi[j]);
	Real k2=k02*(I2+j*j);
	F->f(w2,k2,n,args);
      }
    }
  }
  va_end(args);
};

//***** MDSN derived from MultiProblem *****//
class MDNS : public MultiProblem {
private:
  unsigned mx, my; // size of data arrays
  array1<Grid *> G;
  void ComputeInvariants(Real& E, Real& Z, Real& P);
public:
  friend class Grid;
  MDNS();
  ~MDNS();
  enum Field {OMEGA,Nfields}; // TODO: include TRANSFER and EK as well


  Table<InitialConditionBase> *InitialConditionTable;
  void InitialConditions();

  // Source functions
  void Source(const vector2&, const vector2&, double);
  void ConservativeSource(const vector2& Src, const vector2& Y, double t);
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void ExponentialSource(const vector2& Src, const vector2& Y, double t);
  void Stochastic(const vector2& Y, double t, double dt);
  
  void Project(unsigned g);
  void Prolong(unsigned g);
  
  void FinalOutput();
  void Output(int it);
};

//***** Global problem *****//
MDNS *MDNSProblem;

// ***** Grid functions *****//
Grid::Grid()
{
}

Grid::~Grid()
{
  fftwpp::deleteAlign(block);
}

void Grid::InitialConditions(unsigned g)
{
  // load vocabulary from global variables
  myg=g;
  xpad=::xpad;
  ypad=::ypad;
  nuH=::nuH;
  nuL=::nuL;
  pH=::pH;
  pL=::pL;
  Nx=gN(::Nx,g);
  Ny=gN(::Ny,g);
  //  spectrum=::spectrum;

  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");
  if(Nx != Ny) msg(ERROR,"Nx and Ny must be equal");
  if(radix != 1 && radix != 4) 
    msg(ERROR,"only radix 1 (trivial) or 4 decimations enabled");

  k0=gk(1.0,myg);
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  xorigin=mx-1;

  if(g == 0) {
    Invisible=Invisible2=0;
  } else {
    Invisible=gm(::Nx,g-1)*gk(1.0,myg-1)/gk(1.0,myg);
    Invisible2=Invisible*Invisible;
  }

  NY[OMEGA]=Nx*my;
  NY[TRANSFER]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;

  cout << "\nALLOCATING FFT BUFFERS" << endl;
  align=sizeof(Complex);

  Allocator(align);  // FIXME: refers to kernel.h, but should go somewhere else?

  //Dimension(T,nshells);

  unsigned Nxmy=Nx*my;
  unsigned nbuf=3*Nxmy;
  //  unsigned Nx0=Nx+xpad;//unused so far...
  //  unsigned Ny0=Ny+ypad; //unused so far...
  //  int my0=Ny0/2+1; //unused so far...
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

  InitialCondition=MDNS_Vocabulary.NewInitialCondition(ic);
  w.Set(Y[OMEGA]);
  InitialCondition->Set(w,NY[OMEGA]);
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);

  for(unsigned i=0; i < nshells; i++)
    Y[EK][i]=0.0;
}

void Grid::NonLinearSource(const vector2& Src, const vector2& Y, double t)
{
  DNS::NonLinearSource(Src,Y,t);
}

void Grid::ComputeInvariants(Real& E, Real& Z, Real& P)
{
  FunctRRPtr F = new Invariants;
  loopwF(F,3,&Z,&E,&P);
}

/****** Vocabulary *****/
MDNSVocabulary::MDNSVocabulary()
{
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
  VOCAB(radix,1,INT_MAX,"Radix number for grid decimation");
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

void MDNS::InitialConditions()
{
  nfields=Nfields;

  //***** Vocabulary *****//
  Ngrids=::Ngrids;
  MultiProblem::InitialConditions(Ngrids);

  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  for(unsigned g=0; g < Ngrids; ++g) {
    NY[Nfields*g+OMEGA]=gN(Nx,g)*gm(my,g);
    //NY[Nfields*g+TRANSFER]=
    //NY[Nfields*g+MOMENT]=
  }


  Allocator(align); // allocates MultiIntegrator

  G.Allocate(Ngrids);
  for(unsigned g=0; g < Ngrids; ++g) 
    G[g]=new Grid();
  for(unsigned g=0; g < Ngrids; ++g) 
    G[g]->InitialConditions(g);

}

void MDNS::Project(unsigned g) 
{
  if(radix == 1) {
    // FIXME
  } 
  if(radix == 4) {
    // FIXME
  }
}

void MDNS::Prolong(unsigned g) 
{
  if(radix == 1) {
    // FIXME
  } 
  if(radix == 4) {
    // FIXME
  }
}

void MDNS::Output(int it)
{
  // FIXME
}

void MDNS::ComputeInvariants(Real& E, Real& Z, Real& P)
{
  E=Z=P=0.0;

  Real tempE=0.0, tempZ=0.0, tempP=0.0;
  for(unsigned g=0; g < Ngrids; ++g) {
    // add up the individual invariants
    G[g]->ComputeInvariants(tempE,tempZ,tempP);
    Real scale=pow(radix,g);
    scale=1; // FIXME: temporary
    E += scale*tempE;
    Z += scale*tempZ;
    P += scale*tempP;
  }
}

void MDNS::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

void DNS::Output(int it) // Necessary for a table in triad somewhere
{
}

void MDNS::Source(const vector2& Src, const vector2& Y, double t)
{
  ConservativeSource(Src,Y,t);
  NonConservativeSource(Src,Y,t);
}

void MDNS::ConservativeSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  //Transfer(Src,Y);
  LinearSource(Src,Y,t); // FIXME: is this really a conservative source?
}

void MDNS::NonConservativeSource(const vector2& Src, const vector2& Y, double t)
{
  //Moments(Src,Y,t);
}

void MDNS::ExponentialSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  //Transfer(Src,Y);
  NonConservativeSource(Src,Y,t);
}

void MDNS::LinearSource(const vector2& Src, const vector2& Y, double t) 
{
  G[grid]->LinearSource(Src,Y,t);
}

void MDNS::NonLinearSource(const vector2& Src, const vector2& Y, double t) 
{
  G[grid]->NonLinearSource(Src,Y,t);
}

void MDNS::Stochastic(const vector2& Y, double t, double dt) 
{
  //G[grid]->Stochastic(T,y,t,dt);
}

/*
void MDNS::Moments(const vector2& Src, const vector2& Y, double t) 
{
  G[grid]->Moments(Src,Y,t);
}
  
void MDNS::Transfer(const vector2& Src, const vector2& Y) 
{
  G[grid]->Transfer(Src,Y,t);
}
*/
