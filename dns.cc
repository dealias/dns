#include "options.h"
#include "dns.h"

const double ProblemVersion=1.0;

#ifndef DEPEND
#if !COMPLEX
#error Navier requires COMPLEX=1 in options.h
#endif
#endif

const char *problem="Direct Numerical Simulation of Turbulence";

const char *method="DNS";
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

// other global variables:
int xpad=1;
int ypad=1;
DNS *DNSProblem;

InitialConditionBase *InitialCondition;
ForcingBase *Forcing;

class DNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return problem;}
  const char *Abbrev() {return "DNS";}
  DNSVocabulary();

  Table<InitialConditionBase> *InitialConditionTable;
  InitialConditionBase *NewInitialCondition(const char *& key) {
    return InitialConditionTable->Locate(key);
  }

  Table<ForcingBase> *ForcingTable;
  ForcingBase *NewForcing(const char *& key) {
    return ForcingTable->Locate(key);
  }
};

DNSVocabulary DNS_Vocabulary;

DNSVocabulary::DNSVocabulary()
{
  Vocabulary=this;

  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  VOCAB(spectrum,0,1,"Spectrum flag (0=off, 1=on)");
  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");

  METHOD(DNS);

  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  VOCAB(icalpha,0.0,0.0,"initial condition parameter");
  VOCAB(icbeta,0.0,0.0,"initial condition parameter");
  INITIALCONDITION(Zero);
  INITIALCONDITION(Constant);
  INITIALCONDITION(Equipartition);

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
  FORCING(None);
  FORCING(WhiteNoiseBanded);
}


// DNS setup routines
DNS::DNS()
{
  DNSProblem=this;
  check_compatibility(DEBUG);
  ConservativeIntegrators(DNS_Vocabulary.IntegratorTable,this);
  ExponentialIntegrators(DNS_Vocabulary.IntegratorTable,this);
}

DNS::~DNS()
{
  fftwpp::deleteAlign(block);
}


void DNS::InitialConditions()
{
  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");

  k0=1.0;
  k02=k0*k0;

  mx=(Nx+1)/2;
  my=(Ny+1)/2;

  xorigin=mx-1;
  origin=xorigin*my;
  nshells=spectrum ? (unsigned) (hypot(mx-1,my-1)+0.5) : 0;


  NY[OMEGA]=Nx*my;
  NY[TRANSFER]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;

  cout << "\nALLOCATING FFT BUFFERS" << endl;
  size_t align=sizeof(Complex);

  Allocator(align);

  Dimension(T,nshells);

  w.Dimension(Nx,my);
  f0.Dimension(Nx,my);

  unsigned int Nxmy=Nx*my;
  unsigned int nbuf=3*Nxmy;
  unsigned int Nx0=Nx+xpad;
  unsigned int Ny0=Ny+ypad;
  int my0=Ny0/2+1;
  if(movie)
    nbuf=max(nbuf,Nx0*my0);

  block=fftwpp::ComplexAlign(nbuf);
  f1.Dimension(Nx,my,block);
  g0.Dimension(Nx,my,block+Nxmy);
  g1.Dimension(Nx,my,block+2*Nxmy);

  F[1]=f1;
  G[0]=g0;
  G[1]=g1;

  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,2);

  Allocate(count,nshells);

  if(movie) {
    buffer.Dimension(Nx0,my0,block);
    wr.Dimension(Nx0,2*my0,(Real *) block);
    Padded=new fftwpp::ExplicitHConvolution2(Nx0,Ny0,mx,my,block);
  }

  InitialCondition=DNS_Vocabulary.NewInitialCondition(ic);
  w.Set(Y[OMEGA]);
  InitialCondition->Set(w,NY[OMEGA]);
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);

  for(unsigned i=0; i < nshells; i++)
    Y[EK][i]=0.0;

  Forcing=DNS_Vocabulary.NewForcing(forcing);

  if(dynamic && false) {
    Allocate(errmask,ny);
    for(unsigned i=0; i < ny; ++i)
      errmask[i]=1;

    array2<int> omegamask(Nx,my,(int *) errmask+(Y[OMEGA]-y));
    for(unsigned i=0; i <= xorigin; i++)
      omegamask(i,0)=0;
  }

  tcount=0;
  if(restart) {
    Real t0;
    ftin.open(Vocabulary->FileName(dirsep,"t"));
    while(ftin >> t0, ftin.good()) tcount++;
    ftin.close();
  }

  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");

  if(!restart) {
    remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
    remove_dir(Vocabulary->FileName(dirsep,"transfer"));
  }

  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
  mkdir(Vocabulary->FileName(dirsep,"transfer"),0xFFFF);

  errno=0;

  if(output)
    open_output(fwk,dirsep,"wk");

  if(movie)
    open_output(fw,dirsep,"w");
}

