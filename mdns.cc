#include "dnsbase.h"
#include "MultiIntegrator.h"
#include "Conservative.h"

using namespace fftwpp;

const double ProblemVersion=1.0;

#ifndef DEPEND
#if !COMPLEX
#error Navier requires COMPLEX=1 in options.h
#endif
#endif

const char *problem="Multispectral Direct Numerical Simulation of Turbulence";

const char *method="MDNS";
const char *integrator="MultiIntegrator";
const char *ic="Constant";
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
enum PRTYPE {AREA,POINT};
unsigned prtype=AREA;


// global variables

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
unsigned Ngrids=2;
const char *subintegrator="rk4";

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
  virtual ~FunctRR() {};
  virtual void f(Real w2,Real k2, int n,va_list args) = 0;
};
typedef FunctRR* FunctRRPtr;

unsigned gN(unsigned N, unsigned g) {// FIXME: should this be part of MDNS?
  return radix == 1 ? pow(2,(int) g)*(N+1)-1 : N;
};

unsigned gm(unsigned m, unsigned g) { // FIXME: should this be part of MDNS?
  return radix == 1 ? pow(2,(int) g)*m : m;
};


//***** MDNS derived from MultiProblem *****//
class MDNS : public MultiProblem {
private:
  unsigned mx, my; // size of data arrays

  unsigned glast;
  void ComputeInvariants(const vector2 & Y, Real& E, Real& Z, Real& P);
  array2<Complex> wa, wb; // Vorticity field for projection/prolongation
  unsigned nshells;
  vector spectra;
  ofstream Mfevt;
  oxstream fekvk;
  ofstream ft;
  oxstream fprolog;
  Real k0;
  int tcount;
public:
  MDNS();
  ~MDNS();

  Real gk(unsigned g) {return k0*pow(sqrt((double) radix),g);};

  //***** Grid class based on DNSBase *****//
  class Grid : public DNSBase {
  private:
    unsigned Invisible;
    int Invisible2;
    Complex *wSblock;
    Complex *SrcwSblock;
    array2<Complex> wS, SrcwS;
    unsigned sNx, smx, smy, sxorigin, xdiff;
    void NLDimension();
    MDNS * parent;
    unsigned nfields;
    bool lastgrid;
    
    unsigned mym1;
    
    vector Sp; // Spectrum
    unsigned shellsbelow, myshells;


    DNSBase * smallDNSBase; // for subgrid non-linearity calculation

  public:
    Grid();
    Grid(unsigned, MDNS *, bool);
    Grid(MDNS *prob, unsigned g, const vector2& Y0);
    ~Grid();
    void AttachTo(MDNS *prob, const vector2 & Y);
    void SetParams();
    unsigned myg;

    array2<Real> tildeB;

    array1<unsigned>::opt count;
    void setcount();
    void setcountoverlap(array1<unsigned>::opt &);
    //unsigned extrashells;

    void settolast() {lastgrid=1;};
    unsigned getnshells() {return nshells;};
    unsigned getNx() {return Nx;};
    unsigned getmy() {return my;};
    unsigned getInvisible() {return Invisible;};
    void setshellsbelow(const unsigned i) {shellsbelow=i;};
    unsigned getshellsbelow() {return shellsbelow;};
    void setmyshells(const unsigned i) {myshells=i;};
    unsigned getmyshells() {return myshells;};
    void settow(array2<Complex> & w0) {
      Set(w0,w);
      Dimension(w0,w);
    };
    //unsigned nshellbelow;

    void InitialConditions(unsigned g);
    //Real gk(Real k, unsigned g) {return k*pow(sqrt((double) radix),g);};
  
    void NonLinearSource(const vector & Src, const vector & Y, double t);
    void NonConservativeSource(const vector & Src, const vector & Y, double t);
    void LinearSource(const vector & Src, const vector & Y, double t);
    void Transfer(const vector2 & Src, const vector2 & Y);
    void Spectrum(vector& S, const vector& w0);    
    void SpectrumOverlap(vector& S);

    void ComputeInvariants(const vector2 & Y, Real& E, Real& Z, Real& P);
    
    void coutw() {cout << "mygrid is "<<myg<<endl<< w << endl;}
  
    void loopwF(const FunctRRPtr,int n,...);
    /*
    class Invariants : public FunctRR {
      void f(const Real w2, const Real k02, int n,va_list args) {
	double * Z=va_arg(args,double *);
	double * E=va_arg(args,double *);
	double * P=va_arg(args,double *);
	*Z += w2;
	*E += w2/k2;
	*P += w2*k2;
      }
    };
    */
  };
  array1<Grid *> G;
  

  enum Field {OMEGA,TRANSFER,EK,Nfields};
  unsigned getnfields(unsigned g) {
    return Nfields;
    //return (g == glast && spectrum) ? Nfields : 1;
  };
  unsigned getNfields(unsigned g) {return Nfields;};
  unsigned getnshells(unsigned g) {
    return (unsigned) (hypot(gm(mx,g)-1,gm(my,g)-1)+0.5);
  };

  //array1<unsigned>::opt count;
  //void setcount();

  Table<InitialConditionBase> *InitialConditionTable;
  void InitialConditions();

  // Source functions
  void Source(const vector2&, const vector2&, double);
  void ConservativeSource(const vector2& Src, const vector2& Y, double t);
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void Transfer(const vector2& Src, const vector2& Y);
  void ExponentialSource(const vector2& Src, const vector2& Y, double t);
  void Stochastic(const vector2& Y, double t, double dt);
  void Spectrum(vector& S, const vector& w0);
  
  void Project(unsigned ga);
  void Prolong(unsigned gb);

  void IndexLimits(unsigned int& start, unsigned int& stop,
		   unsigned int& startT, unsigned int& stopT,
		   unsigned int& startM, unsigned int& stopM) {
    // FIXME: this might not be set right for exponential integrators
    start=Start(OMEGA);
    stop=Stop(OMEGA);
    startT=Start(TRANSFER);
    stopT=Stop(TRANSFER);
    startM=Start(EK);
    stopM=Stop(EK);
  }

  // Output functions
  Real getSpectrum(unsigned i) {
    if(radix == 4) {
      for(unsigned g=0; g < Ngrids; ++g) {
	const unsigned nbelow=G[g]->getshellsbelow();
	const unsigned last=nbelow+G[g]->getmyshells();

	if (i <=  nbelow ||  i > last) {
	  // this shell is not part of this grid
	} else { // FIXME
	  // two cases:
	  
	  // 1. the shell is only on grid g.
	  
	  
	  // 2. the shell is between two grids.
	}
      }
      return 0.0;
    }

  //cerr << "radix=1 spectrum still needs work" << endl; exit(1);
  return 0.0;

  };

  void Computek(DynVector<unsigned>&);
  Real getk(unsigned i) {return (Real) i;};

  void FinalOutput();
  void Initialize();
  void Output(int it);
};

//***** Global problem *****//
MDNS *MDNSProblem;

// ***** Grid functions *****//
void MDNS::Grid::loopwF(const FunctRRPtr F,int n,...)
{
  va_list args;
  va_start(args,n);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j >= Invisible || I2 >= Invisible2) {
	Real w2=abs2(wi[j]);
	Real k2=k02*(I2+j*j);
	F->f(w2,k2,n,args);
      }
    }
  }
  va_end(args);
};


MDNS::Grid::Grid()
{
  lastgrid=0;
}

MDNS::Grid::Grid(unsigned g, MDNS * parent0, bool islast=0) 
{
  parent=parent0;
  myg=g;
  lastgrid=islast;
  SetParams();
}

MDNS::Grid::Grid(MDNS *prob, unsigned g, const vector2 & Y)
{
  parent=prob;
  w.Set(Y[OMEGA]);
}

MDNS::Grid::~Grid()
{
  fftwpp::deleteAlign(block);
  if(myg > 0) {
    fftwpp::deleteAlign(wSblock);
    fftwpp::deleteAlign(SrcwSblock);
  }
}

void MDNS::Grid::AttachTo(MDNS *prob, const vector2 & Y)
{
  parent=prob;
  nfields=parent->getNfields(myg);
  w.Set(Y[nfields*myg+OMEGA]);
  Set(Sp,Y[nfields*myg+EK]);
}

void MDNS::Grid::SetParams()
{
  Nx=gN(::Nx,myg);
  Ny=gN(::Ny,myg);
  k0=parent->gk(myg);
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  mym1=my-1;
  xorigin=mx-1;

  w.Dimension(Nx,my);
  origin=xorigin*my;
  
  if(myg == 0) {
    Invisible=Invisible2=0;
  } else {
    if(radix==1) 
      Invisible=my;

    if(radix==2) {
      cerr << "radix two case not implimented" << endl;
      exit(1);
    }

    if(radix==4)
      Invisible=gm(my,myg-1)/2;

    Invisible2=Invisible*Invisible;
  }
}

void MDNS::Grid::InitialConditions(unsigned g)
{
  myg=g;
  //cout << Ngrids << endl; exit(1);

  // load vocabulary from global variables
  nuH=::nuH;
  nuL=::nuL;

  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");
  if(Nx != Ny) msg(ERROR,"Nx and Ny must be equal");
  if(radix != 1 && radix != 4) 
    msg(ERROR,"only radix 1 (trivial) or 4 decimations enabled");

  nshells=0;

  if(spectrum) {
    nshells=lastgrid ? (unsigned) (hypot(mx-1,my-1)+0.5) : my;
    Dimension(Sp,nshells); // FIXME: allocated?
    for(unsigned i=0; i < nshells; i++)  Sp[i]=Complex(0.0,0.0);
  }

  //  unsigned Nx0=Nx+xpad;//unused so far...
  //  unsigned Ny0=Ny+ypad; //unused so far...
  //  int my0=Ny0/2+1; //unused so far...

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;

  cout << "\nALLOCATING FFT BUFFERS" << endl;

  block=ComplexAlign(3*Nx*my);
  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,2);
  NLDimension();

  if(myg > 0) { 
    // buffer for calculating small-small-small calculations
    sNx=Invisible;
    smx=smy=(Invisible+1)/2;
    sxorigin=(sNx+1)/2-1;
    unsigned nS=sNx*smy;
    xdiff=xorigin-sxorigin;

    wSblock=fftwpp::ComplexAlign(nS);
    SrcwSblock=fftwpp::ComplexAlign(nS);
    wS.Dimension(sNx,smy,wSblock);
    SrcwS.Dimension(sNx,smy,SrcwSblock);
    smallDNSBase=new DNSBase(sNx,smy,k0); // TODO: share memory block

    tildeB.Dimension(2*Invisible+1,Invisible+1);
    Allocate(tildeB,(2*Invisible+1)*(Invisible+1));
    for(unsigned i=0; i < 2*Invisible+1; ++i)
      tildeB[i][0]=0.0;
  }

  InitialCondition=MDNS_Vocabulary.NewInitialCondition(ic);
  InitialCondition->Set(w,Nx*my);


  Allocate(count,nshells);
  for(unsigned k=0; k < nshells; ++k) 
    count[k]=0;
    
  // FIXME: initialize spectrum to zero as well?

  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);
}

void MDNS::Grid::NLDimension()
{
  w.Dimension(Nx,my);
  f0.Dimension(Nx,my);
  f1.Dimension(Nx,my,block);
  g0.Dimension(Nx,my,block+Nx*my);
  g1.Dimension(Nx,my,block+2*Nx*my);
  F[0]=f0;
  F[1]=f1;
  G[0]=g0;
  G[1]=g1;
}

void MDNS::Grid::LinearSource(const vector& source, const vector& w0, double t)
{
  // FIXME: overcounts
  DNSBase::LinearSource(source,w0,t);
}

void MDNS::Grid::NonLinearSource(const vector& source, const vector& w0,
				 double t)
{
  w.Set(w0);
  f0.Set(source);

  //  cout << "w from Grid::NonLinearSource, myg=" << myg << endl;
  //  cout << w << endl;

  DNSBase::NonLinearSource(source,w,t); // FIXME: w0 here? uninitialied?

  if(myg > 0) {
    //copy overlapping modes to wS
    for(unsigned i=0; i < sNx; ++i) {
      vector wi=w[i+xdiff];
      vector wSi=wS[i];
      for(unsigned j=0; j < smy; ++j)
	wSi[j]=wi[j];
    }
    fftwpp::HermitianSymmetrizeX(smx,smy,sxorigin,wS);

    // find the nonlinear interaction for small-small-small
    smallDNSBase->NonLinearSource(SrcwS,wS,t); 
 
    // subtract small-small-small source
    for(unsigned i=0; i < sNx; ++i) {
      vector f0i=f0[i+xdiff];
      vector SrcwSi=SrcwS[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < smy; ++j)
	f0i[j] -= SrcwSi[j];
    }
    fftwpp::HermitianSymmetrizeX(mx,my,xorigin,f0);
  }
}

void MDNS::Grid::NonConservativeSource(const vector& source, const vector& w0, 
				       double)
{
  if(spectrum) {
    vector S;
    Set(S,source);
    Dimension(S,nshells);
    Spectrum(S,w0);
  }
}

void MDNS::Grid::Transfer(const vector2 & Src, const vector2 & Y)
{
  // FIXME
}

void MDNS::Grid::setcount()
{
  if(spectrum) {
    Real kbound=lastgrid ? hypot(mx,my)+1: mym1;
    for(unsigned i=0; i < Nx; i++) {
      const int I=(int) i-(int) xorigin;
      const int I2=I*I;
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	if(j >= Invisible || I2 >= Invisible2) {
	  const Real kint=sqrt(I2+j*j);
	  if(kint <= kbound) {
	    count[(unsigned)(kint-0.5)] += 1;
	  }
	}
      }
    }
  }
}

void MDNS::Grid::setcountoverlap(array1<unsigned>::opt &Count)
{
  if(spectrum) {
    Real overlambda=1.0;
    if(radix==4) overlambda=0.5;
    Real kbound=mym1;
    for(unsigned i=0; i < Nx; i++) {
      const int I=(int) i-(int) xorigin;
      const int I2=I*I;
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	if(j >= Invisible || I2 >= Invisible2) {
	  const Real kint=sqrt(I2+j*j);
	  if(kint-0.5 > kbound) {
	    Count[(unsigned)(overlambda*(kint-0.5))] += 1; // FIXME: radix?
	  }
	}
      }
    }
  }
}

void MDNS::Spectrum(vector& S, const vector& w0)
{
  G[grid]->Spectrum(S,w0);

  if(grid != 0) {
    // FIXME: if grid!=0, then this should take the previous grid's w0.
    G[grid-1]->SpectrumOverlap(S);
  }
}

void MDNS::Grid::Spectrum(vector& S, const vector& w0)
{
  w.Set(w0);
  S.Load(Complex(0.0,0.0));

  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j > Invisible || I2 > Invisible2) {
	const Real kint=sqrt(I2+j*j);
	const Real k=k0*kint;
	const Real overk=1.0/k;
	const Real w2=abs2(wi[j]);
	if(kint <= my) {
	  S[(unsigned)(kint-0.5)].re += w2*overk; 
	  //S[(unsigned)(k-0.5)].im += nuk(k2)*w2; // FIXME: nuk set?
	}
      }
    }
  }
}

void MDNS::Grid::SpectrumOverlap(vector& S)
{
  settow(w);
  S.Load(Complex(0.0,0.0));

  Real overlambda=1.0;
  if(radix==4) overlambda=0.5;
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j > Invisible || I2 > Invisible2) {
	const Real kint=sqrt(I2+j*j);
	const Real k=k0*kint;
	const Real overk=1.0/k;
	const Real w2=abs2(wi[j]);
	if(kint > my) {
	  S[(unsigned)(overlambda*kint-0.5)].re += w2*overk; 
	  //S[(unsigned)(overlambda*kint-0.5)].im += nuk(k2)*w2; 
	}
      }
    }
  }
}


void MDNS::Grid::ComputeInvariants(const vector2 & Y, Real& E, Real& Z, Real& P)
{
  w.Set(Y[OMEGA]);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j > Invisible || I2 > Invisible2) {
	Real w2=abs2(wi[j]);
	Real k2=k02*(I2+j*j);
	Z += w2;
	E += w2/k2;
	P += w2*k2;
      }
    }
  }
  //  FunctRRPtr F = new Invariants;
  //  loopwF(F,3,&Z,&E,&P);
}

/****** Vocabulary *****/
MDNSVocabulary::MDNSVocabulary()
{
  Vocabulary=this;

  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  //  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  VOCAB(spectrum,0,1,"Spectrum flag (0=off, 1=on)");
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
  VOCAB(prtype,0,1,"Synchronization scheme (0=area, 1=coincident)");
}

MDNS::MDNS() 
{ 
  nfields = Nfields;
  MDNSProblem=this;
  MProblem=this;
  check_compatibility(DEBUG);
  align=sizeof(Complex);
  grid=0;
  ConservativeIntegrators(MDNS_Vocabulary.IntegratorTable,this);
}

MDNS::~MDNS()
{
}

//***** wrappers for output curves *****//
Real curve_Spectrum(unsigned i) {return MDNSProblem->getSpectrum(i);}
Real curve_nuk(unsigned i) {return 0.0;} // FIXME
Real curve_k(unsigned i) {return MDNSProblem->getk(i);}

void MDNS::InitialConditions()
{
  //***** Vocabulary *****//
  Ngrids=::Ngrids;
  glast=Ngrids-1;
  saveF=OMEGA;
  MultiProblem::InitialConditions(Ngrids);
  MProblem=this;
  k0=1.0; // FIXME: this should eventually come from command-line.

  mx=(Nx+1)/2;
  my=(Ny+1)/2;

  nshells=0;

  //  Allocate(spectra,nshells);
  //  for (unsigned i=0; i < nshells; ++i) spectra[i]=0.0;

  for(unsigned g=0; g < Ngrids; ++g) {
    //nfields=getNfields(g);
    NY[nfields*g+OMEGA]=gN(Nx,g)*gm(my,g);
    NY[nfields*g+TRANSFER]=0; // getnshells(g); FIXME 
    NY[nfields*g+EK]=getnshells(g);
  }

  Allocator(align); // allocate MultiIntegrator
  
  G.Allocate(Ngrids);

  for(unsigned g=0; g < Ngrids; ++g) {
    G[g]=new Grid(g,this,g==glast);
    G[g]->AttachTo(this,Y);
    G[g]->InitialConditions(g);
  }
  for (unsigned g=1; g< Ngrids; g++) 
    Project(g);

  if(spectrum) {
    for(unsigned g=0; g < Ngrids; ++g) {
      G[g]->setcount();
      if(g!=glast) 
	G[g]->setcountoverlap(G[g+1]->count);
    }
  }
    
  if(spectrum) {
    unsigned lambda;
    if(radix==1) lambda=1;
    if(radix==4) lambda=2;
    unsigned nextra=gm(my,0);
    nshells += nextra;
    //cout << nextra << endl;
    G[0]->setshellsbelow(0);
    for(unsigned g=1; g < glast; ++g) {
      nextra=gm(my,g) - gm(my,g-1)/lambda;
      G[g]->setshellsbelow(nshells);
      nshells += nextra;
      //cout << nextra << endl;
    }
    unsigned gmx=gm(mx,glast);
    unsigned gmy=gm(my,glast);
    nextra=(unsigned) (hypot(gmx-1,gmy-1) + 0.5*pow(lambda,glast) 
		       - gm(my,glast-1)/lambda);
    G[glast]->setshellsbelow(nshells);
    nshells += nextra;
    //cout << nextra << endl;
  }

  //***** output *****//
  open_output(fprolog,dirsep,"prolog",0);
  out_curve(fprolog,curve_k,"kc",nshells);
  out_curve(fprolog,curve_k,"kb",nshells);
  fprolog.close();

  tcount=0;
  //setcount();
  
  open_output(Mfevt,dirsep,"evt");
  open_output(ft,dirsep,"t");
  if(!restart) remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
}

void MDNS::Project(unsigned gb) 
{
  //return;
  if(verbose > 2) cout << "project onto " << G[gb]->myg << endl;
  unsigned ga=gb-1;

  const unsigned bInvisible=G[gb]->getInvisible();
  const unsigned axorigin=G[ga]->getxorigin();
  const unsigned bxorigin=G[gb]->getxorigin();
  const unsigned aNx=G[ga]->getNx();
  const unsigned amy=G[ga]->getmy();
  const unsigned bNx=G[gb]->getNx();
  const unsigned bmy=G[gb]->getmy();
  const unsigned dx=bxorigin-axorigin;

  wa.Dimension(aNx,amy);
  wb.Dimension(bNx,bmy);

  G[ga]->settow(wa);
  G[gb]->settow(wb);

  if(radix == 1) {
    for(unsigned int i=0; i < aNx; i++) {
      vector wai=wa[i];
      vector wbi=wb[i+dx];
      for(unsigned j=0; j < amy; ++j) {
	wbi[j]=wai[j];
      }
    }
  }

  //cout << "wa:\n"<< wa<< endl << "wb:\n"<< wb << endl;

  if(radix == 4) {

    if(prtype==AREA) {
      const unsigned xstart=bInvisible;
      const unsigned xstop=bxorigin+bInvisible;
      
      unsigned I=1;
      for(unsigned i=xstart; i < xstop; ++i) {
	//int I= 2*((int)i - (int)bxorigin)+axorigin;
	
	//cout << "i="<<i<<" I="<<I<<endl;
	
	vector wai, waim, waip;
	Set(wai,wa[I]);
	Set(waim,wa[I-1]);
	Set(waip,wa[I+1]);
	Dimension(wai,aNx);
	Dimension(waim,aNx);
	Dimension(waip,aNx);
	
	vector wbi;
	Set(wbi,wb[i]);
	Dimension(wbi,wb[i]);
	
      for(unsigned j= i < axorigin ? 1 : 0 ; j < bInvisible; ++j) {
	//cout << "j="<<j<<endl;
	Real B2=abs2(wbi[j]);

 	const int aJ=2*j;
	const int aJp=aJ+1;
	const int aJm=aJ==0? aJp : aJ-1;
	
	// co-incident point 
	//wai[aJ]=Complex(I,aJ); // 00
	Real A2=abs2(wai[aJ]);
	
	// same row
	//waim[aJ]=Complex(I-1,aJ); // -0
	//waip[aJ]=Complex(I+1,aJ); // +0
	A2 += 0.5*(abs2(waim[aJ]) + abs2(waip[aJ]));

	// same column
	if(aJ) {  // 0-
	  //wai[aJm]=Complex(I,aJm);
	  A2 += 0.5*abs2(wai[aJm]);
	} else {
	  //wa[2*axorigin-I][aJp]=Complex(-I,aJm);
	  A2 += 0.5*abs2(wa[2*axorigin-I][aJp]);
	}
	//wai[aJp]=Complex(I,aJp); // 0+
	A2 += 0.5*abs2(wai[aJp]);

	// quad-projects:
	//waim[aJp]=Complex(I-1,aJp); // -+
	//waip[aJp]=Complex(I+1,aJp); // ++
	A2 += 0.25*(abs2(waim[aJp]) + abs2(waip[aJp]));
	if(aJ) { // -- and +-
	  //waim[aJm]=Complex(I-1,aJm);
	  //waip[aJm]=Complex(I+1,aJm);
	  A2 += 0.25*(abs2(waim[aJm]) + abs2(waip[aJm]));
	} else {
	  //wa[2*axorigin-I+1][aJm]=Complex(-I+1,aJm);
	  //wa[2*axorigin-I-1][aJm]=Complex(-I-1,aJm);
	  A2 += 0.25*(abs2(wa[2*axorigin-I+1][aJm]) +
		      abs2(wa[2*axorigin-I-1][aJm]));
	}

	A2 *= 0.25; // radix-4 correction
	if(B2) {
	  wbi[j] *= sqrt(A2/B2);
	  //cout << "sqrt(A2/B2)=" << sqrt(A2/B2) << endl;
	} else {
	  cerr << "energy for mode (" << i << ","<<j << ") on grid "<< gb
	       << " is zero in project."<<endl;
	  exit(1); // FIXME: work out something better for this case.
	}
      }
      I += 2;
    }

    
    { // copy to tildeB
      // TODO: make this a function?
      G[gb]->tildeB.Load(0.0); // necessary?
      const unsigned istop=2*bInvisible+1;
      const unsigned jstop=bInvisible+1;
      unsigned wi=(int) bxorigin-bInvisible;
      for(unsigned i=0; i < istop; ++i) {
	array1<Real> tildeBi;
	tildeBi.Set(G[gb]->tildeB[i]);
	tildeBi.Dimension(jstop);
	vector wbi;
	Set(wbi,wb[wi]);
	Dimension(wbi,wb[wi]);
	++wi;
	for(unsigned j=0; j < jstop; ++j)
	  tildeBi[j]=abs2(wbi[j]);
	// TODO: some of these are redundant from Hermitian symmerty.
      }
    }
    if(verbose > 2) {
      cout <<"wa:\n"<<wa<<endl << "wb:\n"<< wb;
      cout <<"\n~B:\n" << G[gb]->tildeB << endl;
    }
    /*
    for(unsigned i=0; i < aNx; i += 2) {
      int I= (int) i - (int) bxorigin;
      int bI=I/2;
      unsigned bi=bI + bxorigin;
      for(unsigned j= i < axorigin? 2 : 0; j <  amy; j +=2) {
	unsigned bj=j/2;
	Real bZ=abs2(wb[bi][bj]);
      }
    }
    */
    }
    
    if(prtype==POINT) { // FIXME: check
    unsigned xstart=bInvisible;
    unsigned xstop=bxorigin+bInvisible;
    unsigned I=1;
    for(unsigned i=xstart; i < xstop; ++i) {
      vector wai;
      Set(wai,wa[I]);
      Dimension(wai,aNx);
      vector wbi;
      Set(wbi,wb[i]);
      Dimension(wbi,wb[i]);
      for(unsigned j= i < axorigin ? 1 : 0 ; j < bInvisible; ++j)
	wbi[j]=wai[j+j];
      I += 2;
    }
  }
  }

  //cout << "wa:\n"<< wa<< endl << "wb:\n"<< wb << endl;
  //  exit(1);
  //  HermitianSymmetrizeX(bNx,bmy,bxorigin,wb);  // FIXME: messed up?
}
  
void MDNS::Prolong(unsigned ga)
{
  //return;
  if(verbose > 2) cout << "prolong onto " << G[ga]->myg << endl;
  unsigned gb=ga+1;

  //const unsigned aInvisible=G[ga]->getInvisible();
  const unsigned bInvisible=G[gb]->getInvisible();
  const unsigned axorigin=G[ga]->getxorigin();
  const unsigned bxorigin=G[gb]->getxorigin();
  const unsigned aNx=G[ga]->getNx();
  const unsigned amy=G[ga]->getmy();
  //const unsigned amx=G[ga]->getmx();
  const unsigned bNx=G[gb]->getNx();
  const unsigned bmy=G[gb]->getmy();
  const unsigned dx=bxorigin-axorigin;

  const unsigned dtilX=bInvisible+bInvisible;

  wa.Dimension(aNx,amy);
  wb.Dimension(bNx,bmy);

  G[ga]->settow(wa);
  G[gb]->settow(wb);

  if(radix == 1) {
    for(unsigned int i=0; i < aNx; i++) {
      vector wai=wa[i];
      vector wbi=wb[i+dx];
      for(unsigned j=0; j < amy; ++j) {
  	wai[j]=wbi[j];
      }
    }
  }

  if(radix == 4) {
    if(prtype==AREA) {
      const unsigned xstart=bInvisible;
      const unsigned xstop=bxorigin+bInvisible;

      //cout << "xstart=" << xstart << "\nxstop=" << xstop << endl;
      for(unsigned i=xstart; i < xstop; ++ i) {
	const int I= 2*((int)i - (int)bxorigin)+axorigin;
	const unsigned tildei=i+1-xstart;

	// cout << "i="<<i << " I="<<I << " tildei=" << tildei << endl;

	vector wai;
	Set(wai,wa[I]);
	Dimension(wai,wa[I]);
	vector waip;
	Set(waip,wa[I+1]);
	Dimension(waip,wa[I+1]);
	vector waim;
	Set(waim,wa[I-1]);
	Dimension(waim,wa[I-1]);
      
	array1<Real> tildeBim;
	tildeBim.Set(G[gb]->tildeB[tildei-1]);
	tildeBim.Dimension(G[gb]->tildeB[tildei-1]);
	array1<Real> tildeBi;
	tildeBi.Set(G[gb]->tildeB[tildei]);
	tildeBi.Dimension(G[gb]->tildeB[tildei]);
	array1<Real> tildeBip;
	tildeBip.Set(G[gb]->tildeB[tildei+1]);
	tildeBip.Dimension(G[gb]->tildeB[tildei+1]);

	vector wbi;
	Set(wbi,wb[i]);
	Dimension(wbi,wb[i]);
	vector wbim;
	Set(wbim,wb[i-1]);
	Dimension(wbim,wb[i-1]);
	vector wbip;
	Set(wbip,wb[i+1]);
	Dimension(wbip,wb[i+1]);
      
	// possible optimization: only some of these need to be recalculate.
	// most can be copied from one stage to the next, put into
	// different positions.
	// FIXME: figure out what to deal with division-by-zero case.
	for(unsigned j= i < axorigin ? 1 : 0 ; j < bInvisible; ++j) {
	  //cout << "j=" << j << endl;
	  const unsigned aJ=j+j;
	  const unsigned aJp=aJ+1;
	  const bool jmOK=aJ > 0;
	  const unsigned aJm=jmOK ? aJ-1 : aJp;

	  const Real Bij=abs2(wbi[j]);
	  const Real tildeBij=tildeBi[j];
	
	  const Real Bimj=abs2(wbim[j]);
	  const Real tildeBimj=tildeBim[j];

	  const Real Bipj=abs2(wbip[j]);
	  const Real tildeBipj=tildeBip[j];

	  const Real Bimjp=abs2(wbim[j+1]);
	  const Real tildeBimjp=tildeBim[j+1];
	
	  const Real Bipjm=jmOK ? abs2(wbip[j-1]) :abs2(wb[2*bxorigin-i-1][1]);
	  const Real tildeBipjm=jmOK ? tildeBip[j-1] :
	    G[gb]->tildeB[dtilX-tildei-1][1];

	  const Real Bijm=jmOK ? abs2(wbi[j-1]) : abs2(wb[2*bxorigin-i][1]);
	  const Real tildeBijm=jmOK ? tildeBi[j-1] : 
	    G[gb]->tildeB[dtilX-tildei][1];

	  const Real Bimjm=jmOK ? abs2(wbim[j-1]) : abs2(wb[2*bxorigin-i+1][1]);
	  const Real tildeBimjm=jmOK ? tildeBim[j-1] : 
	    G[gb]->tildeB[dtilX-tildei+1][1];

	  const Real Bipjp=abs2(wbip[j+1]);
	  const Real tildeBipjp=tildeBip[j+1];

	  const Real Bijp=abs2(wbip[j+1]);
	  const Real tildeBijp=tildeBi[j+1];

	  // co-incident point
	  wai[aJ] *= sqrt(Bij/tildeBij);

	  // same row
	  waim[aJ] *= sqrt((Bimj+Bij)/(tildeBimj+tildeBij));
	  //	cout << (Bimj+Bij)/(tildeBimj+tildeBij) << endl;
	  waip[aJ] *= sqrt((Bij+Bipj)/(tildeBij+tildeBipj));

	  // same column
	  wai[aJm] *= sqrt((Bijm+Bij)/(tildeBijm+tildeBij));
	  wai[aJp] *= sqrt((Bij+Bijp)/(tildeBij+tildeBijp));
	
	  // quad-prolongs:
	  // --
	  waim[aJm] *= sqrt((Bij+Bimj+Bimjm+Bijm)
			    /(tildeBij+tildeBimj+tildeBimjm+tildeBijm));
	  //	cout << (Bij+Bimj+Bimjm+Bijm)
	  //	  /(tildeBij+tildeBimj+tildeBimjm+tildeBijm) << endl;

	  // -+
	  waim[aJ] *= sqrt((Bij+Bijp+Bimjp+Bimj)
			   /(tildeBij+tildeBijp+tildeBimjp+tildeBimj));
	  // +-
	  waip[aJm] *= sqrt((Bij+Bijm+Bipjm+Bipj)
			    /(tildeBij+tildeBijm+tildeBipjm+tildeBipj));
	  // ++
	  waip[aJm] *= sqrt((Bij+Bipj+Bipjp+Bijp)
			    /(tildeBij+tildeBipj+tildeBipjp+tildeBijp));
	
	}
      }
    }

    if(prtype==POINT) {
      const unsigned xstart=bInvisible; // FIXME: check
      const unsigned xstop=bxorigin+bInvisible;
      
      unsigned I=1;
      for(unsigned i=xstart; i < xstop; ++i) {
	//cout << "i="<<i<<" I="<<I<<endl;
	vector wai;
	Set(wai,wa[I]);
	Dimension(wai,aNx);
	vector wbi;
	Set(wbi,wb[i]);
	Dimension(wbi,wb[i]);
	
	for(unsigned j= i < axorigin ? 1 : 0 ; j < bInvisible; ++j) {
	  //cout << "j="<<j<<endl;
	  wai[j+j]=wbi[j];
	}
	I += 2;
      }
    }
  }
  
  //HermitianSymmetrizeX(amx,G[ga]->getmy(),axorigin,wa); // FIXME: messed up?
  // maybe only on the overlapping modes?
  if(verbose > 2) {
    cout <<"wa:\n"<<wa<<endl << "wb:\n"<< wb;
    cout <<"\n~B:\n" << G[gb]->tildeB << endl;
  }
  // FIXME: copy stuff from spectra onto just onto lastgrid's Src[EK]
  //  exit(1);
}

void MDNS::Initialize()
{
  Mfevt << "#   t\t\t E\t\t\t Z" << endl;
}

void MDNS::Output(int it)
{
  Real E,Z,P;

  ComputeInvariants(Y,E,Z,P);
  Mfevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;
    
  if(spectrum) {
    ostringstream buf;
    buf << "ekvk" << dirsep << "t" << tcount; 
    open_output(fekvk,dirsep,buf.str().c_str(),0);
    out_curve(fekvk,t,"t");
    out_curve(fekvk,curve_Spectrum,"Ek",nshells);
    out_curve(fekvk,curve_nuk,"nuk*Ek",nshells);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");
  }
  tcount++;
  ft << t << endl;
}

void MDNS::ComputeInvariants(const vector2& Y, Real& E, Real& Z, Real& P)
{
  E=Z=P=0.0;

  Real tempE=0.0, tempZ=0.0, tempP=0.0;
  for(unsigned g=0; g < Ngrids; ++g) {
    // add up the individual invariants
    G[g]->ComputeInvariants(Y,tempE,tempZ,tempP);
    Real scale=pow((double) radix,(double) g);
    E += scale*tempE;
    Z += scale*tempZ;
    P += scale*tempP;
  }
}

void MDNS::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(Y,E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

void MDNS::Source(const vector2& Src, const vector2& Y, double t)
{
  ConservativeSource(Src,Y,t);
  NonConservativeSource(Src,Y,t);
}

void MDNS::ConservativeSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  Transfer(Src,Y);
  LinearSource(Src,Y,t); // is this really a conservative source?
}

/*
void MDNS::setcount() 
{
  count.Allocate(nshells);
  for(unsigned g=0; g < Ngrids; ++g) {
    unsigned gInvisible=G[g]->getInvisible();
    int gInvisible2=gInvisible*gInvisible;
    unsigned gNx=G[g]->getNx();
    unsigned gxorigin=G[g]->getxorigin();
    for(unsigned i=0; i < gNx; i++) {
      int I=(int) i-(int)gxorigin;
      int I2=I*I;
      for(unsigned j=i <= gxorigin ? 1 : 0; j < my; ++j) {
	if(j > gInvisible || I2 > gInvisible2) {
	}
      }
    }
  }
  exit(1);
}
*/

void MDNS::NonConservativeSource(const vector2& Src, const vector2& Y, double t)
{
  vector w0,source;
  Set(w0,Y[OMEGA]);
  Set(source,Src[EK]);
  G[grid]->NonConservativeSource(source,w0,t);
  //Moments(Src,Y,t);
}

void MDNS::ExponentialSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  Transfer(Src,Y);
  NonConservativeSource(Src,Y,t);
}

void MDNS::LinearSource(const vector2& Src, const vector2& Y, double t) 
{
  vector w0,source;
  Set(w0,Y[OMEGA]);
  Set(source,Src[OMEGA]);
  G[grid]->LinearSource(source,w0,t);
}

void MDNS::NonLinearSource(const vector2& Src, const vector2& Y, double t) 
{
  vector w0,source;
  Set(w0,Y[OMEGA]);
  Set(source,Src[OMEGA]);
  G[grid]->NonLinearSource(source,w0,t);
}

void MDNS::Stochastic(const vector2& Y, double t, double dt) 
{
//  G[grid]->Stochastic(T,y,t,dt);
}

/*
void MDNS::Moments(const vector2& Src, const vector2& Y, double t) 
{
  G[grid]->Moments(Src,Y,t);
}
*/
  
void MDNS::Transfer(const vector2& Src, const vector2& Y) 
{
//  G[grid]->Transfer(Src,Y);
}

