#include "dnsbase.h"

const double ProblemVersion=1.0;

#if !COMPLEX
#error DNS requires COMPLEX=1 in options.h
#endif

using namespace utils;

const char *problem="Direct Numerical Simulation of Turbulence";

const char *method="DNS";
const char *integrator="RK5";
const char *ic="Equipartition";
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
const int Nforce=2;
Complex forces[Nforce];
int kxforces[Nforce];
int kyforces[Nforce];
Real deltaf=1.0;
unsigned movie=0;
unsigned rezero=0;
unsigned spectrum=1;
Real icalpha=1.0;
Real icbeta=1.0;
Real k0=1.0;
int randomIC=0;

class DNS : public DNSBase, public ProblemBase {
public:
  DNS();
  ~DNS();
  
  void InitialConditions();
  
  void Output(int);
  void FinalOutput();
  oxstream fprolog;

  void IndexLimits(unsigned& start, unsigned& stop,
		   unsigned& startT, unsigned& stopT,
		   unsigned& startM, unsigned& stopM) {
    start=Start(OMEGA);
    stop=Stop(OMEGA);
    startT=Start(TRANSFER);
    stopT=Stop(TRANSFER);
    startM=Start(EK);
    stopM=Stop(EK);
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    DNSBase::Source(Src,Y,t);
  }
  
  void Stochastic(const vector2&Y, double t, double dt) {
    DNSBase::Stochastic(Y,t,dt);
  }
  
  void Initialize();
};

DNS *DNSProblem;

InitialConditionBase *InitialCondition;
ForcingBase *Forcing;

class Zero : public InitialConditionBase {
public:
  const char *Name() {return "Zero";}
  
  Var Value(Real) {return 0.0;}
};

class Constant : public InitialConditionBase {
public:
  const char *Name() {return "Constant";}
  
  Var Value(Real) {return Complex(icalpha,icbeta);}
};

class Equipartition : public InitialConditionBase {
public:
  const char *Name() {return "Equipartition";}

  Var Value(Real k) {
// Distribute enstrophy evenly between the real and imaginary components.
    Real k2=k*k;
    Real v=icalpha+icbeta*k2;
    v=v ? k/sqrt(2.0*v) : 0.0;
    return randomIC ? v*expi(twopi*drand()) : sqrt(0.5)*Complex(v,v);
  }
};

// forcing
class None : public ForcingBase {
};

class ConstantBanded : public ForcingBase {
protected:  
  double h;
public:
  const char *Name() {return "Constant Banded";}
  
  void Init() {
    h=0.5*deltaf;
  }
  
  bool active(double k) {
    return abs(k-kforce) < h;
  }
  
  void Force(Complex& w, double& T, double k, double, double) {
    if(active(k)) {
      T += realproduct(force,w);
      w += force;
    }
  }
};

class ConstantList : public ForcingBase {
protected:  
public:
  const char *Name() {return "Constant List";}
  
  int active(double i, double j) {
    for(int index=0; index < Nforce; ++index)
      if(i == kxforces[index] && j == kyforces[index]) return index;
    return -1;
  }
  
  void Force(Complex& w, double& T, double k, double i, double j) {
    int index=active(i,j);
    if(index >= 0) {
      Complex force=forces[index];
      T += realproduct(force,w);
      w += force;
    }
  }
};

class WhiteNoiseBanded : public ConstantBanded {
  Complex f;
  Real etanorm;
public:
  const char *Name() {return "White-Noise Banded";}
  
  void Init(unsigned fcount) {
    etanorm=1.0/((Real) fcount);
  }
  
  bool Stochastic(double dt) {
    f=sqrt(2.0*dt*eta*etanorm)*crand_gauss();
    return true;
  }
  
   void ForceStochastic(Complex& w, double& T, double k) {
    if(active(k)) {
      T += realproduct(f,w)+0.5*abs2(f);
      w += f;
    }
  }
};

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
  VOCAB(randomIC,0,1,"randomize the initial conditions?");
  INITIALCONDITION(Zero);
  INITIALCONDITION(Constant);
  INITIALCONDITION(Equipartition);

  VOCAB(k0,0.0,0.0,"spectral spacing coefficient");
  VOCAB(nuH,0.0,REAL_MAX,"High-wavenumber viscosity");
  VOCAB(nuL,0.0,REAL_MAX,"Low-wavenumber viscosity");
  VOCAB(pH,0,0,"Power of Laplacian for high-wavenumber viscosity");
  VOCAB(pL,0,0,"Power of Laplacian for molecular viscosity");

  VOCAB_NOLIMIT(forcing,"Forcing type");
  ForcingTable=new Table<ForcingBase>("forcing");

  VOCAB(eta,0.0,REAL_MAX,"vorticity injection rate");
  VOCAB(force,(Complex) 0.0, (Complex) 0.0,"constant external force");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB_ARRAY(kxforces,"kx force wavenumbers");
  VOCAB_ARRAY(kyforces,"ky force wavenumbers");
  VOCAB_ARRAY(forces,"force ampligudes");
  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");
  FORCING(None);
  FORCING(ConstantBanded);
  FORCING(ConstantList);
  FORCING(WhiteNoiseBanded);
}

// DNSBase setup routines
DNS::DNS()
{
  DNSProblem=this;
  check_compatibility(DEBUG);
  ConservativeIntegrators(DNS_Vocabulary.IntegratorTable,this);
  ExponentialIntegrators(DNS_Vocabulary.IntegratorTable,this);
}

DNS::~DNS()
{
  deleteAlign(block);
}

// wrapper for outcurve routines
class cwrap{
public:
  static Real Spectrum(unsigned int i) {return DNSProblem->getSpectrum(i);}
  static Real Dissipation(unsigned int i) {return DNSProblem->Dissipation(i);}
  static Real Pi(unsigned int i) {return DNSProblem->Pi(i);}
  static Real Eta(unsigned int i) {return DNSProblem->Eta(i);}
  static Real kb(unsigned int i) {return DNSProblem->kb(i);}
  static Real kc(unsigned int i) {return DNSProblem->kc(i);}
};

void DNS::Initialize()
{
  DNSBase::Initialize();
}

void DNS::InitialConditions()
{
  fftw::maxthreads=threads;
  
  // load vocabulary from global variables
  Nx=::Nx;
  Ny=::Ny;
  nuH=::nuH;
  nuL=::nuL;

  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");

  k0=::k0;
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  xorigin=mx-1;
  origin=xorigin*my;
  nshells=spectrum ? (unsigned) (hypot(mx-1,my-1)+0.5) : 0;

  NY[PAD]=my;
  NY[OMEGA]=Nx*my; // Allow for X Hermitian padding
  NY[TRANSFER]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;
  cout << "\nALLOCATING FFT BUFFERS" << endl;
  size_t align=sizeof(Complex);

  Allocator(align);

  Dimension(E,nshells);
  Dimension(T,nshells);

  w.Dimension(Nx,my);
  S.Dimension(Nx,my);

  unsigned int Nx1my=(Nx+1)*my;
  unsigned int nbuf=3*Nx1my;
  
  block=ComplexAlign(nbuf);
  f0.Dimension(Nx+1,my,-1,0);
  f1.Dimension(Nx+1,my,block,-1,0);
  g0.Dimension(Nx+1,my,block+Nx1my,-1,0);
  g1.Dimension(Nx+1,my,block+2*Nx1my,-1,0);

  F[1]=f1;
  F[2]=g0;
  F[3]=g1;

  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,false,true,4,1);

  Allocate(count,nshells);
  setcount();

  if(movie) {
    wr.Dimension(Nx+1,2*my-1,(Real *) g1());
    Backward=new fftwpp::crfft2d(Nx+1,2*my-1,g0,(Real *) block);
  }

  InitialCondition=DNS_Vocabulary.NewInitialCondition(ic);
  
  w.Set(Y[OMEGA]);

  Init(T,Y[TRANSFER]);
  Init(T,Y[EK]);
  
  Forcing=DNS_Vocabulary.NewForcing(forcing);

  tcount=0;
  if(restart) {
    Real t0;
    ftin.open(Vocabulary->FileName(dirsep,"t"));
    while(ftin >> t0, ftin.good()) tcount++;
    ftin.close();
  }

  open_output(fprolog,dirsep,"prolog",0);
  out_curve(fprolog,cwrap::kb,"kb",nshells+1);
  out_curve(fprolog,cwrap::kc,"kc",nshells);
  fprolog.close();

  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");

  if(!restart) {
    remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
    remove_dir(Vocabulary->FileName(dirsep,"transfer"));
  }

  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
  mkdir(Vocabulary->FileName(dirsep,"transfer"),0xFFFF);

  errno=0;

  DNSBase::InitialConditions();
  DNSBase::SetParameters();

  if(output)
    open_output(fwk,dirsep,"wk");

  if(movie)
    open_output(fw,dirsep,"w");
}

void DNS::Output(int it)
{
  Real E,Z,P;

  vector f=Y[PAD];
  for(unsigned int i=0; i < my; ++i)
    f[i]=0.0;
  
  vector y=Y[OMEGA];
  
  w.Set(y);
  ComputeInvariants(w,E,Z,P);
  fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;

  if(output) out_curve(fw,y,"w",NY[OMEGA]);

  if(movie)
    OutFrame(it);

  if(spectrum) {
    ostringstream buf;
    Set(T,Y[EK]);
    buf << "ekvk" << dirsep << "t" << tcount;
    const string& s=buf.str();
    open_output(fekvk,dirsep,s.c_str(),0);
    out_curve(fekvk,t,"t");
    out_curve(fekvk,cwrap::Spectrum,"Ek",nshells);
    out_curve(fekvk,cwrap::Dissipation,"nuk*Ek",nshells);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");

    Set(T,Y[TRANSFER]);
    buf.str("");
    buf << "transfer" << dirsep << "t" << tcount;
    const string& S=buf.str();
    open_output(ftransfer,dirsep,S.c_str(),0);
    out_curve(ftransfer,t,"t");
    out_curve(ftransfer,cwrap::Pi,"Pi",nshells);
    out_curve(ftransfer,cwrap::Eta,"Eta",nshells);
    ftransfer.close();
    if(!ftransfer) msg(ERROR,"Cannot write to file transfer");
  }

  tcount++;
  ft << t << endl;

  if(rezero && it % rezero == 0 && spectrum) {
    vector2 Y=Integrator->YVector();
    vector T=Y[TRANSFER];
    for(unsigned i=0; i < nshells; i++)
      T[i]=0.0;
    vector S=Y[EK];
    for(unsigned i=0; i < nshells; i++)
      S[i]=0.0;
  }
}

void DNS::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(w,E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}
