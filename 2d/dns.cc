#include "dns.h"

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
Real kH=0.0, kL=0.0;
int pH=1;
int pL=0;
unsigned Nx=1023;
unsigned Ny=1023;
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
Real k0=1.0; // Obsolete
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

class ForcingMask {
  DNS *b;
public: 
  ForcingMask(DNS *b) : b(b) {}
  inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
    if(Forcing->active(i,j)) {
      Complex w=0.0;
      double T=0.0;
      Complex S=0.0;
      Forcing->ForceMask(w,S,T,i,j);
      if(S != 0.0)
        b->fprolog << i << j << abs(S);
    }
  }
};
  
DNS *DNSProblem;

InitialConditionBase *InitialCondition;
ForcingBase *Forcing;

class Zero : public InitialConditionBase {
public:
  const char *Name() {return "Zero";}
  
  Var Value(Real,Real) {return 0.0;}
};

class Constant : public InitialConditionBase {
public:
  const char *Name() {return "Constant";}
  
  Var Value(Real,Real) {return Complex(icalpha,icbeta);}
};

class Equipartition : public InitialConditionBase {
public:
  const char *Name() {return "Equipartition";}

  Var Value(Real kx, Real ky) {
    Real k2=kx*kx+ky*ky;
    Real k=sqrt(k2);
    Real v=icalpha+icbeta*k2;
    v=v ? k*sqrt(2.0/v) : 0.0;
    return randomIC ? v*expi(twopi*drand()) : v*sqrt(0.5)*Complex(1,1);
  }
};

class Benchmark : public InitialConditionBase {
public:
  const char *Name() {return "Benchmark";}

  Var Value(Real kx, Real ky) {
    Real k2=kx*kx+ky*ky;
    Real k=sqrt(k2);
    Real v=icalpha+icbeta*k2;
    v=v ? sqrt(2.0/v) : 0.0;
    return randomIC ? k*v*expi(twopi*drand()) : v*sqrt(0.5)*Complex(k,kx+ky);
  }
};

class Power : public InitialConditionBase {
public:
  const char *Name() {return "Power";}

  Var Value(Real kx, Real ky) {
    Real k2=kx*kx+ky*ky;
    Real v=icbeta*pow(k2,-0.5*icalpha);
    return randomIC ? v*expi(twopi*drand()) : v;
  }
};

// forcing
class None : public ForcingBase {
};

class ConstantBanded : public ForcingBase {
protected:  
  double K1,K2;
public:
  const char *Name() {return "Constant Banded";}
  
  void Init() {
    double h=0.5*deltaf;
    K1=kforce-h;
    K1 *= K1;
    K2=kforce+h;
    K2 *= K2;
  }
  
  bool active(int i, int j) {
    int k=i*i+j*j;
    return K1 < k && k < K2;
  }
  
  void Force(Complex& w, Complex& S, double& T, int i, int j) {
    if(active(i,j)) {
      T += realproduct(force,w);
      S += force;
    }
  }
};

class ConstantList : public ForcingBase {
protected:  
public:
  const char *Name() {return "Constant List";}
  
  int Active(int i, int j) {
    for(int index=0; index < Nforce; ++index)
      if(i == kxforces[index] && j == kyforces[index]) return index;
    return -1;
  }
  
  bool active(int i, int j) {
    return Active(i,j) >= 0;
  }
  
  void Force(Complex& w, Complex& S, double& T, int i, int j) {
    int index=Active(i,j);
    if(index >= 0) {
      Complex force=forces[index];
      T += realproduct(force,w);
      S += force;
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
  
  // For forcing diagnostic
  void ForceMask(Complex& w, Complex& S, double&, int i, int j) {
    if(active(i,j)) {
      S=sqrt(2.0*eta*etanorm);
    }
  }
  
  void ForceStochastic(Complex& w, double& T, int i, int j) {
     if(active(i,j)) {
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
  INITIALCONDITION(Benchmark);
  INITIALCONDITION(Power);

  VOCAB_OBSOLETE(k0,0.0,0.0,"spectral spacing coefficient");
  VOCAB(nuH,0.0,REAL_MAX,"High-wavenumber viscosity");
  VOCAB(nuL,0.0,REAL_MAX,"Low-wavenumber viscosity");
  VOCAB(kL,0.0,STD_MAX,"Restrict low wavenumber dissipation to [1,kL]");
  VOCAB(kH,0.0,STD_MAX,"Restrict high wavenumber dissipation to [kH,infinity)");
  VOCAB(pH,0,0,"Power of Laplacian for high-wavenumber viscosity");
  VOCAB(pL,0,0,"Power of Laplacian for molecular viscosity");

  VOCAB_NOLIMIT(forcing,"Forcing type");
  ForcingTable=new Table<ForcingBase>("forcing");

  VOCAB(eta,0.0,REAL_MAX,"vorticity injection rate");
  VOCAB(force,(Complex) 0.0, (Complex) 0.0,"constant external force");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB_ARRAY(kxforces,"kx force wavenumbers");
  VOCAB_ARRAY(kyforces,"ky force wavenumbers");
  VOCAB_ARRAY(forces,"force amplitudes");
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
class cwrap {
public:
  static Real Spectrum(unsigned i) {return DNSProblem->getSpectrum(i);}
  static Real Dissipation(unsigned i) {return DNSProblem->Dissipation(i);}
  static Real Pi(unsigned i) {return DNSProblem->Pi(i);}
  static Real Eta(unsigned i) {return DNSProblem->Eta(i);}
  static Real kb(unsigned i) {return DNSProblem->kb(i);}
  static Real kc(unsigned i) {return DNSProblem->kc(i);}
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
  kH2=kH*kH;
  kL2=kL*kL;

  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");

  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  
  nshells=spectrum ? (unsigned) (hypot(mx-1,my-1)+0.5) : 0;

  NY[PAD]=my;
  NY[OMEGA]=Nx*my;
  NY[TRANSFER]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;
  cout << "\nALLOCATING FFT BUFFERS" << endl;
  size_t align=sizeof(Complex);

  Allocator(align);

  Dimension(E,nshells);
  Dimension(T,nshells);

  w.Dimension(Nx,my,-mx+1,0);
  S.Dimension(Nx,my,-mx+1,0);

  block=ComplexAlign((Nx+1)*my);
  f0.Dimension(Nx+1,my,-mx,0);
  f1.Dimension(Nx+1,my,block,-mx,0);

  vector f=Y[PAD];
  for(int j=0; j < my; ++j)
    f1(j)=f[j]=0.0;
    
  F[1]=f1;

  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,false,true,2,2);

  Allocate(count,nshells);
  setcount();

  if(movie) {
    wr.Dimension(Nx+1,2*my,(Real *) f1());
    Backward=new fftwpp::crfft2d(Nx+1,2*my-1,f1);
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

  open_output(fprolog,dirsep,"prolog",false);
  
  out_curve(fprolog,cwrap::kb,"kb",nshells+1);
  out_curve(fprolog,cwrap::kc,"kc",nshells);
  
  Loop(InitNone(this),ForcingMask(this));
  
  fprolog.close();

  if(output)
    open_output(fwk,dirsep,"wk");

  if(movie)
    open_output(fw,dirsep,"w");
}

void DNS::Output(int it)
{
  Real E,Z,P;

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
