#include "dns.h"

const double ProblemVersion=1.0;

#if !COMPLEX
#error DNS requires COMPLEX=1 in options.h
#endif

using namespace utils;
using namespace parallel;

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
uInt Nx=1023;
uInt Ny=1023;
Real eta=0.0;
Complex force=0.0;
Real kforce=1.0;
const uInt Nforce=2;
Complex forces[Nforce];
Int kxforces[Nforce];
Int kyforces[Nforce];
Real deltaf=1.0;
uInt movie=0;
uInt rezero=0;
uInt spectrum=1;
uInt modalenergies=0;
Real icalpha=1.0;
Real icbeta=1.0;
Real k0=1.0; // Obsolete
uInt randomIC=0;

// This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
// requires only 4 FFTs per stage.
void multadvection2(Complex **F, uInt n, Indices *, uInt threads)
{
  double* F0=(double *) F[0];
  double* F1=(double *) F[1];

  PARALLELIF(
    n > threshold,
    for(uInt j=0; j < n; ++j) {
      double u=F0[j];
      double v=F1[j];
      F0[j]=v*v-u*u;
      F1[j]=u*v;
    });
}

Real *Triplet, *Norm1, *Norm2;

// A=8, B=0 TODO: Reduce to A=7, B=0 using incompressibility
void multTriplet(Complex **F, uInt n, Indices *, uInt threads)
{
  double* F0=(double *) F[0];
  double* F1=(double *) F[1];
  double* F2=(double *) F[2];
  double* F3=(double *) F[3];
  double* F4=(double *) F[4];
  double* F5=(double *) F[5];
  double* F6=(double *) F[6];
  double* F7=(double *) F[7];

  Real triplet=0.0, norm1=0.0, norm2=0.0;
  PARALLELIF(
    n > threshold,
    for(unsigned int j=0; j < n; ++j) {
      double u=F0[j];
      double v=F1[j];
      double ux=F2[j];
      double uy=F3[j];
      double vx=F4[j];
      double vy=F5[j];
      double A2u=F6[j];
      double A2v=F7[j];
      double bx=u*ux+v*uy;
      double by=u*vx+v*vy;

      triplet += bx*A2u+by*A2v;
      norm1 += bx*bx+by*by;
      norm2 += A2u*A2u+A2v*A2v;
    });
  size_t t=parallel::get_thread_num();
  Triplet[t] += triplet;
  Norm1[t] += norm1;
  Norm2[t] += norm2;
}

class DNS : public DNSBase, public ProblemBase {
public:
  DNS();
  ~DNS();

  void InitialConditions();

  void Output(uInt);
  void FinalOutput();
  oxstream fprolog;

  void IndexLimits(uInt& start, uInt& stop,
		   uInt& startT, uInt& stopT,
		   uInt& startM, uInt& stopM) {
    start=Start(OMEGA);
    stop=Stop(OMEGA);
    startT=Start(TRANSFERE);
    stopT=Stop(DISSIPATIONZ);
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

  void ZeroDiagnostics();
};

class ForcingMask {
  DNS *b;
public:
  ForcingMask(DNS *b) : b(b) {}
  inline void operator()(const vector&, const vector&,
                         const nuvector &, Int i, Int j, size_t) {
    if(Forcing->active(i,j)) {
      Complex w=0.0;
      Complex S=0.0;
      Forcing->ForceMask(w,S,i,j);
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

  bool active(Int i, Int j) {
    Int k=i*i+j*j;
    return K1 < k && k < K2;
  }

  double Force(Complex& w, Complex& S, Int i, Int j) {
    if(active(i,j)) {
      S += force;
      return realProduct(force,w);
    }
    return 0.0;
  }
};

class ConstantList : public ForcingBase {
protected:
public:
  const char *Name() {return "Constant List";}

  Int Active(Int i, Int j) {
    for(uInt index=0; index < Nforce; ++index)
      if(i == kxforces[index] && j == kyforces[index]) return index;
    return -1;
  }

  bool active(Int i, Int j) {
    return Active(i,j) >= 0;
  }

  double Force(Complex& w, Complex& S, Int i, Int j) {
    Int index=Active(i,j);
    if(index >= 0) {
      Complex force=forces[index];
      S += force;
      return realProduct(force,w);
    }
    return 0.0;
  }
};

class WhiteNoiseBanded : public ConstantBanded {
  Complex f0;
  Real etanorm;
public:
  const char *Name() {return "White-Noise Banded";}

  void Init(uInt fcount) {
    etanorm=1.0/((Real) fcount);
  }

  bool Stochastic(double dt) {
    f0=sqrt(2.0*dt*eta*etanorm);
    return true;
  }

  // For forcing diagnostic
  void ForceMask(Complex& w, Complex& S, Int i, Int j) {
    if(active(i,j))
      S=sqrt(2.0*eta*etanorm);
  }

  double ForceStochastic(Complex& w, Int i, Int j) {
    if(active(i,j)) {
      Complex f=f0*crand_gauss();
      double eta=realProduct(f,w)+0.5*abs2(f);
      w += f;
      return eta;
    }
    return 0.0;
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
  VOCAB(Nx,1,SIZE_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,SIZE_MAX,"Number of dealiased modes in y direction");
  VOCAB(movie,0,1,"Output movie? (0=no, 1=yes)");
  VOCAB(spectrum,0,1,"Output spectrum? (0=no, 1=yes)");
  VOCAB(modalenergies,0,1,"Output modal energies? (0=no, 1=yes)");
  VOCAB(rezero,0,SIZE_MAX,"Rezero moments every rezero output steps for high accuracy");
  VOCAB_NODUMP(threshold,0,SIZE_MAX,"Multithreading threshold");

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
  deleteAlign(block[0]);
  delete [] block;
}

// wrapper for outcurve routines
class cwrap {
public:
  static Real Spectrum(uInt i) {return DNSProblem->getSpectrum(i);}
  static Real TE(uInt i) {return DNSProblem->TE_(i);}
  static Real TZ(uInt i) {return DNSProblem->TZ_(i);}
  static Real Eps(uInt i) {return DNSProblem->Eps_(i);}
  static Real Eta(uInt i) {return DNSProblem->Eta_(i);}
  static Real Zeta(uInt i) {return DNSProblem->Zeta_(i);}
  static Real DE(uInt i) {return DNSProblem->DE_(i);}
  static Real DZ(uInt i) {return DNSProblem->DZ_(i);}

  static Real kb(uInt i) {return DNSProblem->kb(i);}
  static Real kc(uInt i) {return DNSProblem->kc(i);}
};

void DNS::Initialize()
{
  DNSBase::Initialize();
}

void DNS::InitialConditions()
{
  fftw::maxthreads=threads;
  if(threshold < SIZE_MAX)
    lastThreads=0;
  Threshold(threads);

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

  nshells=spectrum ? (uInt) (hypot(mx-1,my-1)+0.5) : 0;

  NY[OMEGA]=mx*(my-1)+(mx-1)*my;
  NY[TRANSFERE]=nshells;
  NY[TRANSFERZ]=nshells;
  NY[EPS]=nshells;
  NY[ETA]=nshells;
  NY[ZETA]=nshells;
  NY[DISSIPATIONE]=nshells;
  NY[DISSIPATIONZ]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;
  cout << "\nALLOCATING FFT BUFFERS" << endl;
  uInt align=utils::ALIGNMENT;

  Allocator(align);

  Dimension(Sum,nshells);

  uInt Nshells=threads*nshells;
  Allocate(TE,Nshells);
  Allocate(TZ,Nshells);
  Allocate(Eps,Nshells);
  Allocate(Eta,Nshells);
  Allocate(Zeta,Nshells);
  Allocate(DE,Nshells);
  Allocate(DZ,Nshells);
  Allocate(E,Nshells);

  Allocate(enstrophy,threads);
  Allocate(energy,threads);
  Allocate(palinstrophy,threads);

  Dimension(Y0,Y);
  Set(Y0,Y);

  w.Dimension(mx,my);
  S.Dimension(mx,my);

  uInt Mx=3*mx-2;
  uInt My=3*my-2;

  my1=utils::align(my+1);

  uInt A=2, B=2; // 2 inputs, 2 outputs in Basdevant scheme
  uInt N=max(A,B);

// Disable overwrite optimization to allow indexing transformed values.
  Application appx(A,B,multNone,threads,false);
//  Application appx(A,B,multNone,threads,1024,1,1);
//  Application appx(A,B,multNone,threads,8192,1,1);
  auto fftx=new fftPadCentered(Nx+1,Mx,appx,my,my1);

  Application appy(A,B,multadvection2,appx);
//  Application appy(A,B,multadvection2,appx,3072,1,0);
//  Application appy(A,B,multadvection2,appx,24576,1,0);
  auto ffty=new fftPadHermitian(Ny,My,appy);

  Convolve=new Convolution2(fftx,ffty);

  saveWisdom();

  block=ComplexAlign(N,fftx->inputLength());

  u.Dimension(Nx+1,my1,block[0],-mx,0);
  v.Dimension(Nx+1,my1,block[1],-mx,0);

  F[0]=u;
  F[1]=v;

  Allocate(count,nshells);

  if(movie) {
    uInt A=8, B=0; // 8 inputs, 0 outputs in Basdevant scheme

    uInt nx=4*mx-3;
    uInt ny0=4*my-3;

    Triplet=new Real[threads];
    Norm1=new Real[threads];
    Norm2=new Real[threads];

    Application appX(A,B,multNone,threads);
    auto fftX=new fftPadCentered(Nx+1,nx,appX,my,my1);

    Application appY(A,B,multTriplet,appX);
    auto fftY=new fftPadHermitian(Ny,ny0,appY);

    Convolve2=new Convolution2(fftX,fftY);

    saveWisdom();

    ux.Allocate(Nx+1,my1,-mx,0,align);
    uy.Allocate(Nx+1,my1,-mx,0,align);
    vx.Allocate(Nx+1,my1,-mx,0,align);
    vy.Allocate(Nx+1,my1,-mx,0,align);
    A2u.Allocate(Nx+1,my1,-mx,0,align);
    A2v.Allocate(Nx+1,my1,-mx,0,align);

    F[2]=ux;
    F[3]=uy;
    F[4]=vx;
    F[5]=vy;
    F[6]=A2u;
    F[7]=A2v;

    W.Dimension(Nx+1,my,v(),-mx,0);
    wr.Dimension(Nx+1,2*my,(Real *) v());
    Backward=new fftwpp::crfft2d(Nx+1,2*my-1,v);
  }

  InitialCondition=DNS_Vocabulary.NewInitialCondition(ic);

  w.Set(Y[OMEGA]);

  ZeroDiagnostics();

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

  Loop(InitNone(this),ForcingMask(this),1);

  fprolog.close();

  if(modalenergies)
    open_output(fek,dirsep,"ek");

  if(movie)
    open_output(fw,dirsep,"w");
}

void DNS::ZeroDiagnostics()
{
  Zero(Y[TRANSFERE]);
  Zero(Y[TRANSFERZ]);
  Zero(Y[EPS]);
  Zero(Y[ETA]);
  Zero(Y[ZETA]);
  Zero(Y[DISSIPATIONE]);
  Zero(Y[DISSIPATIONZ]);
  Zero(Y[EK]);
}

void DNS::Output(uInt it)
{
  vector y=Y[OMEGA];
  w.Set(y);

  ComputeInvariants();
  fevt << t << "\t" << Energy << "\t" << Enstrophy << "\t" << Palinstrophy << endl;

  if(output) out_curve(fw,y,"w",NY[OMEGA]);

  if(movie)
    OutFrame(it);

  if(modalenergies)
    OutEnergies();

  if(spectrum) {
    ostringstream buf;
    buf << "ekvk" << dirsep << "t" << tcount;
    const string& s=buf.str();
    open_output(fekvk,dirsep,s.c_str(),0);
    out_curve(fekvk,t,"t");
    out_curve(fekvk,cwrap::Spectrum,"Ek",nshells);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");

    buf.str("");
    buf << "transfer" << dirsep << "t" << tcount;
    const string& S=buf.str();
    open_output(ftransfer,dirsep,S.c_str(),0);
    out_curve(ftransfer,t,"t");
    out_curve(ftransfer,cwrap::TE,"TE",nshells);
    out_curve(ftransfer,cwrap::TZ,"TZ",nshells);
    out_curve(ftransfer,cwrap::Eps,"eps",nshells);
    out_curve(ftransfer,cwrap::Eta,"eta",nshells);
    out_curve(ftransfer,cwrap::Zeta,"zeta",nshells);
    out_curve(ftransfer,cwrap::DE,"DE",nshells);
    out_curve(ftransfer,cwrap::DZ,"DZ",nshells);
    ftransfer.close();
    if(!ftransfer) msg(ERROR,"Cannot write to file transfer");
  }

  tcount++;
  ft << t << endl;

  if(rezero && it % rezero == 0 && spectrum) {
    vector2 Y=Integrator->YVector();
    ZeroDiagnostics();
  }
}

void DNS::FinalOutput()
{
  ComputeInvariants();
  cout << endl;
  cout << "Energy = " << Energy << newl;
  cout << "Enstrophy = " << Enstrophy << newl;
  cout << "Palinstrophy = " << Palinstrophy << newl;
}
