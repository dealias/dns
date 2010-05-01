#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "fftw++.h"
#include "Forcing.h"
#include "InitialCondition.h"
#include "convolution.h"
#include "Conservative.h"

#include <sys/stat.h> // On Sun computers this must come after xstream.h

using namespace Array;

const double ProblemVersion=1.0;

const char *method="DNS";
const char *integrator="RK5";
const char *forcing="WhiteNoiseBanded";
const char *ic="Equipartition";
Real icalpha=1.0;
Real icbeta=1.0;

ForcingBase *Forcing;
InitialConditionBase *InitialCondition;

unsigned int nspectrum=1;
// Vocabulary
//Real nu=1.0;
//Real ICvx=1.0;
//Real ICvy=1.0;
unsigned Nx=1;
unsigned Ny=1;
//Real force=1.0;
//Real kforce=1.0;
//Real deltaf=2.0;
int movie=0;
int rezero=0;

int xpad=0;
int ypad=0;

enum Field {OMEGA,EK};

class DNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Direct Numerical Simulation of Turbulence";}
  const char *Abbrev() {return "DNS";}
  DNSVocabulary();
  
  Table<ForcingBase> *ForcingTable;
  Table<InitialConditionBase> *InitialConditionTable;
  
  ForcingBase *NewForcing(const char *& key) {
    return ForcingTable->Locate(key);
  }
  InitialConditionBase *NewInitialCondition(const char *& key) {
    return InitialConditionTable->Locate(key);
  }
};
   
class DNS : public ProblemBase {
  unsigned mx, my; // size of data arrays
  unsigned origin; // linear index of Fourier origin.
  unsigned xorigin; // horizontal index of Fourier origin.

  Real kx0, ky0; // grid spacing factor
  array2<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field;
  
  int tcount;
  array1<unsigned>::opt count;
  
  unsigned nmode;
  unsigned nshells;  // Number of spectral shells
  
  rvector spectrum;
  
  array2<Complex> f0,f1,g0,g1;
  array2<Complex> buffer;
  Complex *block;
  Complex *F[2];
  Complex *G[2];
  
  ImplicitHConvolution2 *Convolution;
  ExplicitHConvolution2 *Padded;
  
public:
  DNS();
  virtual ~DNS() {}
  
  void IndexLimits(unsigned& start, unsigned& stop,
		   unsigned& startT, unsigned& stopT,
		   unsigned& startM, unsigned& stopM) {
    start=Start(OMEGA);
    stop=Stop(OMEGA);
    startT=0;//Start(TRANSFER);
    stopT=0;//Stop(TRANSFER);
    startM=Start(EK);
    stopM=Stop(EK);
  }
  unsigned getNx() {return Nx;}
  unsigned getmx() {return mx;}
  unsigned getmy() {return my;}
  Real getkx0() {return kx0;}
  Real getky0() {return ky0;}
  unsigned getxorigin() {return xorigin;}


  void InitialConditions();
  void Initialize();
  void Output(int it);
  void FinalOutput();
  void OutFrame(int it);
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  
  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    LinearSource(Src,Y,t);
  }
  
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(nspectrum) Spectrum(Src[EK],Y[OMEGA]);
  }
  
  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    NonConservativeSource(Src,Y,t);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    ConservativeSource(Src,Y,t);
    NonConservativeSource(Src,Y,t);
  }
  void ComputeInvariants(Real& E, Real& Z, Real& P);
  void Stochastic(const vector2& Y, double, double);
  
  void Spectrum(vector& S, const vector& y);
};

DNS *DNSProblem;

class Zero : public InitialConditionBase {
public:
  const char *Name() {return "Zero";}
  void Set(Var *w, unsigned n) {
    for(unsigned i=0; i < n; i++) 
      w[i]=0.0;
  }
};

class Constant : public InitialConditionBase {
public:
  const char *Name() {return "Constant";}
  void Set(Var *w, unsigned n) {
    for(unsigned i=0; i < n; i++) 
      w[i]=Complex(icalpha,icbeta);
  }
};

class Equipartition : public InitialConditionBase {
public:
  const char *Name() {return "Equipartition";}
  void Set(Var *w0, unsigned) {
    unsigned Nx=DNSProblem->getNx();
    unsigned my=DNSProblem->getmy();
    unsigned xorigin=DNSProblem->getxorigin();
    Real kx0=DNSProblem->getkx0();
    Real ky0=DNSProblem->getky0();

    array2<Var> w(Nx,my,w0);
    w(xorigin,0)=0;
    for(unsigned i=0; i < Nx; i++) {
      Real kx=kx0*((int) i-(int) xorigin);
      Real kx2=kx*kx;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	Real ky=ky0*j;
	Real k2=kx2+ky*ky;
// Distribute the enstrophy evenly between the real and imaginary components
        Real v=icalpha+icbeta*k2;
        v=v ? sqrt(0.5*k2/v) : 0.0;
	wi[j]=Complex(v,v);
      }
    }
      
// For movie testing:
// w(xorigin,my-1)=2.0;
//    w(mx-1,0)=2.0;
  }
};


DNSVocabulary::DNSVocabulary()
{
  Vocabulary=this;

  //  VOCAB(nu,0.0,REAL_MAX,"Kinematic viscosity");
  //  VOCAB(force,0.0,REAL_MAX,"force coefficient");
  //  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  //  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");
  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  VOCAB(nspectrum,0,1,"Spectrum flag (0=off, 1=on)");

  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");
  VOCAB(icalpha,-REAL_MAX,REAL_MAX,"initial condition parameter");
  VOCAB(icbeta,-REAL_MAX,REAL_MAX,"initial condition parameter");

//  ForcingTable=new Table<ForcingBase>("forcing");
  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  
  METHOD(DNS);

  INITIALCONDITION(Zero);
  INITIALCONDITION(Constant);
  INITIALCONDITION(Equipartition);
}


DNSVocabulary DNS_Vocabulary;

DNS::DNS()
{
  DNSProblem=this;
  check_compatibility(DEBUG);
  ConservativeIntegrators(DNS_Vocabulary.IntegratorTable,this);
}

ifstream ftin;
ofstream ft,fevt,fu;
oxstream fvx,fvy,fw,fekvk;

void DNS::InitialConditions()
{
  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");
  
  kx0=1.0;
  ky0=1.0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;

  xorigin=mx-1;
  origin=xorigin*my;
  nshells=(unsigned) (hypot(mx-1,my-1)+0.5);
  
  NY[OMEGA]=Nx*my;
  NY[EK]=nshells;
  
  size_t align=sizeof(Complex);  
  
  Allocator(align);
  
  w.Dimension(Nx,my);
  f0.Dimension(Nx,my);
  
  unsigned int Nxmy=Nx*my;
  unsigned int nbuf=3*Nxmy;
  if(movie)
    nbuf=max(nbuf,(Nx+xpad)*((Ny+ypad)/2+1));

  block=ComplexAlign(nbuf);
  f1.Dimension(Nx,my,block);
  g0.Dimension(Nx,my,block+Nxmy);
  g1.Dimension(Nx,my,block+2*Nxmy);
  
  F[1]=f1;
  G[0]=g0;
  G[1]=g1;
    
  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl; 

  cout << "\nALLOCATING FFT BUFFERS (" << mx << " x " << my << ")" << endl;
  
  Convolution=new ImplicitHConvolution2(mx,my,2);

  Allocate(count,NY[EK]);
  
  InitialCondition=DNS_Vocabulary.NewInitialCondition(ic);
  // Initialize the vorticity field.
  w.Set(Y[OMEGA]);
  
  if(movie) {
    buffer.Dimension(Nx+xpad,(Ny+ypad)/2+1,block);
    wr.Dimension(Nx+xpad,2*((Ny+ypad)/2+1),(double *) block);
    Padded=new ExplicitHConvolution2(Nx+xpad,Ny+ypad,mx,my,block);
  }
  
  InitialCondition->Set(w,NY[OMEGA]);
  w(xorigin,0)=0;
  HermitianSymmetrizeX(mx,my,xorigin,w);
  
  for(unsigned i=0; i < NY[EK]; i++)
    Y[EK][i]=0.0;
  
  tcount=0;
  if(restart) {
    Real t0;
    ftin.open(Vocabulary->FileName(dirsep,"t"));
    while(ftin >> t0, ftin.good()) tcount++;
    ftin.close();
  }
  
  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");
  
  if(!restart) 
    remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
  
  errno=0;
    
  if(output) 
    open_output(fu,dirsep,"u");
  
  if(movie)
    open_output(fw,dirsep,"w");
}

void DNS::Initialize()
{
  fevt << "#   t\t\t E\t\t\t Z" << endl;
  
  for(unsigned i=0; i < nshells; i++)
    count[i]=0;
  
  for(unsigned i=0; i < Nx; i++) {
    Real kx=kx0*((int) i-(int) xorigin);
    Real kx2=kx*kx;
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real ky=ky0*j;
      Real k2=kx2+ky*ky;
      count[(unsigned)(sqrt(k2)-0.5)]++;
    }
  }
}

void DNS::Spectrum(vector& S, const vector& y)
{
  w.Set(y);
  
  for(unsigned K=0; K < NY[EK]; K++)
    S[K]=0.0;

  // Compute instantaneous angular sum of vorticity squared over circular shell.
		
  for(unsigned i=0; i < Nx; i++) {
    Real kx=kx0*((int) i-(int) xorigin);
    Real kx2=kx*kx;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real ky=ky0*j;
      S[(unsigned)(sqrt(kx2+ky*ky)-0.5)] += abs2(wi[j]);
    }
  }
}

void DNS::Output(int it)
{
  Real E,Z,P;
	
  w.Set(y);
  ComputeInvariants(E,Z,P);
  fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;

  Var *y=Y[0];
  if(output) out_curve(fu,y,"u",NY[0]);
  
  if(movie) OutFrame(it);
	
  if(nspectrum) {
    ostringstream buf;
    buf << "ekvk" << dirsep << "t" << tcount;
    open_output(fekvk,dirsep,buf.str().c_str(),0);
    out_curve(fekvk,t,"t");
    Var *y1=Y[EK];
    fekvk << NY[EK];
    for(unsigned K=0; K < NY[EK]; K++)
      fekvk << y1[K]*twopi/(K*count[K]);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");
  }    

  tcount++;
  ft << t << endl;
  
  if(rezero && it % rezero == 0) {
    vector2 Y=Integrator->YVector();
    vector T=Y[EK];
    for(unsigned i=0; i < NY[EK]; i++)
      T[i]=0.0;
  }
}

void DNS::OutFrame(int)
{
  w.Set(Y[OMEGA]);
  unsigned int offset=(Nx+xpad)/2-mx+1;
  unsigned int stop=2*mx-1;
  for(unsigned int i=0; i < stop; ++i) {
    unsigned int I=i+offset;
    for(unsigned int j=0; j < my; j++)
      buffer(I,j)=w(i,j);
  }
    
  Padded->pad(buffer);
  Padded->backwards(buffer);

  fw << 1 << (Ny+ypad) << (Nx+xpad);

  for(int j=Ny+ypad-1; j >= 0; j--) {
    for(unsigned i=0; i < (Nx+xpad); i++) {
      fw << (float) wr(i,j);
    }
  }
  
  fw.flush();
}	

void DNS::ComputeInvariants(Real& E, Real& Z, Real& P)
{
  E=Z=P=0.0;
  w.Set(Y[OMEGA]);
  for(unsigned i=0; i < Nx; i++) {
    Real kx=kx0*((int) i-(int) xorigin);
    Real kx2=kx*kx;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real w2=abs2(wi[j]);
      Z += w2;
      Real ky=ky0*j;
      Real k2=kx2+ky*ky;
      E += w2/k2;
      P += k2*w2;
    }
  }
}

void DNS::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

void DNS::LinearSource(const vector2& Src, const vector2& Y, double)
{
}

void DNS::NonLinearSource(const vector2& Src, const vector2& Y, double)
{
  w.Set(Y[OMEGA]);
  f0.Set(Src[OMEGA]);
 
  f0(origin)=0.0;
  f1(origin)=0.0;
  g0(origin)=0.0;
  g1(origin)=0.0;
  
  for(unsigned i=0; i < Nx; ++i) {
    Real kx=kx0*((int) i-(int) xorigin);
    Real kx2=kx*kx;
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    vector g0i=g0[i];
    vector g1i=g1[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real ky=ky0*j;
      Complex wij=wi[j];
      Complex kxw=Complex(-kx*wij.im,kx*wij.re);
      Complex kyw=Complex(-ky*wij.im,ky*wij.re);
      f0i[j]=kxw;
      f1i[j]=kyw;
      Real k2inv=1.0/(kx2+ky*ky);
      g0i[j]=k2inv*kyw;
      g1i[j]=-k2inv*kxw;
    }
  }
  
  F[0]=f0;
  Convolution->convolve(F,G);
  
#if 0
  double sum=0.0;
  for(unsigned i=0; i < Nx; ++i) {
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Complex wij=wi[j];
      sum += (f0[i][j]*conj(wij)).re;
    }
  }
  
  cout << sum << endl;
#endif  
}

void DNS::Stochastic(const vector2&Y, double, double)
{
//  u.Set(Y[OMEGA]);
//  Real factor=sqrt(2.0*dt)*rand_gauss();
}
