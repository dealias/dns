#include "options.h"
#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"
#include "fft.h"
#include "Array.h"
#include "Linearity.h"
#include "Forcing.h"
#include "InitialCondition.h"

using namespace Array;

const char *method="DNS";
const char *integrator="RK5";
const char *geometry="Cartesian";
const char *linearity="BandLimited";
const char *forcing="WhiteNoiseBanded";
const char *ic="Equilibrium";

GeometryBase *Geometry;
LinearityBase *Linearity;
ForcingBase *Forcing;
InitialConditionBase *InitialCondition;

Real krmin;
Real krmax;
int reality=1; // Reality condition flag 

// Vocabulary
Real nu=1.0; // Laplacian diffusion
Real rho=1.0;
Real P0=1.0;
Real P1=1.0;
Real ICvx=1.0;
Real ICvy=1.0;
Real xmin=0.0;
Real xmax=1.0;
Real ymin=0.0;
Real ymax=1.0;
Real nonlinear=1.0;

int movie=0;
int pressure=1;

int dumpbins=0;

// Global Variables;
//unsigned int Nmode;
extern unsigned Nxb,Nxb1,Nyb,Nyp;
Real shellmin2;
Real shellmax2;

// Local Variables;
Array1<Real> k2invmask;
Array1<Real> kxmask;
Array1<Real> kymask;

class DNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Direct Numerical Simulation of Turbulence";}
  const char *Abbrev() {return "DNS";}
  DNSVocabulary();
  
  Table<LinearityBase> *LinearityTable;
  Table<GeometryBase> *GeometryTable;
  Table<ForcingBase> *ForcingTable;
  Table<InitialConditionBase> *InitialConditionTable;
  
  GeometryBase *NewGeometry(const char *& key) {
    return GeometryTable->Locate(key);
  }
  LinearityBase *NewLinearity(const char *& key) {
    return LinearityTable->Locate(key);
  }
  ForcingBase *NewForcing(const char *& key) {
    return ForcingTable->Locate(key);
  }
  InitialConditionBase *NewInitialCondition(const char *& key) {
    return InitialConditionTable->Locate(key);
  }
};
   
class DNS : public ProblemBase {
  Real hxinv, hyinv;
  unsigned nmode;
	
  Array3<Real> u,S,f;
  int nspecies;
  
  Array2<Real> us,dudx,dudy;
  Array2<Complex> uk,ikxu,ikyu;
  
  Array3<Real> Ss;
  Array3<Complex> Sk;
  
public:
  DNS();
  virtual ~DNS() {}
  void InitialConditions();
  void Initialize();
  void Output(int it);
  void OutFrame(int it);
  void Source(Var *, Var *, double);
  void ComputeInvariants(Real& E);
  
  Real X(int i) {return xmin+(xmax-xmin)*i/Nxb;}
  Real Y(int i) {return ymin+(ymax-ymin)*i/Nyb;}
};

DNS *DNSProblem;

DNS::DNS()
{
  DNSProblem=this;
}

class CorrectC_PC {
 public:
  void Correct(Real y0, Real y1, Real& y,
	       Real source0, Real source,
	       double dt, double halfdt, int& invertible);
  void Correct(const Complex& y0, const Complex& y1, Complex& y,
	       const Complex& source0, const Complex& source,
	       double dt, double halfdt, int& invertible);
};

class C_PC : public PC, public CorrectC_PC {
 protected:
  int skipmoments;
  int invertible;
 public:
  void Allocate(int n) {
    ny=n; yi=NULL; source=new(n) (Var);
    y=y1=new Var [n]; source0=new(n) (Var); new_y0=1;
    pgrow=0.5/order; pshrink=0.5/(order-1);
    skipmoments=0;
  }
  const char *Name() {return "Conservative Predictor-Corrector";}
	
  void Source(Var *src, Var *Y, double t) {
    DNSProblem->Source(src,Y,t);
  }
	
  void Predictor(double t, double dt, unsigned start,
		 unsigned stop) {
    StandardPredictor(t,dt,start,stop);
  }
	
  int Corrector(double, int, unsigned, unsigned);
};

inline void CorrectC_PC::Correct(Real y0, Real y1, Real& y,
				 Real source0, Real source,
				 double dt, double halfdt, int& invertible)
{
  Real discr=y0*y0+dt*(y0*source0+y1*source);
  if(discr >= 0.0) y=sgn(y1)*sqrt(discr);
  else {
    if(hybrid) y=y0+halfdt*(source0+source);
    else invertible=0;
  }
}

inline void CorrectC_PC::Correct(const Complex& y0, const Complex& y1,
				 Complex& y, const Complex& source0,
				 const Complex& source, double dt,
				 double halfdt, int& invertible)
{
  Correct(y0.re,y1.re,y.re,source0.re,source.re,dt,halfdt,invertible);
  if(!invertible) return;
  Correct(y0.im,y1.im,y.im,source0.im,source.im,dt,halfdt,invertible);
}

int C_PC::Corrector(double dt, int dynamic, unsigned int start,
		    unsigned int stop)
{
  unsigned int j;
  invertible=1;
	
  if(dynamic) {
    for(j=start; j < stop; j++) {
      Var pred=y1[j];
      Correct(y0[j],y1[j],y[j],source0[j],source[j],dt,halfdt,
	      invertible);
      if(!invertible) return 0;
      if(!errmask || errmask[j]) CalcError(y0[j],y[j],pred,y[j]);
    }
  } else for(j=start; j < stop; j++) {
    Correct(y0[j],y1[j],y[j],source0[j],source[j],dt,halfdt,
	    invertible);
    if(!invertible) return 0;
  }		
  return 1;
}

DNSVocabulary::DNSVocabulary()
{
  Vocabulary=this;

  VOCAB(rho,0.0,REAL_MAX,"");
  VOCAB(nu,0.0,REAL_MAX,"");
    
  VOCAB(k0,0.0,STD_MAX,"");
  
  VOCAB(P0,-REAL_MAX,REAL_MAX,"");
  VOCAB(P1,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICvx,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICvy,-REAL_MAX,REAL_MAX,"");

  VOCAB(nonlinear,0.0,1.0,"");
	
  VOCAB(Nx,1,INT_MAX,"");
  VOCAB(Ny,1,INT_MAX,"");
  
  VOCAB(xmin,-REAL_MAX,REAL_MAX,"");
  VOCAB(xmax,-REAL_MAX,REAL_MAX,"");
	
  VOCAB(ymin,-REAL_MAX,REAL_MAX,"");
  VOCAB(ymax,-REAL_MAX,REAL_MAX,"");

  VOCAB(movie,0,1,"");
  VOCAB(pressure,0,1,"");
  
  GeometryTable=new Table<GeometryBase>("geometry");
//  LinearityTable=new Table<LinearityBase>("linearity");
//  ForcingTable=new Table<ForcingBase>("forcing");
//  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  
  BASIS(Cartesian);
  
  METHOD(DNS);
  
  INTEGRATOR(C_PC);
}

DNSVocabulary DNS_Vocabulary;

ofstream ft,fevt,fu;
oxstream fvx,fvy;

void DNS::InitialConditions()
{
  nspecies=2;
  set_fft_parameters();
  
  cout << endl << "GEOMETRY: (" << Nxb << " X " << Nyb << ")" << endl; 
	
  Geometry=DNS_Vocabulary.NewGeometry(geometry);
  Geometry->Create(0);
  nmode=Geometry->nMode();
  
  ny=Nxb*Nyb*nspecies;
  y=new Var[ny];
	
  u.Dimension(Nxb,Nyb,nspecies);
  S.Dimension(Nxb,Nyb,nspecies);
  u.Set(y);
		
  // Initialize arrays with zero boundary conditions
  u=0.0;
	
  for(unsigned i=0; i < Nxb; i++) {
    for(unsigned j=0; j < Nyb; j++) {
      // Incompressible initial velocity
      Real x=X(i);
      Real y=Y(j);
      Real vx=-cos(twopi*x)*sin(twopi*y);
      Real vy=sin(twopi*x)*cos(twopi*y);
      u(i,j,0)=ICvx*vx;
      u(i,j,1)=ICvy*vy;
    }
  }
	
//  cout << u << endl;
  
  us.Allocate(Nxb,2*Nyp);
  uk.Dimension(Nxb,Nyp,(Complex *) us());
  
  dudx.Allocate(Nxb,2*Nyp);
  ikxu.Dimension(Nxb,Nyp,(Complex *) dudx());
  
  dudy.Allocate(Nxb,2*Nyp);
  ikyu.Dimension(Nxb,Nyp,(Complex *) dudy());
  
  Ss.Allocate(nspecies,Nxb,2*Nyp);
  Sk.Dimension(nspecies,Nxb,Nyp,(Complex *) Ss());
  
  
//  BoundaryConditions(u);
  
  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");
  if(output) open_output(fu,dirsep,"u");
  if(movie) {
    open_output(fvx,dirsep,"vx");
    open_output(fvy,dirsep,"vy");
  }
}

void DNS::Initialize()
{
  fevt << "#   t\t\t E\t\t\t Z\t\t\t I\t\t\t C" << endl;
}

void Basis<Cartesian>::Initialize()
{
  cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << " x " << Nyp
       << ")." << endl;
  
  Array1<Complex> temp(nmode);
  Array1<Complex> mask(nfft);
  
  k2invmask.Allocate(nfft);
  kxmask.Allocate(nfft);
  kymask.Allocate(nfft);
  
  for(unsigned int k=0; k < nmode; k++) temp[k]=I*CartesianMode[k].X();
  CartesianPad(mask,temp);
  for(unsigned int k=0; k < nfft; k++) kxmask[k]=mask[k].im;
  
  for(unsigned int k=0; k < nmode; k++) temp[k]=I*CartesianMode[k].Y();
  CartesianPad(mask,temp);
  for(unsigned int k=0; k < nfft; k++) kymask[k]=mask[k].im;
  
  for(unsigned int k=0; k < nfft; k++) {
    Real k2=kxmask[k]*kxmask[k]+kymask[k]*kymask[k];
    k2invmask[k]=k2 ? 1.0/sqrt(k2) : 0.0;
  }
}

void DNS::Output(int it)
{
  Real E;
	
  u.Set(y);
  ComputeInvariants(E);
  fevt << t << "\t" << E << endl;

  if(movie) OutFrame(it);
	
  ft << t << endl;
}

void DNS::OutFrame(int it)
{
  fvx << Nxb << Nyb << 1;
  fvy << Nxb << Nyb << 1;
  
  for(int j=Nyb-1; j >= 0; j--) {
    for(unsigned i=0; i < Nxb; i++) {
      fvx << (float) u(i,j,0);
      fvy << (float) u(i,j,1);
    }
  }
  
  fvx.flush();
  fvy.flush();
}	

void DNS::ComputeInvariants(Real& E)
{
  E=0.0;
	
  for(unsigned i=0; i < Nxb; i++) {
    Array2<Real> ui=u[i];
    for(unsigned j=0; j < Nyb; j++) {
    Array1(Real) uij=ui[j];
      for(int s=0; s < nspecies; s++) {
	E += uij[s]*uij[s];
      }
    }
  }
  
  Real factor=0.5;
  E *= factor;
}

void DNS::Source(Var *source, Var *Y, double)
{
  unsigned nspecies=2;
	
  u.Set(Y);
  S.Set(source);
  
  for(unsigned s=0; s < nspecies; s++) {
    
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	us(i,j)=u(i,j,s);
      }
    }
  
    rcfft2d(uk,log2Nxb,log2Nyb,-1);
  
    for(unsigned i=0; i < nfft; i++) {
      ikxu(i).re=-uk(i).im*kxmask[i];
      ikxu(i).im=uk(i).re*kxmask[i];
    }
    
    crfft2d(ikxu,log2Nxb,log2Nyb,1);
  
    for(unsigned i=0; i < nfft; i++) {
      ikyu(i).re=-uk(i).im*kymask[i];
      ikyu(i).im=uk(i).re*kymask[i];
    }
    crfft2d(ikyu,log2Nxb,log2Nyb,1);
  
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	Ss(s,i,j)=u(i,j,0)*dudx(i,j)+u(i,j,1)*dudy(i,j);
      }
    }
    
    rcfft2d(Sk[s],log2Nxb,log2Nyb,-1);
    
  }
  
  Real Nxybinv2=Nxybinv*Nxybinv;
  
  for(unsigned i=0; i < nfft; i++) {
    Real kx=kxmask[i];
    Real ky=kymask[i];
    if(k2invmask[i] == 0.0) Sk[0](i)=Sk[1](i)=0.0;
    else {
    // Calculate -i*P
      Complex miP=(kx*Sk[0](i)+ky*Sk[1](i))*k2invmask[i];

      Sk[0](i)=(kx*miP-Sk[0](i))*Nxybinv2;
      Sk[1](i)=(ky*miP-Sk[1](i))*Nxybinv2;
    }
  }
  
  for(unsigned s=0; s < nspecies; s++) {
    Array2<Real> Sss=Ss[s];
    Array2<Complex> Sks=Sk[s];
    crfft2d(Sks,log2Nxb,log2Nyb,1);
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	S(i,j,s)=Sss(i,j);
      }
    }
  }
  
}
