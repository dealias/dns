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
unsigned int reality=1; // Reality condition flag 

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
unsigned int Nmode;
extern unsigned Nxb,Nxb1,Nyb,Nyp;
Real shellmin2;
Real shellmax2;

// Local Variables;
Real *k2invfactor;
Real *kxmask;
Real *kymask;

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
  DNS() {}
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

class C_PC : public PC {
public:
  const char *Name() {return "Conservative Predictor-Corrector";}
  void Predictor(double, double, unsigned int, unsigned int);
  int Corrector(double, int, unsigned int, unsigned int);
};

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
  Nmode=Geometry->Create(0);
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
  
  k2invfactor=new Real[nmode];
  kxmask=new Real[nmode];
  kymask=new Real[nmode];
  for(unsigned int k=0; k < nmode; k++) {
    k2invfactor[k]=1.0/CartesianMode[k].K2();
    kxmask[k]=CartesianMode[k].X();
    kymask[k]=CartesianMode[k].Y();
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
  
    cout << us << endl;
    
    rcfft2d(uk,log2Nxb,log2Nyb,-1);
  
    for(unsigned i=0; i < nmode; i++) {
      ikxu(i).re=-uk(i).im*kxmask[i];
      ikxu(i).im=uk(i).re*kxmask[i];
    }
    crfft2d(ikxu,log2Nxb,log2Nyb,1);
  
    for(unsigned i=0; i < nmode; i++) {
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
  
  for(unsigned i=0; i < nmode; i++) {
    Real kx=CartesianMode[i].X();
    Real ky=CartesianMode[i].Y();
    // Calculate -i*P
    Complex miP=(kx*Sk[0](i)+ky*Sk[1](i))*k2invfactor[i];
    Sk[0](i)=(kx*miP-Sk[0](i))*Nxybinv;
    Sk[1](i)=(ky*miP-Sk[1](i))*Nxybinv;
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
