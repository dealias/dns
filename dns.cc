#include "options.h"
#include "kernel.h"
#include "Cartesian.h"
#include "fft.h"
#include "Array.h"

using namespace Array;

const char *method="DNS";
const char *integrator="RK5";

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

// Global Variables;
extern unsigned Nxb,Nxb1,Nyb,Nyp;

// Local Variables;
Real *k2invfactor;

class DNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Direct Numerical Simulation of Turbulence";}
  const char *Abbrev() {return "DNS";}
  DNSVocabulary();
};
   
class DNS : public ProblemBase {
  Real hxinv, hyinv;
  unsigned Nx, Ny;
  unsigned Nxi, Nyi;
  unsigned nmode;
	
  Array3<Real> u,S,f;
  int nspecies;
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
	
  nmode=Nxb*Nyb;
  ny=nmode*nspecies;
  y=new Var[ny];
	
  u.Dimension(Nxb,Nyb,nspecies);
  S.Dimension(Nxb,Nyb,nspecies);
  u.Set(y);
		
  // Initialize arrays with zero boundary conditions
  u=0.0;
	
  for(unsigned i=0; i < Nx; i++) {
    for(unsigned j=0; j < Ny; j++) {
      // Incompressible initial velocity
      Real x=X(i);
      Real y=Y(j);
      Real vx=-cos(twopi*x)*sin(twopi*y);
      Real vy=sin(twopi*x)*cos(twopi*y);
      u(i,j,0)=ICvx*vx;
      u(i,j,1)=ICvy*vy;
    }
  }
	
  k2invfactor=new Real[nmode];
  for(unsigned int k=0; k < nmode; k++) {
//    k2invfactor[k]=1.0/K2();
  }
  
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

void DNS::Output(int it)
{
  Real E;
	
  ComputeInvariants(E);
  fevt << t << "\t" << E << endl;

  if(movie) OutFrame(it);
	
  ft << t << endl;
}

void DNS::OutFrame(int it)
{
  fvx << Nxi << Nyi << 1;
  fvy << Nxi << Nyi << 1;
  
  for(unsigned j=Nyi; j >= 1; j--) {
    for(unsigned i=1; i <= Nxi; i++) {
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
	
  for(unsigned i=1; i <= Nxi; i++) {
    Array2<Real> um=u[i-1], ui=u[i], up=u[i+1];
    for(unsigned j=1; j <= Nyi; j++) {
      for(int s=0; s < nspecies; s++) {
	E += ui(j,s)*ui(j,s);
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
  
  Array2<Real> us(Nxb,2*Nyp);
  Array2<Complex> uk(Nxb,Nyp,(Complex *) us());
  
  Array2<Real> dudx(Nxb,2*Nyp);
  Array2<Complex> ikxu(Nxb,Nyp,(Complex *) dudx());
  
  Array2<Real> dudy(Nxb,2*Nyp);
  Array2<Complex> ikyu(Nxb,Nyp,(Complex *) dudy());
  
  Array3<Real> Ss(nspecies,Nxb,2*Nyp);
  Array3<Complex> Sk(nspecies,Nxb,Nyp,(Complex *) Ss());
  
  for(unsigned s=0; s < nspecies; s++) {
    
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	us(i,j)=u(i,j,s);
      }
    }
  
    rcfft2d(uk,log2Nxb,log2Nyb,-1);
  
    for(unsigned i=0; i < nmode; i++) {
      ikxu(i).re=-uk(i).im*kxmask(i);
      ikxu(i).im=uk(i).re*kxmask(i);
    }
    crfft2d(ikxu,log2Nxb,log2Nyb,1);
  
    for(unsigned i=0; i < nmode; i++) {
      ikyu(i).re=-uk(i).im*kymask(i);
      ikyu(i).im=uk(i).re*kymask(i);
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
    Sk[0]=(kx*miP-Sk[0](i))*Nxybinv;
    Sk[1]=(ky*miP-Sk[1](i))*Nxybinv;
  }
  
  for(unsigned s=0; s < nspecies; s++) {
    crfft2d(Sk[s],log2Nxb,log2Nyb,1);
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	S(i,j,s)=Ss(s,i,j);
      }
    }
  }

}
