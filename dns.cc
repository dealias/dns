#include "options.h"
#include "kernel.h"

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

// Local Variables;
	
class DNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Direct Numerical Simulation of Turbulence";}
  const char *Abbrev() {return "DNS";}
  DNSVocabulary();
};
   
class DNS : public ProblemBase {
  Real hxinv, hyinv;
  int Nx, Ny;
  int Nxi, Nyi;
	
  Array2<U> u,S,f;
  
  Array2<Real> Div;	// Divergence of H
  Array2<Real> P;	// Pressure
  
public:
  DNS() {}
  virtual ~DNS() {}
  void InitialConditions();
  void Initialize();
  void Output(int it);
  void OutFrame(int it);
  void Source(Var *, Var *, double);
  void ComputeInvariants(Real& E);
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
oxstream fvx,fvy,fP;

void DNS::InitialConditions()
{
  Array2<Real> f0;
	
  if(convrate || verbose > 2) defectstat=1;
  if(convrate && niterations == 1) {niterations=100;}
	
  rmin=rmax=0.0;
  imin=imax=0;
	
  cout << endl << "GEOMETRY: (" << Nx << " X " << Ny << ")" << endl; 
	
  ny=Nx*Ny;
  y=new Var[ny];
	
  u.Dimension(Nx,Ny);
  S.Dimension(Nx,Ny);
  u.Set(y);
		
  f.Allocate(Nx,Ny);
  
  P.Allocate(NxP,NyP,oxP,oyP);
  Div.Allocate(NxP,NyP,oxP,oyP);
	
  // Initialize arrays with zero boundary conditions
  u=0.0;
	
  for(int i=0; i < Nx; i++) {
    for(int j=0; j < Ny; j++) {
      // Incompressible initial velocity
      Real x=X(i);
      Real y=Y(j);
      Real vx=-cos(twopi*x)*sin(twopi*y);
      Real vy=sin(twopi*x)*cos(twopi*y);
      u(i,j)=U(ICvx*vx,ICvy*vy);
    }
  }
	
  BoundaryConditions(u);
  
  P=0.0;

  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");
  if(output) open_output(fu,dirsep,"u");
  if(movie) {
    open_output(fvx,dirsep,"vx");
    open_output(fvy,dirsep,"vy");
    open_output(fP,dirsep,"P");
  }
}

void DNS::Initialize()
{
  fevt << "#   t\t\t E\t\t\t Z\t\t\t I\t\t\t C" << endl;
}

void DNS::Output(int it)
{
  Real E,Z,I,C;
	
  ComputeInvariants(E);
  fevt << t << "\t" << E << endl;

  if(movie) OutFrame(it);
	
  ft << t << endl;
}

void DNS::OutFrame(int it)
{
  fvx << Nxi << Nyi << 1;
  fvy << Nxi << Nyi << 1;
  fP << Nxi << Nyi << 1;
  
  for(int j=Nyi; j >= 1; j--) {
    for(int i=1; i <= Nxi; i++) {
      fvx << (float) u(i,j).vx;
      fvy << (float) u(i,j).vy;
      fP << (float) P(i,j);
    }
  }
  
  fvx.flush();
  fvy.flush();
  fP.flush();
}	

void DNS::ComputeInvariants(Real& E)
{
  E=0.0;
	
  Real coeffx=0.5*hxinv;
  Real coeffy=0.5*hyinv;
  for(int i=1; i <= Nxi; i++) {
    Array1(U) um=u[i-1], ui=u[i], up=u[i+1];
    for(int j=1; j <= Nyi; j++) {
      E += ui[j].vx*ui[j].vx+ui[j].vy*ui[j].vy;
    }
  }
  
  Real factor=0.5*volume;
  E *= factor;
}

void PS::FFT(Var *v, Var *u, unsigned int stride)
{
  CartesianPad(v,u,work,stride);
  crfft3d(v,log2Nxb,log2Nyb,log2Nzb,1);
}

void PS::XDerivative(Var *buffer, Var *u, Cartesian *K, unsigned int stride)
{
//#pragma ivdep	
  for(unsigned int i=0; i < nmode; i++) {
    Real k=K[i].X();
    int istride=i*stride;
    buffer[i].re=-u[istride].im*k;
    buffer[i].im=u[istride].re*k;
  }
}

void PS::YDerivative(Var *buffer, Var *u, Cartesian *K, unsigned int stride)
{
//#pragma ivdep	
  for(unsigned int i=0; i < nmode; i++) {
    Real k=K[i].Y();
    int istride=i*stride;
    buffer[i].re=-u[istride].im*k;
    buffer[i].im=u[istride].re*k;
  }
}

void DNS::Source(Var *source, Var *Y, double)
{
  unsigned int i;
  unsigned nspecies=2;
	
  Array3<Real> u(Nxb,Nyb,nspecies);
  Array3<Real> S(Nxb,Nyb,nspecies);
  
  u.Set(Y);
  S.Set(source);
  
  Array2<Real> us(Nxb,2*Nyp);
  Array2<Complex> uk(Nxb,Nyp,us);
  
  Array2<Real> dudx(Nxb,2*Nyp);
  Array2<Complex> ikxu(Nxb,Nyp,dudx);
  
  Array2<Real> dudy(Nxb,2*Nyp);
  Array2<Complex> ikyu(Nxb,Nyp,dudy);
  
  Array3<Real> Ss(nspecies,Nxb,2*Nyp);
  Array3<Complex> Sk(nspecies,Nxb,Nyp,Ss);
  
  for(unsigned s=0; s < nspecies; s++) {
    
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	us(i,j)=u(i,j,s);
      }
    }
  
    rcfft2d(uk,log2Nxb,log2Nyb,-1);
  
    for(i=0; i < nmode; i++) {
      dudx[i].re=-uk[i].im*kxmask[i];
      dudx[i].im=uk[i].re*kxmask[i];
    }
    crfft2d(dudx,log2Nxb,log2Nyb,1);
  
    for(i=0; i < nmode; i++) {
      dudy[i].re=-uk[i].im*kymask[i];
      dudy[i].im=uk[i].re*kymask[i];
    }
    crfft2d(dudy,log2Nxb,log2Nyb,1);
  
  
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	Ss(s,i,j)=u(i,j,0)*dudx(i,j)+u(i,j,1)*dudy(i,j);
      }
    }
    rcfft2d(Sk[s],log2Nxb,log2Nyb,-1);
    
  }
  
  for(i=0; i < nmode; i++) {
    Real kx=CartesianMode[i].X();
    Real ky=CartesianMode[i].Y();
    // Calculate -i*P
    Complex miP=(kx*Sk[0]+ky*Sk[1])*k2invfactor[i];
    Sk[0]=(kx*miP-Sk[0])*Nxybinv;
    Sk[1]=(ky*miP-Sk[1])*Nxybinv;
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
