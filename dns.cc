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

//new parameters: need more test

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


void DNS::Source(Var *source, Var *Y, double)
{
  u.Set(Y);
  BoundaryConditions(u);
  S.Set(source);
	
  Real coeffx=0.5*hxinv;
  Real coeffy=0.5*hyinv;

// The following code is used to calculate fluxes without hyperviscosity.

  for(int i=1; i <= Nxi; i++) {
    Array1(U) Gxi=Gx[i], Gyi=Gy[i];
    Array1(Real) Exi=Ex[i], Eyi=Ey[i];
    Array1(U) ui=u[i];

// Calculate the fluxes......

    for(int j=1; j <= Nyi; j++) {

      // Flux Gx=vx * u:
      Gxi[j].vx=ui[j].vx*ui[j].vx;
      Gxi[j].vy=ui[j].vx*ui[j].vy;
      Gxi[j].C=(ui[j].vx+Exi[j])*ui[j].C;
      // Flux Gy=vy * u:
      Gyi[j].vx=ui[j].vy*ui[j].vx;
      Gyi[j].vy=ui[j].vy*ui[j].vy;
      Gyi[j].C=(ui[j].vy+Eyi[j])*ui[j].C;
    }
  }
		
  BoundaryConditions(Gx);
  BoundaryConditions(Gy);

  Real coeffx0=coeffx*nonlinear;
  Real coeffy0=coeffy*nonlinear;
  
  for(int i=1; i <= Nxi; i++) {
    Array1(U) Gxm=Gx[i-1], Gxp=Gx[i+1];
    Array1(U) Gyi=Gy[i];
    
    Array1(Real) Hxi=Hx[i], Hyi=Hy[i];
    Array1(Real) Fxi=Fx[i], Fyi=Fy[i];
    Array1(U) Si=S[i];

// center in space for momentum equation; optionally upwind C equation

    for(int j=1; j <= Nyi; j++) {
      // div (v u)
      U advection=-coeffx0*(Gxm[j]-Gxp[j])-coeffy0*(Gyi[j-1]-Gyi[j+1]);
      
      // H=F - div (v u)
      Hxi[j]=Fxi[j]-advection.vx;
      Hyi[j]=Fyi[j]-advection.vy;
      
    }
  }

  if(pressure) {
    for(int i=1; i <= Nxi; i++) {
      Array1(Real) Hxm=Hx[i-1], Hxp=Hx[i+1];
      Array1(Real) Hyi=Hy[i];
      Array1(Real) Di=Div[i];

      for(int j=1; j <= Nyi; j++) {
	// Div=div H
	Di[j]=-coeffx*(Hxm[j]-Hxp[j])-coeffy*(Hyi[j-1]-Hyi[j+1]);
      }
    }
  }
	
  Real rhoinv=1.0/rho;
  Real coeffrx=coeffx*rhoinv;
  Real coeffry=coeffy*rhoinv;
	
  for(int i=1; i <= Nxi; i++) {
    Array1(Real) Pm=P[i-1];
    Array1(Real) Pi=P[i];
    Array1(Real) Pp=P[i+1];
    Array1(U) Si=S[i];
    Array1(Real) Hxi=Hx[i];
    Array1(Real) Hyi=Hy[i];
    for(int j=1; j <= Nyi; j++) {
      // H - grad P
      Si[j].vx=Hxi[j]+coeffrx*(Pm[j]-Pp[j]);
      Si[j].vy=Hyi[j]+coeffry*(Pi[j-1]-Pi[j+1]);
    }
  }
}
