#include "options.h"
#include "kernel.h"
#include "Grid2.h"
#include "Poisson.h"
#include "Helmholtz.h"
#include "Parcel.h"

using namespace Array;

const char *method="OS";
const char *integrator="PC";

#define PERIODIC 0
#define ZERODIV 1
#define LAGRANGIAN 1
#define UPWIND 0

// Vocabulary
Real kappa2=1.0;
Real nuv=1.0, nuC=1.0; // Laplacian diffusion
Real alpha=1.0;
Real P0=1.0;
Real P1=1.0;
Real rho=1.0;
Real mu=0.0;
Real ICvx=1.0;
Real ICvy=1.0;
Real ICC=1.0;
Real xmin=0.0;
Real xmax=1.0;
Real ymin=0.0;
Real ymax=1.0;
Real implicit_factor=0.5;
Real explicit_factor=0.5;
Real nonlinear=1.0;

Real delta=0.0; //OBSOLETE

//new parameters: need more test

Real zeta=1.0;
Real Phi0=1.0;
Real Phi1=0.0;

int xlevel=4;
int ylevel=4;
int exactpsi=1;
int niterations=1;
int ncrank=1;
int nfirst=1000;
int npresmooth=0;
int igamma=1;
int npostsmooth=1;
int defectstat=0;
Real convrate=0.0;
int movie=0;
int pressure=1;

// Local Variables;
static U beta;
static const Real twelfth=1.0/12.0;
static int nlevel;
static U rmin=REAL_MAX, rmax=0.0;
static U urmin=REAL_MAX, urmax=0.0;
static U nu;
static int imin=INT_MAX, imax=0;
static int uimin=INT_MAX, uimax=0;
	
class OSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "ELECTRO-OSMOSIS MULTIGRID";}
  const char *Abbrev() {return "OS";}
  OSVocabulary();
};
   
class OS : public ProblemBase {
  Real hxinv, hyinv, volume;
  int Nx, Ny;
  int Nxi, Nyi;
	
  Array2<U> u,S,f;
  Array2<U> Gx; 	// Flux
  Array2<U> Gy;

  Array2<Real> Ex; 	// mu * grad(phi) , dimensionless
  Array2<Real> Ey;
  
  Grid2<Real> *gridP;
  Grid2<Var> *gridCrank;
  Array2<Real> phi,psi;
  Array2<Real> Fx;	// Body force
  Array2<Real> Fy;
  Array2<Real> Hx;	// Body force minus advection
  Array2<Real> Hy;

  Array2<Real> Div;	// Divergence of H
  Array2<Real> P;	// Pressure
  
  Array2<Cell> grid;
  Array2<U> u0,uL;
  DynVector<Parcel> parcel;
  DynVector<unsigned int> Inactive;

  Real defect0,defect;
  U udefect0,udefect;
  int first;
public:
  OS() {}
  virtual ~OS() {}
  void InitialConditions();
  void Initialize();
  void Setup();
  void Output(int it);
  void OutFrame(int it);
  void Source(Var *, Var *, double);
  void ComputeInvariants(Real& E, Real& Z, Real& I, Real& C);
  void Transform(Var *Y, double, double dt, Var *&yi);
  void BackTransform(Var *Y, double, double dt, Var *yi);
  template<class T>
  inline void BoundaryConditions(const Array2<T>&);
  inline void BoundaryConditions(const Array2<Real>&);
};

class C_PC : public PC {
public:
  const char *Name() {return "Conservative Predictor-Corrector";}
  void Predictor(double, double, unsigned int, unsigned int);
  int Corrector(double, int, unsigned int, unsigned int);
};

OSVocabulary::OSVocabulary()
{
  Vocabulary=this;

  VOCAB(mu,0.0,REAL_MAX,"");
  VOCAB(zeta,0.0,REAL_MAX,"");
  VOCAB(Phi0,-REAL_MAX,REAL_MAX,"");
  VOCAB(Phi1,-REAL_MAX,REAL_MAX,"");
    
  VOCAB(nuv,0.0,REAL_MAX,"");
  VOCAB(nuC,0.0,REAL_MAX,"");
  VOCAB(kappa2,0.0,REAL_MAX,"");
  VOCAB(alpha,-REAL_MAX,REAL_MAX,"");

  VOCAB_NODUMP(delta,-REAL_MAX,REAL_MAX,""); // Obsolete
  
  VOCAB(P0,-REAL_MAX,REAL_MAX,"");
  VOCAB(P1,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICvx,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICvy,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICC,-REAL_MAX,REAL_MAX,"");
  VOCAB(nonlinear,0.0,1.0,"");
  VOCAB(explicit_factor,0.0,REAL_MAX,"");
  VOCAB(implicit_factor,0.0,REAL_MAX,"");
	
  VOCAB(xmin,-REAL_MAX,REAL_MAX,"");
  VOCAB(xmax,-REAL_MAX,REAL_MAX,"");
	
  VOCAB(ymin,-REAL_MAX,REAL_MAX,"");
  VOCAB(ymax,-REAL_MAX,REAL_MAX,"");

  VOCAB(xlevel,1,INT_MAX,"");
  VOCAB(ylevel,1,INT_MAX,"");

  VOCAB(exactpsi,1,1,"");
  
  VOCAB(niterations,1,INT_MAX,"");
  VOCAB(ncrank,1,INT_MAX,"");
  VOCAB(nfirst,1,INT_MAX,"");
  VOCAB(npresmooth,0,INT_MAX,"");
  VOCAB(igamma,1,INT_MAX,"");
  VOCAB(npostsmooth,0,INT_MAX,"");
  VOCAB(defectstat,0,1,"");
  VOCAB(convrate,0.0,REAL_MAX,"");
  VOCAB(movie,0,1,"");
  VOCAB(pressure,0,1,"");
  
  METHOD(OS);
}

OSVocabulary OS_Vocabulary;

ofstream ft,fevt,fdefect,fudefect,fu,ftest,fc,fvel;
oxstream fvx,fvy,fC,fP,fphi,fpsi;

Limits xlimphi,ylimphi;
Limits xlimpsi,ylimpsi;
Limits xlimP,ylimP;
Limits xlimP0,ylimP0;
Limits xlimu,ylimu;

template<class T>
class Phi : public Poisson2<T> {
public:	
  Limits XMeshRange() {return xlimphi;}
  Limits YMeshRange() {return ylimphi;}
  void BoundaryConditions(const Array2<T>& u) {
#if PERIODIC
    XPeriodic(u);
    YPeriodic(u);
#else
    XDirichlet2(u,Phi0,Phi1);    
    YNeumann(u);
#endif
  }
};

template<class T>
class Psi : public Helmholtz2<T> {
public:	
  T A() {return 1.0;}
  T B() {return -kappa2;}
  Limits XMeshRange() {return xlimpsi;}
  Limits YMeshRange() {return ylimpsi;}
	
  void BoundaryConditions(const Array2<T>& u) {
#if PERIODIC
    XPeriodic(u);
    YPeriodic(u);
#else
    XNeumann(u);
    YDirichlet2(u,-zeta,-zeta);    
#endif
  }
};

template<class T>
class Crank : public Helmholtz2<T> {
public:	
  T A() {return beta;}
  T B() {return 1.0;}
  Limits XMeshRange() {return xlimu;}
  Limits YMeshRange() {return ylimu;}
	
  void BoundaryConditions(const Array2<T>& u) {
#if PERIODIC
    XPeriodic(u);
    YPeriodic(u);
#else
    Array1(T) u0=u[i1-1], u2=u[i1+1];
    Array1(T) unxm1=u[i2-1], unx1=u[i2+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      // xmin: Neumann BC on v, Constant BC on C
      u0[j].vx=u2[j].vx;
      u0[j].vy=u2[j].vy;
      u0[j].C=0.0;
      // xmax: Neumann BC on v, Constant BC on C
      unx1[j].vx=unxm1[j].vx;
      unx1[j].vy=unxm1[j].vy;
      unx1[j].C=0.0;
    }
    
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1(T) ui=u[i];
      // ymin: Dirichlet BC on v, Neumann BC on C
      ui[j1-2].vx=ui[j1-1].vx=0.0;
      ui[j1-2].vy=ui[j1-1].vy=0.0;
      ui[j1-2].C=ui[j1-1].C=ui[j1].C;
      // ymax: Dirichlet BC on v, Neumann BC on C
      ui[j2+2].vx=ui[j2+1].vx=0.0;
      ui[j2+2].vy=ui[j2+1].vy=0.0;
      ui[j2+2].C=ui[j2+1].C=ui[j2].C;
    }
#endif
  }
};

template<class T>
inline void OS::BoundaryConditions(const Array2<T>& H)
{
  gridCrank->BoundaryConditions(H);
}

// Neumann boundary conditions on fluxes at all boundaries
inline void OS::BoundaryConditions(const Array2<Real>& H)
{
#if PERIODIC
  // Periodic BC's
  Array1(Real) H0=H[0], H1=H[1];
  Array1(Real) HN=H[Nxi], HNp1=H[Nxi+1];
  for(int j=0; j < Ny; j++) {
    H0[j]=HN[j];
    HNp1[j]=H[1][j];
  }
  for(int i=0; i < Nx; i++) {
    Array1(Real) Hi=H[i];
    Hi[0]=Hi[Nyi];
    Hi[Nyi+1]=Hi[1];
  }
#else
  Array1(Real) H0=H[0], H2=H[2];
  Array1(Real) HNm1=H[Nxi-1], HNp1=H[Nxi+1];
  for(int j=0; j < Ny; j++) {
    // xmin: Neumann BC
    H0[j]=H2[j];
    // xmax: Neumann BC
    HNp1[j]=HNm1[j];
  }
  for(int i=0; i < Nx; i++) {
    Array1(Real) Hi=H[i];
    // ymin: Neumann BC
    Hi[0]=Hi[2];
    // ymax: Neumann BC
    Hi[Nyi+1]=Hi[Nyi-1];
  }
#endif
}

template<class T>
class Pressure0 : public Poisson2<T> {
public:	
  Limits XMeshRange() {return xlimP0;}
  Limits YMeshRange() {return ylimP0;}
  void BoundaryConditions(const Array2<T>& u) {
#if PERIODIC
    XPeriodic(u);
    YPeriodic(u);
#else
    XDirichlet2(u,P0,P1);
    YNeumann(u);
#endif
  }
};

template<class T>
class Pressure : public Poisson2<T> {
public:	
  Limits XMeshRange() {return xlimP;}
  Limits YMeshRange() {return ylimP;}
  void Defect(const Array2<T>& d0, const Array2<T>& u, const Array2<T>& f);
  void GaussSeidel(const Array2<T>&, const Array2<T>&, int, int, int, int);
  void Smooth(const Array2<T>& u, const Array2<T>& f);
  void Restrict(const Array2<T>& r, const Array2<T>& u);
	
  void BoundaryConditions(const Array2<T>& u) {
#if PERIODIC
    XPeriodic2(u);
    YPeriodic2(u);
#else
    XDirichlet2(u,P0,P1);
    YNeumann2(u);
#endif
  };
};

template<class T>
void Pressure<T>::Defect(const Array2<T>& d0, const Array2<T>& u,
			 const Array2<T>& f)
{
  T coeffx2y2=-2.0*(hx2inv+hy2inv);
	
  for(int i=i1; i <= i2; i++) {
    Array1(T) di=d0[i];
    Array1(T) fi=f[i];
    Array1(T) uM=u[i-2], ui=u[i], uP=u[i+2];
    for(int j=j1; j <= j2; j++) {
      di[j]=0.25*(hx2inv*(uM[j]+uP[j])+coeffx2y2*ui[j]+
		  hy2inv*(ui[j-2]+ui[j+2]))-fi[j];
    }
  }
}

template<class T>
void Pressure<T>::GaussSeidel(const Array2<T>& u, const Array2<T>& f,
			      int i0, int j0, int istep, int jstep)
{
  T factor=0.5/(hx2inv+hy2inv);
  for(int i=i1+i0; i <= i2; i += istep) {
    Array1(T) fi=f[i];
    Array1(T) uM=u[i-2], ui=u[i], uP=u[i+2];
    for(int j=j1+j0; j <= j2; j += jstep) {
      ui[j]=factor*(hx2inv*(uM[j]+uP[j])+hy2inv*(ui[j-2]+ui[j+2])
		    -4.0*fi[j]);
    }
  }
}

template<class T>
void Pressure<T>::Smooth(const Array2<T>& u, const Array2<T>& f)
{
//  Jacobi(u,f,-1.0/(4.0*(hx2inv+hy2inv)));
  Lexicographical(u,f);
//  RedBlack(u,f);
}

template<class T>
void Pressure<T>::Restrict(const Array2<T>& r, const Array2<T>& u)
{
  if(&r != &u) {
    XDirichlet2(r,u,1);
    YDirichlet2(r,u,1);
  }
  for(int i=i1; i <= i2p; i++) {
    int i2=rx*i+offx;
    Array1(T) ri=r[i];
    Array1(T) um=u[i2-2]+offy, uz=u[i2]+offy, up=u[i2+2]+offy;
    for(int j=j1; j <= j2p; j++) {
      ri[j]=0.25*(0.5*(0.5*(um[ry*j-2]+um[ry*j+2]+up[ry*j-2]+
			    up[ry*j+2])+
		       um[ry*j]+uz[ry*j-2]+uz[ry*j+2]+up[ry*j])+
		  uz[ry*j]);
    }
  }
}

void OS::InitialConditions()
{
  Array2<Real> f0;
	
  if(convrate || verbose > 2) defectstat=1;
  if(convrate && niterations == 1) {niterations=100;}
	
  rmin=rmax=0.0;
  imin=imax=0;
	
  nlevel=max(xlevel,ylevel);
	
#if PERIODIC
  xlimP0=xlimphi=xlimpsi=xlimu=Limits(xmin,xmax,Periodic,nlevel-xlevel);
  ylimP0=ylimphi=ylimpsi=ylimu=Limits(ymin,ymax,Periodic,nlevel-ylevel);
	
  xlimP=Limits(xmin,xmax,Periodic2,nlevel-xlevel);
  ylimP=Limits(ymin,ymax,Periodic2,nlevel-ylevel);
#else
  xlimu=Limits(xmin,xmax,Neumann,nlevel-xlevel);
  ylimu=Limits(ymin,ymax,ExtendedDirichlet,nlevel-ylevel);
	
  xlimphi=xlimP0=Limits(xmin,xmax,ExtendedDirichlet,nlevel-xlevel);
  ylimphi=ylimP0=Limits(ymin,ymax,Neumann,nlevel-ylevel);
	
  xlimpsi=Limits(xmin,xmax,Neumann,nlevel-xlevel);
  ylimpsi=Limits(ymin,ymax,ExtendedDirichlet,nlevel-ylevel);
	
  xlimP=Limits(xmin,xmax,ExtendedDirichlet2,nlevel-xlevel);
  ylimP=Limits(ymin,ymax,Neumann2,nlevel-ylevel);
#endif
	
  MultiGrid< Phi<Real> > MGphi(nlevel);
  MultiGrid< Psi<Real> > MGpsi(nlevel);
  MultiGrid< Crank<Var> > MGCrank(nlevel);
#if ZERODIV	
  MultiGrid< Pressure<Real> > MGpressure(nlevel);
#else	
  MultiGrid< Pressure0<Real> > MGpressure(nlevel);
#endif	
	
  Grid2<Real> *gridphi=&MGphi.Fine();
  Grid2<Real> *gridpsi=&MGpsi.Fine();
  gridCrank=&MGCrank.Fine();
  gridP=&MGpressure.Fine();
	
  Nx=gridCrank->Nxbc();
  Ny=gridCrank->Nybc();
  
  Nxi=Nx-2;
  Nyi=Ny-2;

  int NxP=gridP->Nxbc();
  int NyP=gridP->Nybc();
	
  int oxP=gridP->Ox();
  int oyP=gridP->Oy();
	
  cout << endl << "GEOMETRY: (" << Nxi << " X " << Nyi << ")" << endl; 
	
  hxinv=1.0/gridCrank->Hx();
  hyinv=1.0/gridCrank->Hy();
  volume=gridCrank->Hx()*gridCrank->Hy();
	
  nu.vx=nuv;
  nu.vy=nuv;
  nu.C=nuC;
	
  ny=Nx*Ny;
  y=new Var[ny];
	
  u.Dimension(Nx,Ny);
  S.Dimension(Nx,Ny);
  u.Set(y);
		
  f.Allocate(Nx,Ny);
  f0.Allocate(Nx,Ny);
  phi.Allocate(Nx,Ny);
  psi.Allocate(Nx,Ny);
	
  Ex.Allocate(Nx,Ny);
  Ey.Allocate(Nx,Ny);
  Fx.Allocate(Nx,Ny);
  Fy.Allocate(Nx,Ny);
  Gx.Allocate(Nx,Ny);
  Gy.Allocate(Nx,Ny);
  Hx.Allocate(Nx,Ny);
  Hy.Allocate(Nx,Ny);
  
  P.Allocate(NxP,NyP,oxP,oyP);
  Div.Allocate(NxP,NyP,oxP,oyP);
	
  // Initialize arrays with zero boundary conditions
  Ex=0.0;
  Ey=0.0;
  Fx=0.0;
  Fy=0.0;
  Gx=0.0;
  Gy=0.0;
  Hx=0.0;
  Hy=0.0;
  u=0.0;
	
  for(int i=0; i < Nx; i++) {
    for(int j=0; j < Ny; j++) {
      // Incompressible initial velocity
#if PERIODIC
      Real x=gridCrank->X(i);
      Real y=gridCrank->Y(j);
      Real C=cos(twopi*x)*cos(twopi*x)*cos(twopi*y)*cos(twopi*y);
      Real vx=-cos(twopi*x)*sin(twopi*y);
      Real vy=sin(twopi*x)*cos(twopi*y);
#else
      // Real x=gridCrank->X(i);
      // Real y=gridCrank->Y(j);
      Real vx=0.0; // 0.5*(1.0-y*y);	
      Real vy=0.0;
      Real C=(6 <= i & i <= 10) ? 1.0 : 0.0;
#endif
      u(i,j)=U(ICvx*vx,ICvy*vy,ICC*C);
    }
  }
	
  BoundaryConditions(u);
  
#if LAGRANGIAN
  parcel.Allocate((Nx+1)*(Ny+1));
  grid.Allocate(Nx,Ny);
  u0.Dimension(Nx,Ny);
  uL.Allocate(Nx,Ny);
#endif	

  P=0.0;
  phi=0.0;
  psi=0.0;
  
  gridP->BoundaryConditions(P);	
  
  f=0.0; // Source for Poisson and Helmholtz equations
  first=1;
	
  gridphi->BoundaryConditions(phi);
  if(verbose > 3) {
    gridphi->ReportHeader();
    defect0=gridphi->DefectNorm(phi,f0);
    gridphi->Report(defect0,0);
  }
	
  for(int it=1; it < nfirst; it++) {
    gridphi->Solve(phi,f0,2,2,0);
    if(verbose > 3) { 
      defect=gridphi->DefectNorm(phi,f0);
      gridphi->Report(defect,defect0,it);
      defect0=defect;
    }
  }
	
  Real kappa=sqrt(kappa2);
  if(exactpsi) {
    for(int j=1; j <= Nyi; j++) {
      Real y=gridCrank->Y(j);
      Real solution;
      if(kappa > 1.9*log(REAL_MAX))
	solution=-exp(kappa*(y-1))-exp(-kappa*y);
      else solution=-cosh(kappa*(y-0.5))/cosh(0.5*kappa);
      for(int i=1; i <= Nxi; i++) psi(i,j)=solution;
      gridpsi->BoundaryConditions(psi);
    }
  } else {
    gridpsi->BoundaryConditions(psi);
    if(verbose > 3) {
      gridpsi->ReportHeader();
      defect0=gridpsi->DefectNorm(psi,f0);
      gridpsi->Report(defect0,0);
    }
	
    for(int it=1; it < nfirst; it++) {
      gridpsi->Solve(psi,f0,2,2,0);
      if(verbose > 3) {
	defect=gridpsi->DefectNorm(psi,f0);
	gridpsi->Report(defect,defect0,it);
	defect0=defect;
      }
    }
  }


// Electric field E = -mu*grad(phi)

  Real coeffx=0.5*hxinv;
  Real coeffy=0.5*hyinv;
  
  for(int i=1; i <= Nxi; i++) {
    for(int j=1; j <= Nyi; j++) {
      Real phix=coeffx*(phi(i+1,j)-phi(i-1,j));
      Real phiy=coeffy*(phi(i,j+1)-phi(i,j-1));
      Ex(i,j)=-mu*phix;
      Ey(i,j)=-mu*phiy;
      Fx(i,j)=alpha*psi(i,j)*phix;
      Fy(i,j)=alpha*psi(i,j)*phiy;
    }
  }	
  
  // Neumann boundary conditions on E
  
  BoundaryConditions(Ex);
  BoundaryConditions(Ey);
  
  open_output(fc,dirsep,"c"); 	 // Distribution of C
  open_output(fvel,dirsep,"Vx"); // Fluid velocity
  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");
  open_output(fdefect,dirsep,"defect");
  open_output(fudefect,dirsep,"udefect");
  if(output) open_output(fu,dirsep,"u");
  if(movie) {
    open_output(fvx,dirsep,"vx");
    open_output(fvy,dirsep,"vy");
    open_output(fC,dirsep,"C");
    open_output(fP,dirsep,"P");
    open_output(fphi,dirsep,"phi");
    open_output(fpsi,dirsep,"psi");
  }
}

void OS::Initialize()
{
  fevt << "#   t\t\t E\t\t\t Z\t\t\t I\t\t\t C" << endl;
  fdefect << "t\t\t rmin\t\t rmax\t\t imin\t\t imax" << endl;
  fudefect << "t\t\t rmin\t\t rmax\t\t imin\t\t imax" << endl;
}

void OS::Setup()
{  
#if LAGRANGIAN
// This loop supplies the size and position to each cell and initializes
// the parcels
  int k=0;
  Real hx=gridCrank->Hx();
  Real hy=gridCrank->Hy();
  
  for(int i=0; i < Nx; i++) {
    for(int j=0; j < Ny; j++) {
      grid(i,j).CellSize(hx,hy);
      grid(i,j).Position(i,j);
      parcel[k].Parameters(i*hx,j*hy,hxinv,hyinv,u(i,j));
      k++;
    }
  }
  
  int size=parcel.Size();
  for(int k=0; k < size; k++) parcel[k].Locate(u);
  	
  // Specify the outward normal vector on each boundary.
  for(int j=0; j < Ny; j++) {
    grid[0][j].Boundary(-1.0,0.0);
    grid[Nx-1][j].Boundary(1.0,0.0);
  }
  for(int i=0; i < Nx; i++) {
    grid[i][0].Boundary(0.0,-1.0);
    grid[i][Ny-1].Boundary(0.0,1.0);
  }
#endif
}

void OS::Output(int it)
{
  Real E,Z,I,C;
  u.Set(y);
  BoundaryConditions(u);
	
  ComputeInvariants(E,Z,I,C);
  fevt << t << "\t" << E << "\t" << Z << "\t" 
       << I << "\t" << C << endl;

  // Compute the distribution of C at different locations
//  int x0=(int)(0.025*Nxi);
  int x0=15;
  int x1=(int)(0.25*Nxi);
  int x3=(int)(0.75*Nxi);  
  int y1=(int)(0.5*Nyi);


  Real C0=0.0;
  Real C3=0.0;
  Real strip=1.0/(5.0*41.0);

  for(int i=x0-2; i <= x0+2; i++){
    for(int j=y1-20; j <= y1+20; j++){
      C0 += u[i][j].C;
    }
  }

  C0*=strip;

  for(int i=x3-2; i <= x3+2; i++){
    for(int j=y1-20; j <= y1+20; j++){
      C3 += u[i][j].C;
    }
  }

  C3*=strip;

//  fc << t << "\t" << C0 << "\t" << C3 << endl;
  fc << t << "\t" << C3 << endl;

// Compute the electro-osmotic velocity profile
  for(int j=1; j <= Nyi; j++) {
    Real y=gridCrank->Y(j);
    if(it == 30 || it == 60 || it == 90)
      fvel << y << "\t" << u(x1,j).vx << "\t" << u(x3,j).vx << endl;
  }



  if(defectstat) {
    fdefect << t << "\t" << rmin << "\t"  << rmax << "\t" << imin << "\t" 
	    << imax << endl;
    fudefect << t << "\t" << urmin << "\t" << urmax << "\t" << uimin
	     << "\t" << uimax << endl;
  }
	
  rmin=REAL_MAX; rmax=0.0;
  imin=INT_MAX; imax=0;
	
  urmin=REAL_MAX; urmax=0.0;
  uimin=INT_MAX; uimax=0;
	
  if(output) {
    out_curve(fu,u(),"u",ny);
    fu.flush();
  }
  if(movie) OutFrame(it);
	
  ft << t << endl;
}

void OS::OutFrame(int it)
{
  fvx << Nxi << Nyi << 1;
  fvy << Nxi << Nyi << 1;
  fC << Nxi << Nyi << 1;
  fP << Nxi << Nyi << 1;
  if(it == 0) {
    fphi << Nxi << Nyi << 1;
    fpsi << Nxi << Nyi << 1;
  }
  
  for(int j=Nyi; j >= 1; j--) {
    for(int i=1; i <= Nxi; i++) {
      fvx << (float) u(i,j).vx;
      fvy << (float) u(i,j).vy;
      fC << (float) u(i,j).C;
      fP << (float) P(i,j);
      if(it == 0) {
	fphi << (float) phi(i,j);
	fpsi << (float) psi(i,j);
      }
    }
  }
  
  if(it == 0) {
    fphi.close();
    fpsi.close();
  }
  
  fvx.flush();
  fvy.flush();
  fC.flush();
  fP.flush();
}	

void OS::ComputeInvariants(Real& E, Real& Z, Real& I, Real& C)
{
  E=Z=I=C=0.0;
	
  Real coeffx=0.5*hxinv;
  Real coeffy=0.5*hyinv;
  for(int i=1; i <= Nxi; i++) {
    Array1(U) um=u[i-1], ui=u[i], up=u[i+1];
    for(int j=1; j <= Nyi; j++) {
      E += ui[j].vx*ui[j].vx+ui[j].vy*ui[j].vy;
      Real w=-coeffx*(um[j].vy-up[j].vy)+coeffy*(ui[j-1].vx-ui[j+1].vx);
      Z += w*w;
      I += ui[j].C*w;
      C += ui[j].C*ui[j].C;
    }
  }
  
  Real factor=0.5*volume;
  E *= factor;
  Z *= factor;
  I *= factor;
  C *= factor;
}

void OS::Transform(Var *Y, double, double dt, Var *&yi)
{
  if(!yi) yi=new Var[ny];
  set(yi,Y,ny);
	
  beta=explicit_factor*dt*nu;
	
  u.Set(Y);
  gridCrank->Lu(u,f);
  u=f;
}

void OS::BackTransform(Var *Y, double, double dt, Var *yi)
{
  beta=-implicit_factor*dt*nu;
	
  u.Set(Y);
  f=y;
  u.Load(yi);
	
  int it=gridCrank->Solve(u,f,npresmooth,igamma,npostsmooth,0,
			  ncrank,convrate,
			  defectstat,&udefect0,&udefect);

  if(verbose > 2) {
    cout << endl;
    gridCrank->Report(udefect0,0);
    gridCrank->Report(udefect,udefect0,it);
  }
	
  if(!it) 
    msg(ERROR,"Maximum number of iterations (%d) exceeeded in Solve",
	niterations);
	
  if(defectstat) {
    U ratio=divide0(udefect0,udefect);
    if(ratio < urmin) urmin=ratio;
    if(ratio > urmax) urmax=ratio;
  }
	
  BoundaryConditions(u);
    
#if LAGRANGIAN
  
  u0.Set(yi);
  
  for(int i=0; i < Nx; i++) {
    Array1(Cell) gridi=grid[i];
    Array1(U) uLi=uL[i], ui=u[i], u0i=u0[i];
    Array1(Real) Exi=Ex[i], Eyi=Ey[i];
    for(int j=0; j < Ny; j++) {
      uLi[j].vx=0.5*(ui[j].vx+u0i[j].vx)+Exi[j];
      uLi[j].vy=0.5*(ui[j].vy+u0i[j].vy)+Eyi[j];
      uLi[j].C=ui[j].C;
      grid(i,j).Input(uLi[j]);
      gridi[j].Source(ui[j].C-u0i[j].C); // Parcel diffusion
    }
  }
  
  int size=parcel.Size();
  for(int k=0; k < size; k++) parcel[k].Source(grid);
      
  for(int i=0; i < Nx; i++) {
    Array1(Cell) gridi=grid[i];
    Array1(U) uLi=uL[i];
    for(int j=0; j < Ny; j++) {
      if(gridi[j].Boundary()) InjectParcel(gridi[j],parcel,Inactive);
      uLi[j].C=0.0;
      gridi[j].Area(0.0);
    }
  }
      
  size=parcel.Size();
  for(int k=0; k < size; k++) {
    parcel[k].Advect(dt);
    if(parcel[k].Locate(uL)) Inactive.Push(k);
    else parcel[k].Contribute(uL,grid);
  }
		
  for(int i=0; i < Nx; i++) {
    Array1(U) uLi=uL[i], ui=u[i];
    for(int j=0; j < Ny; j++) {
      ui[j].C=uLi[j].C;
    }
  }
  
  BoundaryConditions(u);
#endif

  if(it < uimin) uimin=it;
  if(it > uimax) uimax=it;
}

void OS::Source(Var *source, Var *Y, double)
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

#if 0
  Div=0.0;
  for(int i=1; i <= Nxi; i++) {
    Array1(U) um=u[i-1], up=u[i+1], ui=u[i];
    for(int j=1; j <= Nyi; j++) {
      // Div = div u
      Div(i,j)=-coeffx*(um[j].vx-up[j].vx)-
	coeffy*(ui[j-1].vy-ui[j+1].vy);
    }
  }

  cout << Div << endl;
#endif	

  Real coeffx0=coeffx*nonlinear;
  Real coeffy0=coeffy*nonlinear;
  
  for(int i=1; i <= Nxi; i++) {
    Array1(U) Gxm=Gx[i-1], Gxp=Gx[i+1];
#if UPWIND
    Array1(U) Gxi=Gx[i];
#endif    
    Array1(U) Gyi=Gy[i];
    
    Array1(Real) Hxi=Hx[i], Hyi=Hy[i];
    Array1(Real) Fxi=Fx[i], Fyi=Fy[i];
    Array1(U) Si=S[i];

// center in space for momentum equation; optionally upwind C equation

    for(int j=1; j <= Nyi; j++) {
      // div (v u)
      U advection=-coeffx0*(Gxm[j]-Gxp[j])-coeffy0*(Gyi[j-1]-Gyi[j+1]);
      
#if UPWIND
      // div(v C + mu E C)
      advection.C=
	hxinv*(u(i,j).vx > 0 ? (Gxi[j].C-Gxm[j].C) : (Gxp[j].C-Gxi[j].C))+
	hyinv*(u(i,j).vy > 0 ? (Gyi[j].C-Gyi[j-1].C) : (Gyi[j+1].C-Gyi[j].C));
#endif      
      
      // H=F - div (v u)
      Hxi[j]=Fxi[j]-advection.vx;
      Hyi[j]=Fyi[j]-advection.vy;
      
#if LAGRANGIAN
      Si[j].C=0.0;
#else
      Si[j].C=-advection.C;
#endif      
    }
  }

  BoundaryConditions(Hx);
  BoundaryConditions(Hy);

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

    int it=gridP->Solve(P,Div,npresmooth,igamma,npostsmooth,0,
			first ? max(nfirst,niterations) : niterations,
			convrate,defectstat,&defect0,&defect);
	
    first=0;
	
    if(verbose > 2) {
      cout << endl;
      gridP->Report(defect0,0);
      gridP->Report(defect,defect0,it);
    }
		
    if(!it) msg(ERROR,"Maximum number of iterations (%d) exceeeded in Solve",
		niterations);
		
    if(defectstat) {
      U ratio=divide0(defect0,defect);
      if(ratio < rmin) rmin=ratio;
      if(ratio > rmax) rmax=ratio;
    }
	
    if(it < imin) imin=it;
    if(it > imax) imax=it;
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
