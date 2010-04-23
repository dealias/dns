#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "fftw++.h"
#include "Forcing.h"
#include "InitialCondition.h"
#include "convolution.h"

using namespace Array;

const double ProblemVersion=1.0;

const char *method="DNS";
const char *integrator="RK5";
const char *geometry="Cartesian";
const char *forcing="WhiteNoiseBanded";
const char *ic="Equilibrium";

ForcingBase *Forcing;
InitialConditionBase *InitialCondition;

// Vocabulary
Real nu=1.0;
Real rho=1.0;
Real ICvx=1.0;
Real ICvy=1.0;
unsigned int Nx=1;
unsigned int Ny=1;
Real xmin=0.0;
Real xmax=1.0;
Real ymin=0.0;
Real ymax=1.0;
Real force=1.0;
Real kforce=1.0;
Real deltaf=2.0;
int movie=0;
int rezero=0;

typedef array1<Complex>::opt cvector;
typedef array2<Complex>::opt cvector2;

//enum Field {OMEGA,EK};
enum Field {OMEGA};

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
  unsigned int mx, my; // size of data arrays
  Real kx0, ky0; // grid spacing factor
  unsigned Nx2, Ny2; // half of Nx, Ny
  
  vector w; // array pointer for vorticity

  
  unsigned Nx1;
  unsigned Nyp;
  unsigned nfft;
  
  Real Nxybinv;
  
  int tcount;
  
  array1<unsigned>::opt count;
  
  Real hxinv, hyinv;
  unsigned nmode;
	
  array3<Complex> u,S;
  
  unsigned nshells;  // Number of spectral shells
  
  rvector spectrum;
  
  array3<Real> ForceMask;
  array3<Complex> FMk;
  
  Real coeffx,coeffy;
  array3<Complex> f,g;
  array2<Complex> f0,f1,g0,g1;
  
public:
  DNS();
  virtual ~DNS() {}
  void InitialConditions();
  ImplicitHConvolution2 Convolution;
  void Initialize();
  void Output(int it);
  void FinalOutput();
  void OutFrame(int it);
  void Source(const vector2& Src, const vector2& Y, double t);
  void ComputeInvariants(Real& E, Real& Z);
  void Stochastic(const vector2& Y, double, double);
  void SetFFTparameters();
  
  void Spectrum(vector& S, const vector& y);
  
  Real Xof(int i) {return xmin+i*(xmax-xmin)/Nx;}
  Real Yof(int i) {return ymin+i*(ymax-ymin)/Ny;}

  /*  
  Real Vorticity(unsigned int i, unsigned int j) {
    return coeffx*(u(i+1 < Nx ? i+1 : 0,j,1)-u(i > 0 ? i-1 : Nx-1,j,1))-
      coeffy*(u(i,j+1 < Ny ? j+1 : 0,0)-u(i,j > 0 ? j-1 : Ny-1,0));
  }
  */

};

DNS *DNSProblem;

DNS::DNS()
{
  DNSProblem=this;
  check_compatibility(DEBUG);
}

DNSVocabulary::DNSVocabulary()
{
  Vocabulary=this;

  VOCAB(rho,0.0,REAL_MAX,"density");
  VOCAB(nu,0.0,REAL_MAX,"Kinematic viscosity");
  VOCAB(force,0.0,REAL_MAX,"force coefficient");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");
  
  VOCAB(ICvx,0.0,0.0,"Initial condition multiplier for vx");
  VOCAB(ICvy,0.0,0.0,"Initial condition multiplier for vy");

  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  
  VOCAB(xmin,0.0,0.0,"Minimum x value");
  VOCAB(xmax,0.0,0.0,"Maximum x value");
	
  VOCAB(ymin,0.0,0.0,"Minimum y value");
  VOCAB(ymax,0.0,0.0,"Maximum y value");

  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  
  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");
  
//  ForcingTable=new Table<ForcingBase>("forcing");
//  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  
  METHOD(DNS);
}

DNSVocabulary DNS_Vocabulary;

ifstream ftin;
ofstream ft,fevt,fu;
oxstream fvx,fvy,fw,fekvk;

void DNS::InitialConditions()
{
  kx0=1.0;
  ky0=1.0;
  mx=Nx;
  my=Ny/2;  // FIXME: plus one?

  Nx2=Nx/2, Ny2=Ny/2;
  
  unsigned Nx2=(Nx-1)/2, Ny2=(Ny-1)/2;
  nshells=(unsigned) (sqrt(Nx2*Nx2+Ny2*Ny2)+0.5);
  
  NY[OMEGA]=Nx*Ny;
  //  NY[EK]=nshells;
  
  f0.Dimension(mx,my);
  f1.Allocate(mx,my);
  g0.Allocate(mx,my);
  g1.Allocate(mx,my);
  
  Allocator();
  
  unsigned int align=sizeof(Complex);  
  f.Allocate(2,2*Nx,my,align);
  g.Allocate(2,2*Nx,my,align);
  
  //u.Dimension(Nx,Ny);
  ForceMask.Allocate(Nx,2*Nyp,align);
  //FMk.Dimension(Nx,Nyp,(Complex *) ForceMask());
  
  cout << endl << "GEOMETRY: (" << Nx << " X " << Ny << ")" << endl; 

  //Convolution.init(Nx,my,2);
  //Convolution=new ImplicitHConvolution2(Nx,my,2);

  Allocate(count,nshells);
  
  //u.Set(Y[OMEGA]);
		
  // Initialize arrays with zero boundary conditions
  u=0.0;
	
  //  for(unsigned i=0; i < NY[EK]; i++) Y[EK][i]=0.0;
  
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
  //mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
  
  errno=0;
    
  if(output) 
    open_output(fu,dirsep,"u");
  if(movie) {
    open_output(fvx,dirsep,"vx");
    open_output(fvy,dirsep,"vy");
    open_output(fw,dirsep,"w");
  }
}

void DNS::Initialize()
{
  fevt << "#   t\t\t E\t\t\t Z" << endl;
}

/*
void DNS::Spectrum(vector& S, const vector& y)
{
  u.Set(y);
  
  for(unsigned K=0; K < nshells; K++) {
    count[K]=0;
    S[K]=0.0;
  }
				
  // Compute instantaneous angular average over circular shell.
		
  for(unsigned i=0; i < Nx; i++) {
    int kx=i-Nx2;
    rvector k2maski=k2maskij[i];
    for(unsigned j=(kx > 0) ? 0 : 1; j < Nyp; j++) {
      if(k2maski[j]) {
        int K2=kx*kx+j*j;
        unsigned K=(unsigned)(sqrt(K2)-0.5);
        count[K]++;
        S[K] += abs2(Sk(i,j));
      }
    }
  }
  
  for(unsigned K=0; K < nshells; K++)
    if(count[K]) S[K] *= twopi*K/count[K];
}
*/

void DNS::Output(int it)
{
  /*
  Real E,Z;
	
  u.Set(y);
  ComputeInvariants(E,Z);
  fevt << t << "\t" << E << "\t" << Z << endl;

  Var *y=Y[0];
  if(output) out_curve(fu,y,"u",NY[0]);
  
  if(movie) OutFrame(it);
	
  ostringstream buf;
  buf << "ekvk" << dirsep << "t" << tcount;
  open_output(fekvk,dirsep,buf.str().c_str(),0);
  out_curve(fekvk,t,"t");
  Var *y1=Y[EK];
  out_curve(fekvk,y1,"ekvk",nshells);
  fekvk.close();
  if(!fekvk) msg(ERROR,"Cannot write to file ekvk");
    
  tcount++;
  ft << t << endl;
  
  if(rezero && it % rezero == 0) {
    vector2 Y=Integrator->YVector();
    vector T=Y[EK];
    for(unsigned int i=0; i < NY[EK]; i++)
      T[i]=0.0;
  }
  */
}

void DNS::OutFrame(int)
{
  /*
  fvx << 1 << Ny << Nx;
  fvy << 1 << Ny << Nx;
  fw << 1 << Ny << Nx;  
  
  for(int j=iNy-1; j >= 0; j--) {
    for(unsigned int i=0; i < Nx; i++) {
      fvx << (float) u(i,j,0);
      fvy << (float) u(i,j,1);
      fw << (float) Vorticity(i,j);
    }
  }
  
  fvx.flush();
  fvy.flush();
  fw.flush();
  */
}	

void DNS::ComputeInvariants(Real& E, Real& Z)
{
  /*
  E=Z=0.0;
	
  for(unsigned i=0; i < Nx; i++) {
    array1<Complex>::opt ui=u[i];
    for(unsigned j=0; j < Ny; j++) {
      Complex uij=ui[j];
      E += uij*uij;
      Real w=Vorticity(i,j);
      Z += w*w;
    }
  }
  
  E *= factor;
  Z *= factor;
  */
}

void DNS::FinalOutput()
{
  /*
  Real E,Z;
  ComputeInvariants(E,Z);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  */
}

void DNS::Source(const vector2& Src, const vector2& Y, double)
{
  /*
  //w.Set(Y[OMEGA]);
  //S.Set(Src[OMEGA]);
  
  f0.Set(Src[OMEGA]); // f0 is complex, but Src[OMEGA] is real. FIXME
  
  unsigned int xorigin=mx-1;
  f0[xorigin][0]=0.0; // Move out later
  f1[xorigin][0]=0.0;
  g0[xorigin][0]=0.0;
  g1[xorigin][0]=0.0;
  for(unsigned int i=0; i < Nx; ++i) {
    Real kx=kx0*(i-xorigin);
    Real kx2=kx*kx;
    vector wi=w[i];
    for(unsigned int j=i == xorigin ? 1 : 0; j < my; ++j) {
      Real ky=ky0*j;
      Complex w0=wi[j];
      Complex kxw=kx*w0;
      Complex kyw=ky*w0;
      f0[i][j]=kxw;
      f1[i][j]=-kyw;
      Real k2inv=1.0/(kx2+ky*ky);
      g0[i][j]=k2inv*kyw;
      g1[i][j]=k2inv*kxw;
    }
  }
  
  Complex *F[2]={f0,f1};
  Complex *G[2]={g0,g1};
  Convolution.convolve(F,G);
  
  //  Spectrum(Src[EK],Y[OMEGA]);
  */
}

void DNS::Stochastic(const vector2&Y, double, double)
{
//  u.Set(Y[OMEGA]);
//  Real factor=sqrt(2.0*dt)*rand_gauss();
}
