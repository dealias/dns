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
Real nu=1.0;
Real rho=1.0;
Real P0=1.0;
Real P1=1.0;
Real ICvx=1.0;
Real ICvy=1.0;
Real xmin=0.0;
Real xmax=1.0;
Real ymin=0.0;
Real ymax=1.0;
Real force=1.0;
Real kforce=1.0;
Real deltaf=2.0;

int movie=0;

int dumpbins=0;

// Global Variables;
//unsigned int Nmode;
extern unsigned Nxb,Nxb1,Nyb,Nyp;
Real shellmin2;
Real shellmax2;

// Local Variables;
Array1<Real> kxmask;
Array1<Real> kymask;
Array1<Real> k2mask, k2invmask;

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
  int iNxb,iNyb;
  Real hxinv, hyinv;
  unsigned nmode;
	
  Array3<Real> u,S,f;
  unsigned nspecies;
  
  Array2<Real> us,dudx,dudy;
  Array2<Complex> uk,ikxu,ikyu;
  
  Array3<Real> Ss;
  Array3<Complex> Sk;
  
  Array3<Real> ForceMask;
  Array3<Complex> FMk;
  
  Real coeffx,coeffy;
  
public:
  DNS();
  virtual ~DNS() {}
  void InitialConditions();
  void Initialize();
  void Output(int it);
  void OutFrame(int it);
  void Source(Var *, Var *, double);
  void ComputeInvariants(Real& E, Real& Z);
  void Stochastic(Var *y, double, double dt);
  void Spectrum(Array3<Complex>& Sk, ofstream& s, const char *text);
  
  Real X(int i) {return xmin+i*(xmax-xmin)/Nxb;}
  Real Y(int i) {return ymin+i*(ymax-ymin)/Nyb;}
  
  Real Vorticity(int i, int j) {
    return coeffx*(u(i+1 < iNxb ? i+1 : 0,j,1)-u(i > 0 ? i-1 : iNxb-1,j,1))-
      coeffy*(u(i,j+1 < iNyb ? j+1 : 0,0)-u(i,j > 0 ? j-1 : iNyb-1,0));
  }

};

DNS *DNSProblem;

DNS::DNS()
{
  DNSProblem=this;
}

DNSVocabulary::DNSVocabulary()
{
  Vocabulary=this;

  VOCAB(rho,0.0,REAL_MAX,"density");
  VOCAB(nu,0.0,REAL_MAX,"Kinematic viscosity");
  VOCAB(force,0.0,REAL_MAX,"force coefficient");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");
  
  VOCAB(P0,-REAL_MAX,REAL_MAX,"");
  VOCAB(P1,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICvx,-REAL_MAX,REAL_MAX,"Initial condition multiplier for vx");
  VOCAB(ICvy,-REAL_MAX,REAL_MAX,"Initial condition multiplier for vy");

  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  
  VOCAB(xmin,-REAL_MAX,REAL_MAX,"Minimum x value");
  VOCAB(xmax,-REAL_MAX,REAL_MAX,"Maximum x value");
	
  VOCAB(ymin,-REAL_MAX,REAL_MAX,"Minimum y value");
  VOCAB(ymax,-REAL_MAX,REAL_MAX,"Maximum y value");

  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  
  GeometryTable=new Table<GeometryBase>("geometry");
//  LinearityTable=new Table<LinearityBase>("linearity");
//  ForcingTable=new Table<ForcingBase>("forcing");
//  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  
  BASIS(Cartesian);
  
  METHOD(DNS);
}

DNSVocabulary DNS_Vocabulary;

ofstream ft,fevt,fekvk,fu;
oxstream fvx,fvy,fw;

void DNS::InitialConditions()
{
  nspecies=2;
  set_fft_parameters();
  
  cout << endl << "GEOMETRY: (" << Nxb << " X " << Nyb << ")" << endl; 
	
  Geometry=DNS_Vocabulary.NewGeometry(geometry);
  Geometry->Create(0);
  nmode=Geometry->nMode();
  
  iNxb=Nxb;
  iNyb=Nyb;
  coeffx=0.5*(xmax-xmin)/Nxb;
  coeffy=0.5*(ymax-ymin)/Nyb;
  
  ny=Nxb*Nyb*nspecies;
  y=new Var[ny];
  
  u.Dimension(Nxb,Nyb,nspecies);
  S.Dimension(Nxb,Nyb,nspecies);
  
  us.Allocate(Nxb,2*Nyp);
  uk.Dimension(Nxb,Nyp,(Complex *) us());
  
  dudx.Allocate(Nxb,2*Nyp);
  ikxu.Dimension(Nxb,Nyp,(Complex *) dudx());
  
  dudy.Allocate(Nxb,2*Nyp);
  ikyu.Dimension(Nxb,Nyp,(Complex *) dudy());
  
  Ss.Allocate(nspecies,Nxb,2*Nyp);
  Sk.Dimension(nspecies,Nxb,Nyp,(Complex *) Ss());
  
  ForceMask.Allocate(nspecies,Nxb,2*Nyp);
  FMk.Dimension(nspecies,Nxb,Nyp,(Complex *) ForceMask());
  
  u.Set(y);
		
  // Initialize arrays with zero boundary conditions
  u=0.0;
	
  for(unsigned i=0; i < Nxb; i++) {
    for(unsigned j=0; j < Nyb; j++) {
#if 0     
      Real x=X(i);
      Real y=Y(j);
      Real vx=-cos(twopi*x)*sin(twopi*y);
      Real vy=sin(twopi*x)*cos(twopi*y);
#else      
      Real vx=drand();
      Real vy=drand();
#endif
      u(i,j,0)=ICvx*vx;
      u(i,j,1)=ICvy*vy;
    }
  }
	
  // Filter high-components from velocity and enforce solenoidal condition
  
  for(unsigned s=0; s < nspecies; s++) {
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	Ss(s,i,j)=u(i,j,s);
      }	
    }
    rcfft2d(Sk[s],log2Nxb,log2Nyb,-1);
  }

  for(unsigned i=0; i < nfft; i++) {
    Real kx=kxmask[i];
    Real ky=kymask[i];
    if(k2invmask[i]) {
      // Calculate -i*P
      Complex miP=(kx*Sk[0](i)+ky*Sk[1](i))*k2invmask[i];
      Sk[0](i)=(Sk[0](i)-kx*miP)*Nxybinv;
      Sk[1](i)=(Sk[1](i)-ky*miP)*Nxybinv;
    } else {
      Sk[0](i)=Sk[1](i)=0.0;
    }
  }
  
  for(unsigned s=0; s < nspecies; s++) {
    Array2<Real> Sss=Ss[s];
    Array2<Complex> Sks=Sk[s];
    crfft2d(Sks,log2Nxb,log2Nyb,1);
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	u(i,j,s)=Sss(i,j);
      }
    }
  }
  
  // Apply band-limited forcing, enforcing solenoidal condition
  
  for(unsigned i=0; i < nfft; i++) {
    Real kx=kxmask[i];
    Real ky=kymask[i];
    if(k2invmask[i]) {
      // Calculate -i*P
      Real K=sqrt(kx*kx+ky*ky);
      if(K >= kforce-0.5*deltaf && K <= kforce+0.5*deltaf) {
	FMk[0](i)=force*crand_gauss();
	FMk[1](i)=force*crand_gauss();
      }
      else FMk[0](i)=FMk[1](i)=0.0;
      Complex miP=(kx*FMk[0](i)+ky*FMk[1](i))*k2invmask[i];
      FMk[0](i)=(FMk[0](i)-kx*miP);
      FMk[1](i)=(FMk[1](i)-ky*miP);
    } else {
      FMk[0](i)=FMk[1](i)=0.0;
    }
  }
  
  for(unsigned s=0; s < nspecies; s++) crfft2d(FMk[s],log2Nxb,log2Nyb,1);
  
  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");
  
  if(output) open_output(fu,dirsep,"u");
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


void Basis<Cartesian>::Initialize()
{
  cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << " x " << Nyp
       << ")." << endl;
  
  Array1<Complex> temp(nmode);
  Array1<Complex> mask(nfft);
  
  kxmask.Allocate(nfft);
  kymask.Allocate(nfft);
  k2mask.Allocate(nfft);
  k2invmask.Allocate(nfft);
  
  for(unsigned int i=0; i < nmode; i++) temp[i]=I*CartesianMode[i].X();
  CartesianPad(mask,temp);
  for(unsigned int i=0; i < nfft; i++) kxmask[i]=mask[i].im;
  
  for(unsigned int i=0; i < nmode; i++) temp[i]=I*CartesianMode[i].Y();
  CartesianPad(mask,temp);
  for(unsigned int i=0; i < nfft; i++) {
    kymask[i]=mask[i].im;
    Real k2=kxmask[i]*kxmask[i]+kymask[i]*kymask[i];
    k2mask[i]=k2;
    k2invmask[i]=k2 ? 1.0/k2 : 0.0;
  }
}

void DNS::Spectrum(Array3<Complex>& Sk, ofstream& os, const char *text)
{
  unsigned Nxb2=Nxb/2, Nyb2=Nyb/2;
  unsigned Kbmax=(unsigned) (sqrt(Nxb2*Nxb2+Nyb2*Nyb2)+0.5);
  
  unsigned Nx2=(Nx-1)/2, Ny2=(Ny-1)/2;
  unsigned Kmax=(unsigned) (sqrt(Nx2*Nx2+Ny2*Ny2)+0.5);
  
  Real *sum=new Real[Kbmax+1];
  unsigned *count=new unsigned[Kbmax+1];
		
  for(unsigned K=0; K <= Kbmax; K++) {
    count[K]=0;
    sum[K]=0.0;
  }
				
  // Compute instantaneous angular average over circular shell.
		
  for(unsigned s=0; s < nspecies; s++) {
    for(unsigned i=0; i < Nxb; i++) {
      int kx=i-Nxb2;
      for(unsigned j=0; j < Nyb; j++) {
	int jp=j-Nyb2;
	int K2=kx*kx+jp*jp;
	int K=(int)(sqrt(K2)+0.5);
	int ip=i;
	if(jp < 0) {
	  if(i > 0) ip=Nxb-i;
	  if(j > 0) jp=-jp;
	  else jp=Nyb2;
	}

	count[K]++;
	sum[K] += abs2(Sk(s,ip,jp));
      }
    }
  }
		
  for(unsigned K=0; K <= Kmax; K++)
    if(count[K]) sum[K] *= 0.5/count[K] * 2.0*pi*K;
				
  open_output(os,dirsep,text,0);
  out_curve(os,t,"t");
  out_curve(os,sum,text,Kmax+1);
  os.close();
}

void DNS::Output(int it)
{
  Real E,Z;
	
  u.Set(y);
  ComputeInvariants(E,Z);
  fevt << t << "\t" << E << "\t" << Z << endl;

  if(movie) OutFrame(it);
	
  for(unsigned s=0; s < nspecies; s++) {
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	Ss(s,i,j)=u(i,j,s);
      }
    }
    rcfft2d(Sk[s],log2Nxb,log2Nyb,-1);
  }
  
  Spectrum(Sk,fekvk,"ekvk");
    
  ft << t << endl;
}

void DNS::OutFrame(int it)
{
  fvx << Nxb << Nyb << 1;
  fvy << Nxb << Nyb << 1;
  fw << Nxb << Nyb << 1;  
  
  for(int j=iNyb-1; j >= 0; j--) {
    for(int i=0; i < iNxb; i++) {
      fvx << (float) u(i,j,0);
      fvy << (float) u(i,j,1);
      fw << (float) Vorticity(i,j);
    }
  }
  
  fvx.flush();
  fvy.flush();
  fw.flush();
}	

void DNS::ComputeInvariants(Real& E, Real& Z)
{
  E=Z=0.0;
	
  for(unsigned i=0; i < Nxb; i++) {
    Array2<Real> ui=u[i];
    for(unsigned j=0; j < Nyb; j++) {
    Array1(Real) uij=ui[j];
      for(unsigned s=0; s < nspecies; s++) {
	E += uij[s]*uij[s];
      }
      
      Real w=Vorticity(i,j);
      Z += w*w;
    }
  }
  
  Real volume=(xmax-xmin)*(ymax-ymin)/(Nxb*Nyb);
  Real factor=0.5*volume;
  
  E *= factor;
  Z *= factor;
}

void DNS::Source(Var *source, Var *y, double t)
{
  u.Set(y);
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
	Ss(s,i,j)=-(u(i,j,0)*dudx(i,j)+u(i,j,1)*dudy(i,j))*Nxybinv;
      }
    }
    
    rcfft2d(Sk[s],log2Nxb,log2Nyb,-1);
    
    for(unsigned i=0; i < nfft; i++) {
      if(k2mask[i])
	Sk[s](i)=(Sk[s](i)-nu*k2mask[i]*uk(i))*Nxybinv;
      else
	Sk[s](i)=0.0;
    }
  }
  
  for(unsigned i=0; i < nfft; i++) {
    Real kx=kxmask[i];
    Real ky=kymask[i];
    // Calculate -i*P
    Complex miP=(kx*Sk[0](i)+ky*Sk[1](i))*k2invmask[i];
    Sk[0](i) -= kx*miP;
    Sk[1](i) -= ky*miP;
  }
  
  for(unsigned s=0; s < nspecies; s++) {
    Array2<Real> Sss=Ss[s];
    crfft2d(Sk[s],log2Nxb,log2Nyb,1);
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	S(i,j,s)=Sss(i,j);
      }
    }
  }
}

void DNS::Stochastic(Var *y, double, double dt)
{
  u.Set(y);
  Real factor=sqrt(2.0*dt)*rand_gauss();
  
  for(unsigned s=0; s < nspecies; s++) {
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	u(i,j,s) += factor*ForceMask(s,i,j);
      }
    }
  }
}
