#include "options.h"
#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"
#include "Array.h"
#include "fftw++.h"
#include "Linearity.h"
#include "Forcing.h"
#include "InitialCondition.h"

using namespace Array;

const double ProblemVersion=1.0;

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
int nparticles=1;
int movie=0;

typedef array1<Complex>::opt cvector;

// Unused variables:

Real shellmin2;
Real shellmax2;
int dumpbins=0;
Real krmin;
Real krmax;
int reality=1;

enum Field {VEL,EK,PX,PV};

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
  unsigned Nxb2, Nyb2; // half of Nxb, Nyb
  
  unsigned Nxb,Nxb1,log2Nxb;
  unsigned Nyb,Nyp,log2Nyb;
  unsigned nfft;
  
  Real Nxybinv;
  
  rvector kxmask;
  rvector kymask;
  rvector k2mask, k2invmask;
  array2<Real> k2maskij;
  int tcount;
  
  array1<unsigned>::opt count;
  
  Real hxinv, hyinv;
  unsigned nmode;
	
  array3<Real> u,S,f;
  unsigned nspecies;
  
  unsigned nshells;  // Number of spectral shells
  
  array2<Real> us,dudx,dudy;
  array2<Complex> uk,ikxu,ikyu;
  
  array3<Real> Ss;
  array3<Complex> Sk;
  
  rvector spectrum;
  
  array3<Real> ForceMask;
  array3<Complex> FMk;
  
  Real coeffx,coeffy;
  crfft2d *cr;
  rcfft2d *rc;
  
public:
  DNS();
  virtual ~DNS() {}
  void InitialConditions();
  void Initialize();
  void Output(int it);
  void FinalOutput();
  void OutFrame(int it);
  void Source(const vector2& Src, const vector2& Y, double t);
  void ComputeInvariants(Real& E, Real& Z);
  void Stochastic(const vector2& Y, double, double);
  void SetFFTparameters();
  
  void Spectrum(vector& S, const vector& y);
  
  Real Xof(int i) {return xmin+i*(xmax-xmin)/Nxb;}
  Real Yof(int i) {return ymin+i*(ymax-ymin)/Nyb;}
  
  Real Vorticity(int i, int j) {
    return coeffx*(u(i+1 < iNxb ? i+1 : 0,j,1)-u(i > 0 ? i-1 : iNxb-1,j,1))-
      coeffy*(u(i,j+1 < iNyb ? j+1 : 0,0)-u(i,j > 0 ? j-1 : iNyb-1,0));
  }

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
  
  VOCAB(P0,0.0,0.0,"Pressure at xmin (not yet implemented)");
  VOCAB(P1,0.0,0.0,"Pressure at xmax (not yet implemented)");
  VOCAB(ICvx,0.0,0.0,"Initial condition multiplier for vx");
  VOCAB(ICvy,0.0,0.0,"Initial condition multiplier for vy");

  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  VOCAB(nparticles,1,INT_MAX,"Number of particles");
  
  VOCAB(xmin,0.0,0.0,"Minimum x value");
  VOCAB(xmax,0.0,0.0,"Maximum x value");
	
  VOCAB(ymin,0.0,0.0,"Minimum y value");
  VOCAB(ymax,0.0,0.0,"Maximum y value");

  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  
  GeometryTable=new Table<GeometryBase>("geometry");
//  LinearityTable=new Table<LinearityBase>("linearity");
//  ForcingTable=new Table<ForcingBase>("forcing");
//  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  
  BASIS(Cartesian);
  
  METHOD(DNS);
}

DNSVocabulary DNS_Vocabulary;

ifstream ftin;
ofstream ft,fevt,fu;
oxstream fvx,fvy,fw,fekvk;

void DNS::SetFFTparameters()
{
  Geometry->SetFFTparameters();
  Nxb=Geometry->getNxb();
  Nxb1=Geometry->getNxb1();
    
  Nyb=Geometry->getNyb();
  Nyp=Geometry->getNyp();
    
  log2Nxb=Geometry->getlog2Nxb();
  log2Nyb=Geometry->getlog2Nyb();
    
  Nxybinv=1.0/(Nxb*Nyb);
  nfft=Nxb1*Nyp;
  
//  array2<Complex> Sk0=Sk[0];
  cr=new crfft2d(Sk[0]);
  rc=new rcfft2d(Sk[0]);
  
  cout << endl << "GEOMETRY: (" << Nxb << " X " << Nyb << ")" << endl; 
}

void DNS::InitialConditions()
{
  nspecies=2;
  
  Geometry=DNS_Vocabulary.NewGeometry(geometry);
  Geometry->Create(0);
  nmode=Geometry->nMode();
  
  SetFFTparameters();
	
  iNxb=Nxb;
  iNyb=Nyb;
  coeffx=0.5*Nxb/(xmax-xmin);
  coeffy=0.5*Nyb/(ymax-ymin);
  
  Nxb2=Nxb/2, Nyb2=Nyb/2;
  
  unsigned Nx2=(Nx-1)/2, Ny2=(Ny-1)/2;
  unsigned Kmax=(unsigned) (sqrt(Nx2*Nx2+Ny2*Ny2)+0.5);
  nshells=Kmax+1;
  
  NY[VEL]=Nxb*Nyb*nspecies;
  NY[EK]=nshells;
//  NY[PX]=nparticles*nspecies;
//  NY[PV]=nparticles*nspecies;
  
  Allocator();
  
//  position.Dimension(nparticles,nspecies);
//  velocity.Dimension(nparticles,nspecies);
  
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
  
  cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << " x " << Nyp
       << ")." << endl;
  
  cvector temp;
  cvector mask;
  
  Allocate(kxmask,nfft);
  Allocate(kymask,nfft);
  Allocate(k2mask,nfft);
  Allocate(k2invmask,nfft);
  Allocate(temp,nmode);
  Allocate(mask,nfft);
  Allocate(count,nshells);
  
  k2maskij.Dimension(Nxb,Nyp,k2mask);
  
  for(unsigned i=0; i < nmode; i++) temp[i]=I*CartesianMode[i].X();
  Geometry->Pad(mask,temp);
  for(unsigned i=0; i < nfft; i++) kxmask[i]=mask[i].im;
  
  for(unsigned i=0; i < nmode; i++) temp[i]=I*CartesianMode[i].Y();
  Geometry->Pad(mask,temp);
  for(unsigned i=0; i < nfft; i++) {
    kymask[i]=mask[i].im;
    Real k2=kxmask[i]*kxmask[i]+kymask[i]*kymask[i];
    k2mask[i]=k2;
    k2invmask[i]=k2 ? 1.0/k2 : 0.0;
  }
  
  u.Set(Y[VEL]);
		
  // Initialize arrays with zero boundary conditions
  u=0.0;
	
  for(unsigned i=0; i < Nxb; i++) {
    for(unsigned j=0; j < Nyb; j++) {
#if 0     
      Real x=Xof(i);
      Real y=Yof(j);
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
    rc->fft(Sk[s]);
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
    array2<Real> Sss=Ss[s];
    cr->fft(Sk[s]);
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
  
  for(unsigned s=0; s < nspecies; s++) cr->fft(FMk[s]);
  
  
  for(unsigned i=0; i < NY[EK]; i++) Y[EK][i]=0.0;
//  for(unsigned i=0; i < NY[PX]; i++) Y[PX][i]=0.0;
//  for(unsigned i=0; i < NY[PV]; i++) Y[PV][i]=0.0;
  
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
  }
  
  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
  errno=0;
    
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
}

void DNS::Spectrum(vector& S, const vector& y)
{
  u.Set(y);
  
  for(unsigned s=0; s < nspecies; s++) {
    for(unsigned i=0; i < Nxb; i++) {
      for(unsigned j=0; j < Nyb; j++) {
	Ss(s,i,j)=u(i,j,s)*Nxybinv;
      }
    }
    rc->fft(Sk[s]);
  }
  
  for(unsigned K=0; K < nshells; K++) {
    count[K]=0;
    S[K]=0.0;
  }
				
  // Compute instantaneous angular average over circular shell.
		
  for(unsigned s=0; s < nspecies; s++) {
    for(unsigned i=0; i < Nxb; i++) {
      int kx=i-Nxb2;
      rvector k2maski=k2maskij[i];
      for(unsigned j=(kx > 0) ? 0 : 1; j < Nyp; j++) {
	if(k2maski[j]) {
	  int K2=kx*kx+j*j;
	  unsigned K=(unsigned)(sqrt(K2)+0.5);
	  count[K]++;
	  S[K] += abs2(Sk(s,i,j));
	}
      }
    }
  }
  
  for(unsigned K=0; K < nshells; K++)
    if(count[K]) S[K] *= twopi*K/count[K];
}

void DNS::Output(int it)
{
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
}

void DNS::OutFrame(int)
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
    array2<Real> ui=u[i];
    for(unsigned j=0; j < Nyb; j++) {
    rvector uij=ui[j];
      for(unsigned s=0; s < nspecies; s++) 
	E += uij[s]*uij[s];
      Real w=Vorticity(i,j);
      Z += w*w;
    }
  }
  
  Real volume=(xmax-xmin)*(ymax-ymin)/(Nxb*Nyb);
  // Real volume=Nxybinv <- to compute energy density instead of total energy
  Real factor=0.5*volume;
  
  E *= factor;
  Z *= factor;
}

void DNS::FinalOutput()
{
  Real E,Z;
  ComputeInvariants(E,Z);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
}

void DNS::Source(const vector2& Src, const vector2& Y, double)
{
  u.Set(Y[VEL]);
  S.Set(Src[VEL]);
  
  for(unsigned s=0; s < nspecies; s++) {
    
    for(unsigned i=0; i < Nxb; i++) {
      array2<Real> ui=u[i];
      rvector usi=us[i];
      for(unsigned j=0; j < Nyb; j++) {
	usi[j]=ui(j,s);
      }
    }
  
    rc->fft(uk);
  
    for(unsigned i=0; i < nfft; i++) {
      ikxu(i).re=-uk(i).im*kxmask[i];
      ikxu(i).im=uk(i).re*kxmask[i];
    }
    
    cr->fft(ikxu);
  
    for(unsigned i=0; i < nfft; i++) {
      ikyu(i).re=-uk(i).im*kymask[i];
      ikyu(i).im=uk(i).re*kymask[i];
    }
    cr->fft(ikyu);
  
    for(unsigned i=0; i < Nxb; i++) {
      array2<Real> ui=u[i];
      rvector dudxi=dudx[i];
      rvector dudyi=dudy[i];
      for(unsigned j=0; j < Nyb; j++) {
	Ss(s,i,j)=-(ui(j,0)*dudxi[j]+ui(j,1)*dudyi[j])*Nxybinv;
      }
    }
    
    array2<Complex> Sks=Sk[s];
    rc->fft(Sks);
    
    for(unsigned i=0; i < nfft; i++) {
      if(k2mask[i]) {
	Sks(i)=(Sks(i)-nu*k2mask[i]*uk(i))*Nxybinv;
      } else Sks(i)=0.0;
    }
  }
  
  array2<Complex> Sk0=Sk[0];
  array2<Complex> Sk1=Sk[1];
  for(unsigned i=0; i < nfft; i++) {
    Real kx=kxmask[i];
    Real ky=kymask[i];
    // Calculate -i*P
    Complex miP=(kx*Sk0(i)+ky*Sk1(i))*k2invmask[i];
    Sk0(i) -= kx*miP;
    Sk1(i) -= ky*miP;
  }
  
  for(unsigned s=0; s < nspecies; s++) {
    array2<Real> Sss=Ss[s];
    cr->fft(Sk[s]);
    for(unsigned i=0; i < Nxb; i++) {
      array2<Real> Si=S[i];
      rvector Sssi=Sss[i];
      for(unsigned j=0; j < Nyb; j++) {
	Si(j,s)=Sssi[j];
      }
    }
  }
 
  Spectrum(Src[EK],Y[VEL]);
  
//  Src[PX]=0.0;
//  Src[PV]=0.0;
}

void DNS::Stochastic(const vector2&Y, double, double)
{
  u.Set(Y[VEL]);
  Real factor=sqrt(2.0*dt)*rand_gauss();
  
  for(unsigned s=0; s < nspecies; s++) {
    array2<Real> ForceMasks=ForceMask[s];
    for(unsigned i=0; i < Nxb; i++) {
      array2<Real> ui=u[i];
      rvector ForceMasksi=ForceMasks[i];
      for(unsigned j=0; j < Nyb; j++) {
	ui(j,s) += factor*ForceMasksi[j];
      }
    }
  }
}
