#include "options.h"
#include "kernel.h"
#include "phin3.h"
#include "Grid3N.h"
#include "fft.h"

using namespace Array;

const char *method="HW";
const char *integrator="PC";

// Vocabulary
Real Dphi=0.0, Dn=0.0; // Laplacian diffusion
Real mu=0.25; // Velocity-dependent diffusion
Real shear=0.0;
Real ICphi=0.25, ICn=0.25;
int mx=1, my=1, mz=1;
int mx2=0, my2=0, mz2=0;
int randomIC=0;
int xperiodic=0;
int yperiodic=0;
int zperiodic=0;
Real xmin=0.0,xmax=1.0;
Real ymin=0.0,ymax=1.0;
Real zmin=0.0,zmax=1.0;
Real nmin=1.0, nmax=1.0;
Real nonlinear=1.0;
Real C=1.0;
Real implicit_factor=0.5;
Real explicit_factor=0.5;
int interpolate=0; 	// Linearly interpolate solver initial guess temporally
int nonlocal=1;
int nhyperviscosity=1;
int pseudospectral=0;
int FDdissipation=0;
int fftout=0;
int spectrum=0;
int spectrum2d=0;
int Nbox=1;
int xlevel=4;
int ylevel=4;
int zlevel=4;
int byte=1;
int movie=1;
int twisted=0;
int nsections=0;	// Requested number of xy cross sections
int niterations=1;
int npresmooth=0;
int npostsmooth=1;
int niterations2=1;
int npresmooth2=0;
int npostsmooth2=1;
int ziterations=1;
int defectstat=0;
Real convrate=0.0;

// Local Variables;
static int nlevel;
static int skip=1;
static int sections;
static int singular;
static int periodic3;
Real vstar=0.0;
static U rmin=REAL_MAX, rmax=0.0;
static int imin=INT_MAX, imax=0;
static const Real twelfth=1.0/12.0;
Array3<Real> *rjoffsetm, *rjoffsetp;
Array4<Real> rjoffset;
	
static U step;
static U minval=REAL_MAX;
static U maxval=-REAL_MAX;

template<class T>
class FineGrid : public Grid3<T> {
public:	
  Limits XMeshRange(),YMeshRange(),ZMeshRange();
  void BoundaryConditions(const Array3<T>& u);
  void Defect(const Array3<T>& d0, const Array3<T>& u, const Array3<T>& f)
  {};
  void Smooth(const Array3<T>& u, const Array3<T>& f) {};
};

class HWVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Hasagawa-Wakatani Multigrid";}
  const char *Abbrev() {return "hw";}
  HWVocabulary();
};

class HW : public ProblemBase {
  Real hxinv,hyinv,hzinv,volume;
  Real hxyinv;
  int Nx,Ny,Nz;
  int Nxi,Nyi,Nzi,Nzp;
  int Nxi2,Nyi2,Nzi2;
  int log2Nxi,log2Nyi,log2Nzi;
  int Nfft;
  Real kxmin,kymin,kzmin;
	
  Complex *phifft,*nfft;
  Real *phifftr,*nfftr;
	
  Array4<U> u,S,f;
  Array3<U> w,L2w,u1;
	
  Array3<Real> Sphi,Sn;
  Array3<Real> Hx,Hy;
  Array3<Real> F1x,F1y;
  Array3<Real> F2x,F2y;
  Array3<Real> F3x,F3y;
  MultiGrid< Poisson3Np<U> > *Multigrid;
  Poisson3N<U> *grid;
  FineGrid<Real> gridreal;
  Real *z;
  U defect0,defect;
  Array1<Real> phi_prof,n_prof;
  Array1<Real> phi_corr,n_corr;
  Real dtold;
public:
  virtual ~HW() {}
  void InitialConditions();
  void Initialize();
  void Output(int it);
  void Spectrum(Complex *fft, ofstream& s, const char *text, int norm);
  void OutFrame(int b, int k, int pass);
  void Source(Var *, Var *, double);
  void Transform(Var *, double, double, Var *&);
  void BackTransform(Var *, double, double, Var *);
  void ComputeInvariants(Real& E, Real& Z, Real& N, Real & I, Real& G,
			 Real& M);
  void SetBeta(double dt);
  void Setuold(Var *yi);
  void FFT(Var *y, Array4<U> &u);
  void FFTinv(Array4<U> &u, Var *y);
  void Laplacian(const Array3<U>& u0, const Array3<U>& f, int phionly=0);
  void Divergence(const Array3<Real>& Fx, const Array3<Real>& Fy,
		  const Array3<Real>& Div);
  inline void VorticityBoundaryConditions(const Array3<U>& v);
  void ShearOffset(const Array2<Real>& rjoffset, int level, Real sign);
};

// Note: this function evaluates the negative of the Laplacian
// (a positive semi-definite operator) to third-order accuracy
// (fourth-order in the absence of shear).
// u1 is a work array of same dimensions as u

static const int PHI_ONLY=1;

void HW::Laplacian(const Array3<U>& u0, const Array3<U>& f0, int phionly)
{
  Real coeffx1=twelfth*hxinv;
	
  if(alpha) {
    for(int i=1; i <= Nxi; i++) {
      Array2<U> um=u0[i-1], up=u0[i+1];
      Array2<U> uP=(i+2 < Nx) ? u0[i+2] : (xperiodic ? u0[2] : u0[Nx-1]);
      Array2<U> uM=(i-2 >= 0) ? u0[i-2] : (xperiodic ? u0[Nx-3] : u0[0]);
      Array2<U> u1i=u1[i];
      for(int j=1; j <= Nyi; j++) {
	Array1<U>::opt umj=um[j], upj=up[j], uPj=uP[j], uMj=uM[j], u1ij=u1i[j];
	for(int k=1; k <= Nzi; k++) {
	  u1ij[k].phi=coeffx1*
	    (uMj[k].phi-8.0*(umj[k].phi-upj[k].phi)-uPj[k].phi);
	  if(!phionly) u1ij[k].n=coeffx1*
			 (uMj[k].n-8.0*(umj[k].n-upj[k].n)-uPj[k].n);
	}
      }
    }
    VorticityBoundaryConditions(u1);
  }
	
  Real alpha2=alpha*alpha;
  Real hy2inv=hyinv*hyinv;
	
  Real coeffx2=-twelfth*hxinv*hxinv;
  for(int k=1; k <= Nzi; k++) 
    coeffy2[k]=-twelfth*hy2inv*(1.0+alpha2*z[k]*z[k]);
  Real coeffy1=-twelfth*hyinv*2.0*alpha;
		
  for(int i=1; i <= Nxi; i++) {
    Array2<U> um=u0[i-1], ui=u0[i], up=u0[i+1];
    Array2<U> uP=(i+2 < Nx) ? u0[i+2] : (xperiodic ? u0[2] : u0[Nx-1]);
    Array2<U> uM=(i-2 >= 0) ? u0[i-2] : (xperiodic ? u0[Nx-3] : u0[0]);
		
    Array2<U> u1i=u1[i];
    Array2<U> fi=f0[i];
    for(int j=1; j <= Nyi; j++) {
      Array1<U>::opt umj=um[j], uim=ui[j-1], uij=ui[j], uip=ui[j+1];
      Array1<U>::opt upj=up[j], uPj=uP[j], uMj=uM[j];
      Array1<U>::opt uiP=(j+2 < Ny) ? ui[j+2] :
	(yperiodic ? ui[2] : ui[Ny-1]);
      Array1<U>::opt uiM=(j-2 >= 0) ? ui[j-2] :
	(yperiodic ? ui[Ny-3] : ui[0]);
			
      Array1<U>::opt u1im=u1i[j-1], u1ip=u1i[j+1];
      Array1<U>::opt u1iP=(j+2 < Ny) ? u1i[j+2] :
	(yperiodic ? u1i[2] : u1i[Ny-1]);
      Array1<U>::opt u1iM=(j-2 >= 0) ? u1i[j-2] :
	(yperiodic ? u1i[Ny-3] : u1i[0]);
			
      Array1<U>::opt fij=fi[j];
      for(int k=1; k <= Nzi; k++) {
	fij[k].phi=coeffx2*(-uMj[k].phi+16.0*umj[k].phi-
			    30.0*uij[k].phi+16.0*upj[k].phi-
			    uPj[k].phi)+
	  coeffy2[k]*(-uiM[k].phi+16.0*uim[k].phi-30.0*uij[k].phi+
		      16.0*uip[k].phi-uiP[k].phi)+
	  coeffy1*z[k]*(u1iM[k].phi-8.0*(u1im[k].phi-u1ip[k].phi)-
			u1iP[k].phi);
	fij[k].n=phionly ? uij[k].n : 
	  coeffx2*(-uMj[k].n+16.0*umj[k].n-
		   30.0*uij[k].n+16.0*upj[k].n-
		   uPj[k].n)+
	  coeffy2[k]*(-uiM[k].n+16.0*uim[k].n-30.0*uij[k].n+
		      16.0*uip[k].n-uiP[k].n)+
	  coeffy1*z[k]*(u1iM[k].n-8.0*(u1im[k].n-u1ip[k].n)-
			u1iP[k].n);
      }
    }
  }
}


HWVocabulary::HWVocabulary()
{
  Vocabulary=this;
	
  VOCAB(Dphi,0.0,REAL_MAX,"");
  VOCAB(Dn,0.0,REAL_MAX,"");
  VOCAB(mu,0.0,REAL_MAX,"");
  VOCAB(shear,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICphi,-REAL_MAX,REAL_MAX,"");
  VOCAB(ICn,-REAL_MAX,REAL_MAX,"");
  VOCAB(nonlinear,0.0,REAL_MAX,"");
  VOCAB(C,0.0,REAL_MAX,"");
	
  VOCAB(xmin,-REAL_MAX,REAL_MAX,"");
  VOCAB(xmax,-REAL_MAX,REAL_MAX,"");
  VOCAB(ymin,-REAL_MAX,REAL_MAX,"");
  VOCAB(ymax,-REAL_MAX,REAL_MAX,"");
  VOCAB(zmin,-REAL_MAX,REAL_MAX,"");
  VOCAB(zmax,-REAL_MAX,REAL_MAX,"");
	
  VOCAB(nmin,0.0,REAL_MAX,"");
  VOCAB(nmax,0.0,REAL_MAX,"");
	
  VOCAB(Nbox,1,INT_MAX,"");
  VOCAB(mx,0,INT_MAX,"");
  VOCAB(my,0,INT_MAX,"");
  VOCAB(mz,0,INT_MAX,"");
  VOCAB(mx2,0,INT_MAX,"");
  VOCAB(my2,0,INT_MAX,"");
  VOCAB(mz2,0,INT_MAX,"");
  VOCAB(randomIC,0,1,"");
  VOCAB(xperiodic,0,1,"");
  VOCAB(yperiodic,0,1,"");
  VOCAB(zperiodic,0,1,"");
	
  VOCAB(explicit_factor,0.0,REAL_MAX,"");
  VOCAB(implicit_factor,0.0,REAL_MAX,"");
  VOCAB(interpolate,0,1,"");
  VOCAB(nonlocal,0,1,"");
  VOCAB(nhyperviscosity,1,INT_MAX,"");
  VOCAB(pseudospectral,0,1,"");
  VOCAB(FDdissipation,0,1,"");
  VOCAB(fftout,0,1,"");
  VOCAB(spectrum,0,1,"");
  VOCAB(spectrum2d,0,1,"");
  VOCAB(xlevel,1,INT_MAX,"");
  VOCAB(ylevel,2,INT_MAX,"");
  VOCAB(zlevel,1,INT_MAX,"");
  VOCAB(byte,0,1,"");
  VOCAB(movie,0,1,"");
  VOCAB(twisted,0,1,"");
  VOCAB(nsections,0,INT_MAX,"");
  VOCAB(niterations,1,INT_MAX,"");
  VOCAB(niterations2,0,INT_MAX,"");
  VOCAB(ziterations,1,INT_MAX,"");
  VOCAB(npresmooth,0,INT_MAX,"");
  VOCAB(npostsmooth,0,INT_MAX,"");
  VOCAB(npresmooth2,0,INT_MAX,"");
  VOCAB(npostsmooth2,0,INT_MAX,"");
  VOCAB(defectstat,0,1,"");
  VOCAB(convrate,0.0,REAL_MAX,"");
	
  VOCAB_NODUMP(vstar,0.0,REAL_MAX,""); // Obsolete
	
  METHOD(HW);
}

HWVocabulary HW_Vocabulary;

ofstream ft,fevt,fprofile,fcorr,fdefect,fekvk,fnkvk,fu;
oxstream fphi,fn,fphiRe,fphiIm;

Limits xlim,ylim,zlim;

template<class T>
Limits Poisson3<T>::XMeshRange() {return xlim;}
template<class T>
Limits Poisson3<T>::YMeshRange() {return ylim;}
template<class T>
Limits Poisson3<T>::ZMeshRange() {return zlim;}
template<class T>
Limits Poisson2<T>::XMeshRange() {return xlim;}
template<class T>
Limits Poisson2<T>::YMeshRange() {return ylim;}

template<class T>
Limits FineGrid<T>::XMeshRange() {return xlim;}
template<class T>
Limits FineGrid<T>::YMeshRange() {return ylim;}
template<class T>
Limits FineGrid<T>::ZMeshRange() {return zlim;}

template<class T>
void Poisson3N<T>::SetNumber() {ndomain=Nbox;}

template<class T>
void Poisson3N<T>::Smooth(const Array4<T>& u, const Array4<T>& f) {
  for(int n=0; n < ndomain; n++) {
    box=n;
    grid3.XYZebra(u[n],f[n]);
  }
		
  for(int l=0; l < ziterations; l++) {
    BoundaryConditions(u);
    for(int n=0; n < ndomain; n++) {
      box=n;
      grid3.ZLexicographical(u[n],f[n]);
    }
  }
}		
	
template<class T>
void Poisson3<T>::BoundaryConditions(const Array3<T>& u)
{
  // NB: This routine is used only to update perpendicular BC's in Smooth().
  if(xperiodic) XPeriodic(u);
  if(yperiodic) YPeriodic(u);
}

template<class T>
void Poisson3N<T>::BoundaryConditions(const Array4<T>& u)
{
  int ndomainm1=ndomain-1;
	
  for(int b=0; b < ndomain; b++) {
    Array3<T> ub=u[b];
    if(xperiodic) grid3.XPeriodic(ub);
    if(yperiodic) grid3.YPeriodic(ub);
	
    int nxbc=Nxbc();
    int nybc=Nybc();
    int nz=Nz();
    Array3<T> up=u[(b > 0) ? b-1 : ndomain-1];
    Array3<T> un=u[(b < ndomainm1) ? b+1 : 0];

    if(alpha) {
      Array2<Real> Rjoffsetp=rjoffsetp[level][(b > 0) ? 0 : 1];
      Array2<Real> Rjoffsetm=rjoffsetm[level][(b < ndomainm1) ? 0 : 1];
      for(int i=0; i < nxbc; i++) {
	Array2<T> upi=up[i], ubi=ub[i], uni=un[i];
	Array1<Real>::opt Rjoffsetmi=Rjoffsetm[i], Rjoffsetpi=Rjoffsetp[i];
	for(int j=0; j < nybc; j++) {
	  Array1<T>::opt ubij=ubi[j];
					
	  Real r=Rjoffsetpi[j];
	  int joff=(int) r;
	  Array1<T>::opt upij1=upi[joff], upij2=upi[joff+1];
	  if(zperiodic || b > 0)
	    ubij[0]=upij1[nz]+(upij2[nz]-upij1[nz])*(r-joff);
					
	  r=Rjoffsetmi[j];
	  joff=(int) r;
	  Array1<T>::opt unij1=uni[joff], unij2=uni[joff+1];
	  if(zperiodic || b < ndomainm1)
	    ubij[nz+1]=unij1[1]+(unij2[1]-unij1[1])*(r-joff);
	}
      }
    } else {
      for(int i=0; i < nxbc; i++) {
	Array2<T> upi=up[i], ubi=ub[i], uni=un[i];
	for(int j=0; j < nybc; j++) {
	  Array1<T>::opt unij=uni[j], ubij=ubi[j], upij=upi[j];
	  if(zperiodic || b > 0) ubij[0]=upij[nz];
	  if(zperiodic || b < ndomainm1) ubij[nz+1]=unij[1];
	}
      }
    }	
  }	
}

inline void FineGrid<Real>::BoundaryConditions(const Array3<Real>& u)
{
  xperiodic ? XPeriodic(u) : XConstant(u);
  yperiodic ? YPeriodic(u) : YConstant(u);
}

inline void HW::VorticityBoundaryConditions(const Array3<U>& v)
{
  xperiodic ? grid->Poisson3().XPeriodic(v) : grid->Poisson3().XConstant(v);
  yperiodic ? grid->Poisson3().YPeriodic(v) : grid->Poisson3().YConstant(v);
}

void HW::ShearOffset(const Array2<Real>& rjoffset, int level, Real deltaz)
{
  Poisson3N<U> *Grid=&Multigrid->Grid(level);
  int nxbc=Grid->Nxbc();
  int nybc=Grid->Nybc();
  int mody=Grid->Ny();
  Real Hyinv=1.0/Grid->Hy();
  Real factor=alpha*deltaz*Hyinv;
  for(int i=0; i < nxbc; i++) {
    Real offset=factor*Grid->X(i);
    Array1<Real>::opt rjoffseti=rjoffset[i];
    for(int j=0; j < nybc; j++) {
      Real r=j+offset;
      while (r < 0 ) r += mody;
      while (r > mody) r -= mody;
      rjoffseti[j]=r;
    }
  }
}

void dispersion(int mx, int my, int mz, Real& kx, Real& ky, Real& kz,
		Complex& Omegap, Complex& Omegam,
		Complex& phasep, Complex& phasem)
{
  const Complex I(0.0,1.0);
  kx=mx*twopi/(xmax-xmin);
  ky=my*twopi/(ymax-ymin);
  kz=mz*twopi/(Nbox*(zmax-zmin));
  Real kperp2=(kx*kx+ky*ky);
  if(kperp2 == 0.0) {Omegap=Omegam=phasep=phasem=0.0; return;}
  Real a=kz*kz*C;
  Real b=1.0+1.0/kperp2;
  Real kperp4=kperp2*kperp2;
	
  Real khyper;
  if(nhyperviscosity == 1) khyper=1.0;
  else {
    khyper=kperp2;
    for(int n=2; n < nhyperviscosity; n++) khyper *= kperp2;
  }
	
  Real Rphi=khyper*Dphi;
  Real Rn=khyper*Dn;
	
  Real B=a*b+(Rphi+Rn)*kperp2;
  Complex discr=sqrt(4*(a*((Rphi*kperp2+Rn)+I*ky*vstar/kperp2)+
			Rphi*Rn*kperp4)-B*B);
  Omegap=0.5*(-I*B+discr);
  Omegam=0.5*(-I*B-discr);
  cout << endl << "LINEAR EIGENVALUES: Omega(" << mx << ", " << my << ", "
       << mz << ") = " << Omegap;
  cout << endl << "                    Omega(" << mx << ", " << my << ", "
       << mz << ") = " << Omegam << endl; 
	
  phasep=phasem=1.0;
  if(kz != 0.0) {
    phasep += (-I*Omegap+Rphi*kperp2)*kperp2/(kz*kz);
    phasem += (-I*Omegam+Rphi*kperp2)*kperp2/(kz*kz);
  }
}

void HW::InitialConditions()
{
  int i;
	
  if(convrate || verbose > 2) defectstat=1;
  if(convrate && niterations == 1) {niterations=100;}
  zmax=zmin+(zmax-zmin)/Nbox;
	
  rmin=rmax=0.0;
  imin=imax=0;
	
  nlevel=max(max(xlevel,ylevel),zlevel);
	
  xlim=xperiodic ? Limits(xmin,xmax,Periodic,nlevel-xlevel) : 
    Limits(xmin,xmax,Dirichlet,nlevel-xlevel);
  ylim=yperiodic ? Limits(ymin,ymax,Periodic,nlevel-ylevel) :
    Limits(ymin,ymax,Dirichlet,nlevel-ylevel);
  zlim=zperiodic ? Limits(zmin,zmax,Periodic,nlevel-zlevel) :
    Limits(zmin,zmax,Dirichlet,nlevel-zlevel);
		
  singular=(xperiodic && yperiodic);
  periodic3=(singular && zperiodic);
	
  if(xperiodic) {
    if(nonlocal) msg(ERROR,"Nonlocal option requires xperiodic=0");
    if(shear) msg(OVERRIDE,"Shear specified with xperiodic=1");
  }
	
  if(pseudospectral && !periodic3)
    msg(ERROR, "Pseudospectral option requires periodic BC's");
	
  if(!pseudospectral) FDdissipation=1;
						   
  if(FDdissipation && nhyperviscosity > 1)
    msg(WARNING,
	"nhyperviscosity > 1 may be inaccurate with finite differencing"); 
	
  kxmin=twopi/(xmax-xmin);
  kymin=twopi/(ymax-ymin);
  kzmin=twopi/(zmax-zmin);
	
  Multigrid=new MultiGrid< Poisson3Np<U> >(nlevel);
  Multi2=new MultiGrid< Poisson2<Real> >(nlevel);
  gridreal.Initialize(nlevel-1);
	
  grid=&Multigrid->Fine();
	
  Nx=grid->Nxbc();
  Ny=grid->Nybc();
  Nz=grid->Nzbc();
	
  Nxi=grid->Nx();
  Nyi=grid->Ny();
  Nzi=grid->Nz();
	
  cout << endl << "GEOMETRY: (" << Nxi << " x " << Nyi << " x " << Nzi*Nbox
       << ")" << endl; 
	
  if(fftout || spectrum || pseudospectral) {
    for(log2Nxi=0; Nxi > (1 << log2Nxi); log2Nxi++);
    for(log2Nyi=0; Nyi > (1 << log2Nyi); log2Nyi++);
    for(log2Nzi=0; Nzi > (1 << log2Nzi); log2Nzi++);
		
    if(Nxi != (1 << log2Nxi) || Nyi != (1 << log2Nyi) || 
       Nzi != (1 << log2Nzi)) {
      const char *s=
	"Spectral data available only for power-of-2 geometries";
      if(pseudospectral) msg(ERROR,s);
      else msg(WARNING,s);
      fftout=0;
      spectrum=0;
    } else {
      Nxi2=Nxi/2, Nyi2=Nyi/2, Nzi2=Nzi/2;
		
      Nzp=Nzi/2+1;
      Nfft=Nxi*Nyi*Nzp;
		
      phifft=new Complex[Nfft];
      phifftr=(Real *) phifft;
      nfft=new Complex[Nfft];
      nfftr=(Real *) nfft;
		
      for(i=0; i < Nfft; i++) phifft[i]=nfft[i]=0.0;
    }
  }
	
  hxinv=1.0/grid->Hx();
  hyinv=1.0/grid->Hy();
  hzinv=1.0/grid->Hz();
	
  volume=grid->Hx()*grid->Hy()*grid->Hz();
	
  hxyinv=hxinv*hyinv;
	
  z=grid->Z();
	
  coeffy2=new Real[Nz-1];
		
  if(pseudospectral) {
    if(nonlocal) msg(ERROR,"Pseudospectral option requires nonlocal=0");
    if(shear) msg(ERROR,"Pseudospectral mode requires shear=0");
    if(C == 0.0) msg(ERROR,"Pseudospectral mode requires C != 0");
    ny=Nbox*Nfft*sizeof(UComplex)/sizeof(Var);
    y=new Var[ny];
    u.Allocate(Nbox,Nx,Ny,Nz);
    S.Allocate(Nbox,Nx,Ny,Nz);
    errmask=new int[ny];
  } else {
    ny=Nbox*Nx*Ny*Nz;
    y=new Var[ny];
    u.Dimension(Nbox,Nx,Ny,Nz);
    S.Dimension(Nbox,Nx,Ny,Nz);
    vtemp.Allocate(Nx,Ny);
    ftemp.Allocate(Nx,Ny);
		
    int m=max(Ny,Nz);
    coeffa=new Real[m];
    coeffb=new Real[m];
    coeffc=new Real[m];
    g=new Real[Nz];
    zavg=new Var[Nz];

    u.Set(y);
		
    nold=new Array2<Real>[nlevel];
    for(i=1; i < nlevel; i++) {
      int nx0=Multi2->Grid(i).Nxbc();
      int ny0=Multi2->Grid(i).Nybc();
      nold[i].Allocate(nx0,ny0);
    }
  }
	
  uold=new Array4<U>[nlevel];
	
  for(i=((nlevel==1) ? 0 : 1); i < (nonlocal ? nlevel-1 : nlevel); i++) {
    int nx0=Multigrid->Grid(i).Nxbc();
    int ny0=Multigrid->Grid(i).Nybc();
    int nz0=Multigrid->Grid(i).Nzbc();
    uold[i].Allocate(Nbox,nx0,ny0,nz0);
    if(!nonlocal) uold[i]=1.0;
  }
	
  if(nonlocal) uold[nlevel-1].Dimension(Nbox,Nx,Ny,Nz);
	
  vstar=(nmax-nmin)/(xmax-xmin);
	
  // For linear dispersion relation
	
  cout << endl << "INVERSE DENSITY SCALE LENGTH: " << vstar << endl;

  if((FDdissipation) && (Dphi || Dn)) {
    L2w.Allocate(Nx,Ny,Nz);
    L2w=0.0;
  }
	
  f.Allocate(Nbox,Nx,Ny,Nz);
  w.Dimension(Nx,Ny,Nz,f);
	
  u1.Allocate(Nx,Ny,Nz);
	
  Sn.Allocate(Nx,Ny,Nz);
  Sphi.Allocate(Nx,Ny,Nz);
	
  Hx.Allocate(Nx,Ny,Nz);
  Hy.Allocate(Nx,Ny,Nz);
	
  F1x.Dimension(Nx,Ny,Nz,Sphi); // Memory optimization
  F1y.Allocate(Nx,Ny,Nz);
	
  F2x.Allocate(Nx,Ny,Nz);
  F2y.Allocate(Nx,Ny,Nz);
	
  F3x.Allocate(Nx,Ny,Nz);
  F3y.Allocate(Nx,Ny,Nz);
	
  phi_prof.Allocate(Nxi);
  n_prof.Allocate(Nxi);
	
  phi_corr.Allocate(Nz);
  n_corr.Allocate(Nz);
	
  f=0.0;
  u1=0.0;
	
  Sn=0.0;
  Sphi=0.0;
	
  F1y=0.0;
	
  F2x=0.0;
  F2y=0.0;
	
  F3x=0.0;
  F3y=0.0;
	
  alpha=shear;
	
  int twomodes=(mx2 || my2 || mz2);
  Real kx,ky,kz;
  Complex Omega1p,Omega1m,phase1p,phase1m;
  Real kx2,ky2,kz2;
  Complex Omega2p,Omega2m,phase2p,phase2m;
	
  int offset=1;
  u=0.0; // Initialize with zero boundary conditions
	
  if(!randomIC) {
    dispersion(mx,my,mz,kx,ky,kz,Omega1p,Omega1m,phase1p,phase1m);
    dispersion(mx2,my2,mz2,kx2,ky2,kz2,Omega2p,Omega2m,phase2p,phase2m);
  }
	
  Real deltaz=zmax-zmin;
	
  int b;
	
  for(b=0; b < Nbox; b++) {
    Real zoff=b*deltaz;
    for(i=offset; i < Nx-offset; i++) {
      Real X=grid->X(i);
      for(int j=offset; j < Ny-offset; j++) {
	Real Y0=grid->Y(j);
	for(int k=offset; k < Nz-offset; k++) {
	  Real phi,n;
	  if(randomIC) {phi=2.0*drand()-1.0; n=2.0*drand()-1.0;}
	  else {
	    Real Z=z[k]+zoff;
	    Real Y=Y0-alpha*Z*X;
	    Real arg=kx*X+ky*Y+kz*Z;
	    phi=cos(arg);
	    n=phi*phase1p.re-sin(arg)*phase1p.im;
	    if(twomodes) {
	      Real arg=kx2*X+ky2*Y+kz2*Z;
	      Real phi2=cos(arg);
	      phi += phi2;
	      n += phi2*phase2p.re-sin(arg)*phase2p.im;
	    }
	  }
	  u(b,i,j,k)=U(ICphi*phi,ICn*n);
	}
      }
    }
	
    // Add mean density gradient
    if(nonlocal) {
      int Nxm1=Nx-1;
      Real nfactor=(nmax-nmin)/Nxm1;
      for(i=0; i < Nx; i++) {
	Real nmean=nmin+nfactor*(Nxm1-i);
	for(int j=0; j < Ny; j++) {
	  for(int k=0; k < Nz; k++) {
	    u(b,i,j,k).n += nmean;
	    if(u(b,i,j,k).n <= 0)
	      msg(ERROR,"Initial density is not positive");
	  }
	}
      }
    }
  }
	
  if(alpha) {
    if((zperiodic || Nbox > 1) && !yperiodic) 
      msg(ERROR,
	  "Shear with zperiodic=1 or Nbox > 1 requires yperiodic=1");
    rjoffsetm=new Array3<Real>[nlevel];
    rjoffsetp=new Array3<Real>[nlevel];
		
    for(int level=0; level < nlevel; level++) {
      Poisson3N<U> *Grid=&Multigrid->Grid(level);
      rjoffsetm[level].Allocate(2,Grid->Nxbc(),Grid->Nybc());
      rjoffsetp[level].Allocate(2,Grid->Nxbc(),Grid->Nybc());
      Real deltaz=zmax-zmin;
      ShearOffset(rjoffsetm[level][0],level,-deltaz);
      ShearOffset(rjoffsetp[level][0],level,deltaz);
      ShearOffset(rjoffsetm[level][1],level,-Nbox*deltaz);
      ShearOffset(rjoffsetp[level][1],level,Nbox*deltaz);
    }
		
  }
	
  rjoffset.Allocate(Nbox,Nz,Nx,Ny);
  for(b=0; b < Nbox; b++) {
    Real zoff=b*deltaz;
    for(int k=0; k < Nz; k++) {
      ShearOffset(rjoffset[b][k],nlevel-1,z[k]+zoff);
    }
  }
	
  grid->BoundaryConditions(u);
  if(nonlocal) Setuold(u);
		
  if(pseudospectral) FFT(y,u);
	
  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");
  open_output(fdefect,dirsep,"defect");
  if(fftout) {
    open_output(fphiRe,dirsep,"phifft.re");
    open_output(fphiIm,dirsep,"phifft.im");
  }
  if(!spectrum) {
    open_output(fekvk,dirsep,"ekvk");
    fekvk.close();
    open_output(fnkvk,dirsep,"nkvk");
    fnkvk.close();
  }
  if(output) open_output(fu,dirsep,"u");
  if(movie) {
    open_output(fphi,dirsep,"phi");
    open_output(fn,dirsep,"n");
		
    // Compute the actual number of xy cross sections.
    if(nsections == 0) nsections=Nzi;
    if(nsections == 1) sections=1;
    else {
      div_t d;
      d=div(Nzi-1,nsections-1);
      skip=(d.rem == 0) ? d.quot : d.quot+1;
      sections=0;
      for(int k=1; k <= Nzi; k += max(1,min(skip,Nzi-k))) sections++; 
    }
    sections *= Nbox;
  }
}

void HW::Initialize()
{
  fevt << "#   t\t\t E\t\t Z\t\t N\t\t I\t\t G\t\t M" << endl;
  fdefect << "t\t\t rmin\t\t rmax\t\t imin\t\t imax" << endl;
}

void HW::SetBeta(double dt)
{
  beta=C*dt;
}

void HW::FFT(Var *Y, Array4<U> &u0)
{
  int i,Nzip2=Nzi+2;
	
  for(int b=0; b < Nbox; b++) {
    Array3<U> ub=u0[b];
    for(i=0; i < Nxi; i++) {
      Array2<U> ui=ub[i+1];
      int Nyii=Nyi*i;
      for(int j=0; j < Nyi; j++) {
	Array1<U>::opt uij=ui[j+1];
	int Nij=Nzip2*(Nyii+j);
	for(int k=0; k < Nzi; k++) {
	  phifftr[Nij+k]=uij[k+1].phi;
	  nfftr[Nij+k]=uij[k+1].n;
	}
      }
    }

    rcfft3d(phifft,log2Nxi,log2Nyi,log2Nzi,-1);
    rcfft3d(nfft,log2Nxi,log2Nyi,log2Nzi,-1);
	
    if(Y) {
      UComplex *Yc=(UComplex *) Y;
      int bNfft=b*Nfft;
      for(i=0; i < Nfft; i++) {
	Yc[bNfft+i].phi=phifft[i];
	Yc[bNfft+i].n=nfft[i];
      }
    }
  }
}

void HW::FFTinv(Array4<U> &u0, Var *Y)
{
  int i;
	
  for(int b=0; b < Nbox; b++) {
    if(Y) {
      UComplex *Yc=(UComplex *) Y;
      int bNfft=b*Nfft;
      for(i=0; i < Nfft; i++) {
	phifft[i]=Yc[bNfft+i].phi;
	nfft[i]=Yc[bNfft+i].n;
      }
    }	
	
    Real scale=1.0/(Nxi*Nyi*Nzi);
    crfft3d(phifft,log2Nxi,log2Nyi,log2Nzi,1,scale);
    crfft3d(nfft,log2Nxi,log2Nyi,log2Nzi,1,scale);
		
    int Nzip2=Nzi+2;
    Array3<U> ub=u0[b];
    for(i=0; i < Nxi; i++) {
      Array2<U> ui=ub[i+1];
      int Nyii=Nyi*i;
      for(int j=0; j < Nyi; j++) {
	Array1<U>::opt uij=ui[j+1];
	int Nij=Nzip2*(Nyii+j);
	for(int k=0; k < Nzi; k++) {
	  uij[k+1].phi=phifftr[Nij+k];
	  uij[k+1].n=nfftr[Nij+k];
	}
      }
    }
  }
  grid->BoundaryConditions(u0);
}

void HW::Setuold(Var *yi)
{
  uold[nlevel-1].Set(yi);
  for(int i=nlevel-2; i >= 1; i--) {
    Multigrid->Grid(i+1).Restrict(uold[i],uold[i+1]);
    Multigrid->Grid(i).BoundaryConditions(uold[i]);
  }
}

void HW::Transform(Var *Y, double, double dt, Var *&yi)
{
  if(yi) {
    if(interpolate) {
      for(unsigned int i=0; i < ny; i++) 
	yi[i]=Y[i]+dt/dtold*(Y[i]-yi[i]);
      dtold=dt;
    } else set(yi,Y,ny);
  } else {yi=new Var[ny]; set(yi,Y,ny); dtold=dt;}
	
  SetBeta(-explicit_factor*dt);
	
  if(pseudospectral) {
    int i;
    UComplex *Yc=(UComplex *) Y;
    for(int b=0; b < Nbox; b++) {
      int Nxib=Nxi*b;
      for(i=0; i < Nxi; i++) {
	Real kx=(i-Nxi2)*kxmin;
	Real kx2=kx*kx;
	int Nyii=Nyi*(Nxib+i);
	for(int j=0; j < Nyi; j++) {
	  Real ky=(j-Nyi2)*kymin;
	  Real kxy2=kx2+ky*ky;
	  int Nij=Nzp*(Nyii+j);
	  for(int k=0; k < Nzp; k++) {
	    Real kz=k*kzmin;
	    Real kz2=kz*kz;
	    Real betakz2=beta*kz2;
	    int index=Nij+k;
	    Complex phi=Yc[index].phi;
	    Yc[index].phi=phi*(kxy2+betakz2)-betakz2*Yc[index].n;
	    Yc[index].n=Yc[index].n*(1.0+betakz2)-betakz2*phi;
	  }
	}
      }
    }
    for(unsigned int m=0; m < ny; m++) errmask[m]=(Y[m] != 0.0);
  } else {
    if(nonlocal) Setuold(yi);
    u.Set(Y);
    grid->Lu(u,f);
    u=f;
  }
}

void HW::BackTransform(Var *Y, double, double dt, Var *yi)
{
  SetBeta(implicit_factor*dt);
	
  if(pseudospectral) {
    int kinit=(implicit_factor ? 1 : Nzp);
    UComplex *Yc=(UComplex *) Y;
    for(int b=0; b < Nbox; b++) {
      Yc[Nfft*b+Nzp*(Nyi*Nxi2+Nyi2)].phi=0.0;
      int Nxib=Nxi*b;
      for(int i=0; i < Nxi; i++) {
	int i0=i-Nxi2;
	Real kx=i0*kxmin;
	Real kx2=kx*kx;
	int Nyii=Nyi*(Nxib+i);
	for(int j=0; j < Nyi; j++) {
	  int j0=j-Nyi2;
	  Real ky=j0*kymin;
	  Real kxy2=kx2+ky*ky;
	  int Nij=Nzp*(Nyii+j);
	  for(int k=((i0 == 0 && j0 == 0) ? kinit : 0); k < Nzp; k++)
	    {
	      int index=Nij+k;
	      Real kz=k*kzmin;
	      Real kz2=kz*kz;
	      Real betakz2=beta*kz2;
	      Real ndenominv=1.0/(1.0+betakz2);
	      Real factor=betakz2*ndenominv;
	      Yc[index].phi=(Yc[index].phi+Yc[index].n*factor)/
		(kxy2+betakz2-betakz2*factor);
	      Yc[index].n=(Yc[index].n+betakz2*Yc[index].phi)*
		ndenominv;
	    }
	}
      }
    }
  } else {
    u.Set(Y);
    f=y;
    u.Load(yi);
    if(nonlocal) Setuold(yi);
		
    int it=grid->Solve(u,f,npresmooth,1,npostsmooth,singular,
		       niterations,convrate,defectstat,&defect0,&defect);

    if(verbose > 2) {
      cout << endl;
      grid->Report(defect0,0);
      grid->Report(defect,defect0,it);
    }
		
    if(!it)
      msg(ERROR,"Maximum number of iterations (%d) exceeeded in Solve",
	  niterations);
		
    if(periodic3) {
      // Subtract off null solution.
      Real sum=0.0;
      for(int b=0; b < Nbox; b++) {
	Array3<U> ub=u[b];
	for(int i=1; i <= Nxi; i++) {
	  Array2<U> ui=ub[i];
	  for(int j=1; j <= Nyi; j++) {
	    Array1<U>::opt uij=ui[j];
	    for(int k=1; k <= Nzi; k++) sum += uij[k].phi;
	  }
	}
      }				
      sum /= Nxi*Nyi*Nzi*Nbox;
      for(unsigned int i=0; i < ny; i++) Y[i].phi -= sum;
    }
		
    if(defectstat) {
      U ratio=divide0(defect0,defect);
      if(ratio < rmin) rmin=ratio;
      if(ratio > rmax) rmax=ratio;
    }
	
    if(it < imin) imin=it;
    if(it > imax) imax=it;
  }
}

void HW::Spectrum(Complex *fft, ofstream& s, const char *text, int norm)
{
  int K, Kmax=(int) (sqrt(Nxi2*Nxi2+Nyi2*Nyi2+Nzi2*Nzi2)+0.5);
  Real *sum=new Real[Kmax+1];
  int *count=new int[Kmax+1];
		
  for(K=0; K <= Kmax; K++) {
    count[K]=0;
    sum[K]=0.0;
  }
				
  // Compute angular average over circular or spherical shell.
  int kstart,kstop;
  if(spectrum2d) {kstart=Nzi2; kstop=Nzi2+1;} // circular shell at kz=0
  else {kstart=0; kstop=Nzi;}
		
  for(int i=0; i < Nxi; i++) {
    int kx=i-Nxi2;
    for(int j=0; j < Nyi; j++) {
      int ky=j-Nyi2;
      int kperp2=kx*kx+ky*ky;
      for(int k=kstart; k < kstop; k++) {
	int kp=k-Nzi2;
	int K2=kperp2+kp*kp;
	int K=(int)(sqrt(K2)+0.5);
	int ip=i, jp=j;
	if(kp < 0) {
	  if(i > 0) ip=Nxi-i;
	  if(j > 0) jp=Nyi-j;
	  if(k > 0) kp=-kp;
	  else kp=Nzi2;
	}
	count[K]++;
	sum[K] += (norm ? K2 : 1.0)*abs2(fft[Nzp*(Nyi*ip+jp)+kp]);
      }
    }
  }
		
  for(K=0; K <= Kmax; K++)
    if(count[K]) sum[K] *= 0.5/count[K] *
		   (spectrum2d ? 2.0*pi*K : 4.0*pi*(K*K+twelfth));
				
  open_output(s,dirsep,text,0);
  out_curve(s,t,"t");
  out_curve(s,sum,text,Kmax+1);
  out_curve(s,Nxi,"Nxi");
  out_curve(s,Nyi,"Nyi");
  out_curve(s,Nzi,"Nzi");
  s.close();
}

void HW::Output(int)
{
  Real E,Z,N,I,G,M;
  int k;
	
  if(pseudospectral) {
    FFTinv(u,y);
		
    // Output spectrum only for last box.
    UComplex *Yc=(UComplex *) y;
    for(int i=0; i < Nfft; i++) {
      phifft[i]=Yc[i].phi;
      nfft[i]=Yc[i].n;
    }
  } else {
    u.Set(y);
    grid->BoundaryConditions(u);
	
    if(fftout || spectrum) FFT(NULL,u);
  }
		
  ComputeInvariants(E,Z,N,I,G,M);
  fevt << t << "\t" << E << "\t" << Z << "\t" << N << "\t" << I << "\t"
       << G << "\t" << M << endl;
	
  if(defectstat && pseudospectral == 0)
    fdefect << t << "\t" << rmin << "\t" << rmax << "\t" << imin << "\t" 
	    << imax << endl;
	
  rmin=REAL_MAX; rmax=0.0;
  imin=INT_MAX; imax=0;
  unsigned int nu=Nx*Ny*Nz;
  if(output) out_curve(fu,u(),"u",nu);
	
  if(movie) {
    int Nx0=xperiodic ? Nxi : Nx;
    int Ny0=yperiodic ? Nyi : Ny;
		
    fphi << Nx0 << Ny0 << sections;
    fn << Nx0 << Ny0 << sections;
		
    for(int pass=0; pass < (byte ? 2 : 1); pass++) {
      if(pass == 1) {
	if (maxval == minval) step=0.0;
	else step=255.0/(maxval-minval);
      }
      for(int b=0; b < Nbox; b++) {
	if(nsections == 1) OutFrame(b,Nz/2,pass);
	else for(int k=1; k <= Nzi; k += max(1,min(skip,Nzi-k))) {
	  OutFrame(b,k,pass);
	}
      }
    }
  }
	
  if(fftout) {
    // Output spectrum only for last box.
    fphiRe << Nxi << Nyi << Nzp;
    fphiIm << Nxi << Nyi << Nzp;
				
    for(int k=0; k < Nzp; k++) {
      for(int j=Nyi-1; j >= 0; j--) {
	for(int i=0; i < Nxi; i++) {
	  fphiRe << (float) phifft[Nzp*(Nyi*i+j)+k].re;
	  fphiIm << (float) phifft[Nzp*(Nyi*i+j)+k].im;
	}
      }
    }
    fphiRe.flush();
    fphiIm.flush();
  }
	
  // Output profiles.
  for(int i=1; i <= Nxi; i++) {
    U sum=0.0;
    for(int b=0; b < Nbox; b++) {
      Array2<U> ubi=u[b][i];
      for(int j=1; j <= Nyi; j++) {
	Array1<U>::opt ubij=ubi[j];
	for(int k=1; k <= Nzi; k++) {
	  sum += ubij[k];
	}
      }
    }
    sum /= (Nyi*Nzi*Nbox);
    phi_prof[i-1]=sum.phi;
    n_prof[i-1]=sum.n;
  }
	
  open_output(fprofile,dirsep,"profile",0);
  out_curve(fprofile,grid->X()+1,"X",Nxi);
  out_curve(fprofile,phi_prof(),"phi",Nxi);
  out_curve(fprofile,n_prof(),"n",Nxi);
  fprofile.close();
	
  // Output correlation function.
  int count0=Nxi*Nyi*Nbox;
	
  for(k=0; k < Nzi; k++) {
    Real sum_phi=0.0;
    Real sum_n=0.0;
    for(int b=0; b < Nbox; b++) {
      Array3<U> ub=u[b];
      for(int i=1; i <= Nxi; i++) {
	Array2<U> ubi=ub[i];
	for(int j=1; j <= Nyi; j++) {
	  Array1<U>::opt ubij=ubi[j];
	  for(int kp=0; kp < Nz-k; kp++) {
	    sum_phi += ubij[kp].phi*ubij[k+kp].phi;
	    sum_n += ubij[kp].n*ubij[k+kp].n;
	  }
	}
      }
    }
    int count=count0*(Nz-k);
    sum_phi /= count;
    sum_n /= count;
    phi_corr[k] += sum_phi;
    n_corr[k] += sum_n;
  }
	
  for(k=1; k < Nzi; k++) {
    phi_corr[k] /= phi_corr[0];
    n_corr[k] /= n_corr[0];
  }
	
  phi_corr[0]=n_corr[0]=1.0;
	
  open_output(fcorr,dirsep,"corr",0);
  out_curve(fcorr,grid->Z(),"Z",Nzi);
  out_curve(fcorr,phi_corr(),"phi",Nzi);
  out_curve(fcorr,n_corr(),"n",Nzi);
  fcorr.close();
	
  if(spectrum) {
    // Output spectrum only for last box.
    Spectrum(phifft,fekvk,"ekvk",1); // Normalize phi^2 by K^2
    Spectrum(nfft,fnkvk,"nkvk",0);
  }
	
  ft << t << endl;
}

void HW::OutFrame(int b, int k, int pass)
{
  int xoffset=xperiodic ? 2 : 0;
  int yoffset=yperiodic ? 2 : 0;
	
  Array3<Real> rjoffsetb=rjoffset[b];
  Array3<U> ub=u[b];

  for(int j=Ny-1-yoffset; j >= 0; j--) {
    for(int i=0; i < Nx-xoffset; i++) {
      U uinterp;
			
      if(alpha && !twisted) {
	Real r=rjoffsetb(k,i,j);
	int joff=(int) r;
	uinterp=ub(i,joff,k)+(ub(i,joff+1,k)-ub(i,joff,k))*(r-joff);
      } else uinterp=ub(i,j,k);
		
      if(byte) {
	if(pass == 0) {
	  minval=min(minval,uinterp);
	  maxval=max(maxval,uinterp);
	} else {
	  U scaled;
	  if (step == 0.0) scaled=127;
	  else scaled=(uinterp-minval)*step+0.5;
	  fphi << (xbyte) (unsigned int) scaled.phi;
	  fn << (xbyte) (unsigned int) scaled.n;
	}
      } else {
	fphi << (float) uinterp.phi;
	fn << (float) uinterp.n;
      }
    }
  }
  fphi.flush();
  fn.flush();
	
  const char *text="Cannot write to movie file %s";
  if(!fphi) msg(ERROR,text,"phi");
  if(!fn) msg(ERROR,text,"n");
}

void HW::ComputeInvariants(Real& E, Real& Z, Real& N, Real& I, Real& G,
			   Real &M)
{
  E=Z=N=I=G=M=0.0;
  Real coeffy=0.5*hyinv;
	
  int Nxm1=Nx-1;
  Real nfactor=(nmax-nmin)/Nxm1;
  for(int b=0; b < Nbox; b++) {
    Array3<U> ub=u[b];
    Laplacian(ub,w,PHI_ONLY);
    for(int i=1; i <= Nxi; i++) {
      Real nmean=nmin+nfactor*(Nxm1-i);
      Real n0;
      if(nonlocal) n0=0.0;
      else {
	n0=nmean;
	nmean=0.0;
      }
      Array2<U> ui=ub[i], wi=w[i];
      for(int j=1; j <= Nyi; j++) {
	Array1<U>::opt uim=ui[j-1], uij=ui[j], uip=ui[j+1], wij=wi[j];
	for(int k=1; k <= Nzi; k++) {
	  Real vort=wij[k].phi;
	  Real n=uij[k].n;
	  E += vort*uij[k].phi;
	  Z += vort*vort;
	  N += (n-nmean)*(n-nmean);
	  I += n*vort;
	  G += n*(uim[k].phi-uip[k].phi)*coeffy;
	  M += n+n0;
	}
      }
    }
  }
	
  Real factor=0.5*volume;
  E *= factor;
  Z *= factor;
  N *= factor;
  I *= factor;
  G *= volume;
  M *= volume;
}

void HW::Divergence(const Array3<Real>& Fx, const Array3<Real>& Fy,
		    const Array3<Real>& D)
{
  gridreal.BoundaryConditions(Fx);
  gridreal.BoundaryConditions(Fy);
	
  Real coeffx12=twelfth*hxinv;
  Real coeffy12=twelfth*hyinv;
	
  for(int i=1; i <= Nxi; i++) {
    Array2<Real> Fxm=Fx[i-1], Fxp=Fx[i+1];
    Array2<Real> FxP=(i+2 < Nx) ? Fx[i+2] : (xperiodic ? Fx[2] : Fx[Nx-1]);
    Array2<Real> FxM=(i-2 >= 0) ? Fx[i-2] : (xperiodic ? Fx[Nx-3] : Fx[0]);
    Array2<Real> Di=D[i], Fyi=Fy[i];
		
    for(int j=1; j <= Nyi; j++) {
      Array1<Real>::opt Dij=Di[j]; 
      Array1<Real>::opt Fxmj=Fxm[j], Fxpj=Fxp[j], FxPj=FxP[j], FxMj=FxM[j];
			
      Array1<Real>::opt Fyim=Fyi[j-1], Fyip=Fyi[j+1];
      Array1<Real>::opt FyiP=(j+2 < Ny) ? Fyi[j+2] :
	(yperiodic ? Fyi[2] : Fyi[Ny-1]);
      Array1<Real>::opt FyiM=(j-2 >= 0) ? Fyi[j-2] :
	(yperiodic ? Fyi[Ny-3] : Fyi[0]);
			
      for(int k=1; k <= Nzi; k++) {
	Dij[k]=coeffx12*(FxMj[k]-8.0*(Fxmj[k]-Fxpj[k])-FxPj[k])+
	  coeffy12*(FyiM[k]-8.0*(Fyim[k]-Fyip[k])-FyiP[k]);
      }
    }
  }
}				

void HW::Source(Var *source, Var *Y, double)
{
  int i;
	
  if(pseudospectral) FFTinv(u,Y);
  else {
    u.Set(Y);
    grid->BoundaryConditions(u);
    S.Set(source);
  }
	
  Array4<U> uoldl=uold[nlevel-1];
	
  for(int b=0; b < Nbox; b++) {
    Array3<U> ub=u[b], Sb=S[b];
		
    if(!xperiodic) {
      // Enforce Neumann BC on phi in x-direction (not in solver)
      Array2<U> ub0=ub[0], ub1=ub[1], ubnx=ub[Nx-2], ubnx1=ub[Nx-1];
      for(int j=0; j < Ny; j++) {
	Array1<U>::opt ub0j=ub0[j], ub1j=ub1[j];
	Array1<U>::opt ubnxj=ubnx[j], ubnx1j=ubnx1[j];
	for(int k=0; k < Nz; k++) {
	  ub0j[k].phi=ub1j[k].phi;
	  ubnx1j[k].phi=ubnxj[k].phi;
	}
      }
    }
		
    if(FDdissipation || mu) {
      Laplacian(ub,w,PHI_ONLY);
      VorticityBoundaryConditions(w);
    }
	
    Real coeffx12=twelfth*hxinv;
    Real coeffy12=twelfth*hyinv;
    Real coeffvstar=nonlocal ? 0.0 : vstar;
    Real coeffmu=0.5*mu;
	
    for(int k=1; k <= Nzi; k++) coeffy2[k]=alpha*z[k];
			
    for(i=1; i <= Nxi; i++) {
      Array2<U> um=ub[i-1], ui=ub[i], up=ub[i+1];
      Array2<U> uP=(i+2 < Nx) ? ub[i+2] : (xperiodic ? ub[2] : ub[Nx-1]);
      Array2<U> uM=(i-2 >= 0) ? ub[i-2] : (xperiodic ? ub[Nx-3] : ub[0]);
		
      Array2<U> wm=w[i-1], wi=w[i], wp=w[i+1];
      Array2<U> wP=(i+2 < Nx) ? w[i+2] : (xperiodic ? w[2] : w[Nx-1]);
      Array2<U> wM=(i-2 >= 0) ? w[i-2] : (xperiodic ? w[Nx-3] : w[0]);
		
      Array2<Real> Hxi=Hx[i], Hyi=Hy[i];
      Array2<Real> F1xi=F1x[i], F1yi=F1y[i];
      Array2<Real> F2xi=F2x[i], F2yi=F2y[i];
      Array2<Real> F3xi=F3x[i], F3yi=F3y[i];
      for(int j=1; j <= Nyi; j++) {
	Array1<U>::opt umj=um[j], uim=ui[j-1], uij=ui[j], uip=ui[j+1];
	Array1<U>::opt upj=up[j], uPj=uP[j], uMj=uM[j];
	Array1<U>::opt uiP=(j+2 < Ny) ? ui[j+2] :
	  (yperiodic ? ui[2] : ui[Ny-1]);
	Array1<U>::opt uiM=(j-2 >= 0) ? ui[j-2] :
	  (yperiodic ? ui[Ny-3] : ui[0]);
		
	Array1<U>::opt wmj=wm[j], wim=wi[j-1], wip=wi[j+1], wpj=wp[j];
		
	Array1<U>::opt wPj=wP[j], wMj=wM[j];
	Array1<U>::opt wiP=(j+2 < Ny) ? wi[j+2] :
	  (yperiodic ? wi[2] : wi[Ny-1]);
	Array1<U>::opt wiM=(j-2 >= 0) ? wi[j-2] :
	  (yperiodic ? wi[Ny-3] : wi[0]);
			
	Array1<Real>::opt Hxij=Hxi[j], Hyij=Hyi[j];
	Array1<Real>::opt F1xij=F1xi[j], F1yij=F1yi[j];
	Array1<Real>::opt F2xij=F2xi[j], F2yij=F2yi[j];
	Array1<Real>::opt F3xij=F3xi[j], F3yij=F3yi[j];
			
	for(int k=1; k <= Nzi; k++) {
	  Real vx=-coeffy12*
	    (uiM[k].phi-8.0*(uim[k].phi-uip[k].phi)-uiP[k].phi);
	  Real Hyn=-uiM[k].n+2.0*(uim[k].n-uip[k].n)+uiP[k].n;
				
	  Real vy=coeffx12*
	    (uMj[k].phi-8.0*(umj[k].phi-upj[k].phi)-uPj[k].phi);
	  Real Hxn=-uMj[k].n+2.0*(umj[k].n-upj[k].n)+uPj[k].n;
				
	  Real avx=abs(vx)*coeffmu;
	  Real avy=abs(vy)*coeffmu;
					
	  vx *= nonlinear;
	  vy *= nonlinear;
				
	  Real alphaz=coeffy2[k];
	  Real phix=vy-alphaz*vx;
				
	  // Flux F1 = v * (-phi_x):
				
	  F1xij[k]=-vx*phix;
	  F1yij[k]=-vy*phix;
				
	  Real phiy=alphaz*vy-(1.0+alphaz*alphaz)*vx;
				
	  // Flux F2 = v * (-phi_y):
				
	  F2xij[k]=-vx*phiy;
	  F2yij[k]=-vy*phiy;
				
	  // Flux F3 = v * n:
				
	  F3xij[k]=vx*uij[k].n+avx*Hxn;
	  F3yij[k]=vy*uij[k].n+coeffvstar*uij[k].phi+avy*Hyn;
				
	  // Hyperviscosity: Equation 19 of
	  // Guzdar et al. Phys. Fluids B 5, p. 3712 (1993);
	  // implemented here as a flux correction.
				
	  Hxij[k]=avx*
	    (-wMj[k].phi+2.0*(wmj[k].phi-wpj[k].phi)+wPj[k].phi);
	  Hyij[k]=avy*
	    (-wiM[k].phi+2.0*(wim[k].phi-wip[k].phi)+wiP[k].phi);
	}
      }
    }
				
    Divergence(F3x,F3y,Sn);
    Divergence(F1x,F1y,F3x);
    Divergence(F2x,F2y,F3y);
	
    Array3<U> uoldlb=uoldl[b];
		
    for(i=1; i <= Nxi; i++) {
      Array2<Real> F3xi=F3x[i], F3yi=F3y[i];
      Array2<Real> Hxi=Hx[i], Hyi=Hy[i];
      Array2<U> Ui=uoldlb[i];
      for(int j=1; j <= Nyi; j++) {
	Array1<Real>::opt F3xij=F3xi[j], F3yij=F3yi[j];
	Array1<Real>::opt Hxij=Hxi[j], Hyij=Hyi[j];
	Array1<U>::opt Uij=Ui[j]; 
	for(int k=1; k <= Nzi; k++) {
	  Real n=Uij[k].n;
	  F3xij[k]=(F3xij[k]+Hxij[k])*n;
	  F3yij[k]=(F3yij[k]+Hyij[k])*n;
	}
      }
    }
	
    Divergence(F3x,F3y,Sphi);
	
    // Combine fluxes with Laplacian viscosity
	
    if(FDdissipation && (Dphi || Dn)) {
      if(Dphi || Dn) {
	for(int n=1; n <= nhyperviscosity; n++) {
	  if(n > 1) w=L2w;
	  Laplacian(w,L2w);
	  VorticityBoundaryConditions(L2w);
	}
      }
      for(i=1; i <= Nxi; i++) {
	Array2<U> Si=Sb[i], L2wi=L2w[i];
	Array2<Real> Sphii=Sphi[i], Sni=Sn[i];
	for(int j=1; j <= Nyi; j++) {
	  Array1<U>::opt Sij=Si[j]; 
	  Array1<U>::opt L2wij=L2wi[j];
	  Array1<Real>::opt Sphiij=Sphii[j];
	  Array1<Real>::opt Snij=Sni[j];
	  for(int k=1; k <= Nzi; k++) {
	    Sij[k].phi=-Sphiij[k]-Dphi*L2wij[k].phi;
	    Sij[k].n=-Snij[k]-Dn*L2wij[k].n;
	  }
	}
      } 
    } else {
      for(i=1; i <= Nxi; i++) {
	Array2<U> Si=Sb[i];
	Array2<Real> Sphii=Sphi[i], Sni=Sn[i];
	for(int j=1; j <= Nyi; j++) {
	  Array1<U>::opt Sij=Si[j];
	  Array1<Real>::opt Sphiij=Sphii[j], Snij=Sni[j];
	  for(int k=1; k <= Nzi; k++) {
	    Sij[k].phi=-Sphiij[k];
	    Sij[k].n=-Snij[k];
	  }
	}
      }
    }
  }
	
  if(pseudospectral) {
    FFT(source,S);
    for(int b=0; b < Nbox; b++) {
      if(!FDdissipation) {
	UComplex *Sc=(UComplex *) source;
	UComplex *Yc=(UComplex *) Y;
	int Nxib=Nxi*b;
	for(i=0; i < Nxi; i++) {
	  int i0=i-Nxi2;
	  Real kx=i0*kxmin;
	  Real kx2=kx*kx;
	  int Nyii=Nyi*(Nxib+i);
	  for(int j=0; j < Nyi; j++) {
	    int j0=j-Nyi2;
	    Real ky=j0*kymin;
	    Real kxy2=kx2+ky*ky;
	    Real kxy4=kxy2*kxy2;
	    Real khyper;
	    if(nhyperviscosity == 1) khyper=1.0;
	    else {
	      khyper=kxy2;
	      for(int n=2; n < nhyperviscosity; n++)
		khyper *= kxy2;
	    }
	    int Nij=Nzp*(Nyii+j);
	    Real nuphi=kxy4*Dphi*khyper;
	    Real nun=kxy2*Dn*khyper;
	    for(int k=((i0 == 0 && j0 == 0) ? 1 : 0); k < Nzp; k++)
	      {
		int index=Nij+k;
		Sc[index].phi -= nuphi*Yc[index].phi;
		Sc[index].n -= nun*Yc[index].n;
	      }
	  }
	}
      }
    }
  }
}
