#ifndef __dnsbase_h__
#define __dnsbase_h__ 1

#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "fftw++.h"
//#include "convolution.h"
#include "tests/convolve.h"
#include "Forcing.h"
#include "InitialCondition.h"
#include "Conservative.h"
#include "Exponential.h"
#include <sys/stat.h> // On Sun computers this must come after xstream.h
#include "tests/explicit.h"

using namespace Array;
using namespace fftwpp;

using std::ostringstream;
using fftwpp::twopi;

typedef Array1<Var>::opt Vector;
typedef Array1<Real>::opt rVector;

extern unsigned spectrum;

extern int pH;
extern int pL;

// This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
// requires only 4 FFTs per stage.
void multadvection2(Complex **F, unsigned int offset, unsigned int n,
                    /*
                    const unsigned int indexsize,
                    const unsigned int *index,
                    unsigned int r,
                    */
                    unsigned int threads)
{
  double* F0=(double *) (F[0]+offset);
  double* F1=(double *) (F[1]+offset);

  for(unsigned int j=0; j < n; ++j) {
    double u=F0[j];
    double v=F1[j];
    F0[j]=v*v-u*u;
    F1[j]=u*v;
  }
}

double Triplet0, Triplet, Triplet2, Norm1, Norm2;

// A=8, B=0 TODO: Reduce to A=7, B=0 using incompressibility
void multTriplet(Complex **F, unsigned int offset, unsigned int n,
                 unsigned int threads)
{
  double* F0=(double *) (F[0]+offset);
  double* F1=(double *) (F[1]+offset);
  double* F2=(double *) (F[2]+offset);
  double* F3=(double *) (F[3]+offset);
  double* F4=(double *) (F[4]+offset);
  double* F5=(double *) (F[5]+offset);
  double* F6=(double *) (F[6]+offset);
  double* F7=(double *) (F[7]+offset);

  for(unsigned int j=0; j < n; ++j) {
    double u=F0[j];
    double v=F1[j];
    double ux=F2[j];
    double uy=F3[j];
    double vx=F4[j];
    double vy=F5[j];
    double A2u=F6[j];
    double A2v=F7[j];
    double bx=u*ux+v*uy;
    double by=u*vx+v*vy;

    double w=vx-uy;
    double E=u*u+v*v;
    double Z=w*w;
    Triplet0 += E;
    Triplet += Z;
    Triplet2 += A2u*u+A2v*v;
    F0[j]=E;
    F1[j]=Z;

//    Triplet += u*u+v*v;
//    F0[j]=u*u+v*v;

//    Triplet += bx*A2u+by*A2v;
//    Norm1 += bx*bx+by*by;
//    Norm2 += A2u*A2u+A2v*A2v;
  }
}

class DNSBase {
protected:
  // Vocabulary:
  unsigned Nx;
  unsigned Ny;
  unsigned nx,ny0;
  Real nuH,nuL;
  Real kH2,kL2;

  // Contiguous: TRANSFERE,TRANSFERZ,EPS,ETA,ZETA,DISSIPATIONE,DISSIPATIONZ
  //
  enum Field {PAD,OMEGA,TRANSFERE,TRANSFERZ,EPS,ETA,ZETA,DISSIPATIONE,
              DISSIPATIONZ,EK};

  int mx,my; // size of data arrays

  Array2<Complex> w; // Vorticity field
  Array2<Complex> w0; // Temporary expanded vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field

  Array2<Complex> u,v,ux,uy,vx,vy,A2u,A2v;
  array2<Real> ur,vr,uxr,uyr,vxr,vyr,A2ur,A2vr; // Inverse Fourier transforms

  int tcount;
  unsigned fcount;

  unsigned nmode;
  unsigned nshells;  // Number of spectral shells

  Array2<Complex> f0,f1;
  Array2<Complex> S;
  array2<Complex> buffer;
  Complex *F[8];
  Complex *block;
//  ImplicitHConvolution2 *Convolution;
  ConvolutionHermitian2 *Convolution;
  ConvolutionHermitian2 *Convolution2;
  crfft2d *Backward2;

  ifstream ftin;
  oxstream fek,fw,fekvk,ftransfer;
  ofstream ft,fevt;

  uvector count;

  vector TE,TZ; // Energy and enstrophy transfers
  vector Eps,Eta,Zeta; // Energy, enstrophy, and palenstrophy injection rates
  vector DE,DZ; // Energy and enstrophy dissipation rates
  vector E; // Energy spectrum
  Real Energy,Enstrophy,Palinstrophy,Palinstrophy2;
  Array2<Real> k2inv;

public:
  void Initialize() {
    fevt << "# t\tE\tZ\tP\tP2" << endl;
  }

  void InitialConditions() {
    w[0][0]=0.0; // Enforce no mean flow
    Loop(Initw(this),InitializeValue(this));
    fftwpp::HermitianSymmetrizeX(mx,my,mx-1,w);
  }

  void SetParameters() {
    setcount();
    fcount=0;

    Forcing->Init();

    Loop(InitNone(this),ForcingCount(this));

    fcount *= 2; // Account for Hermitian conjugate modes.

    Forcing->Init(fcount);

    double scale=sqrt(Convolution->scale);

    k2inv.Allocate(Nx,my,-mx+1,0);
    for(int i=-mx+1; i < mx; ++i) {
      int i2=i*i;
      rVector k2invi=k2inv[i];
      for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
        k2invi[j]=scale/(i2+j*j);
      }
    }
  }

  virtual void setcount() {
#pragma omp parallel for num_threads(threads)
    for(unsigned i=0; i < nshells; i++)
      count[i]=0;

    if(spectrum)
      Loop(InitNone(this),Count(this));
  }

  virtual void OutEnergies() {
    fftwpp::HermitianSymmetrizeX(mx,my,mx-1,w);
    fek << 2*mx-1 << my;
    for(int i=-mx+1; i < mx; ++i) {
      const Vector& wi=w[i];
      for(int j=0; j < my; ++j) {
        Real k2=i*i+j*j;
        Real k2inv=k2 > 0.0 ? 1.0/k2 : 0.0;
        fek << 0.5*abs2(wi[j])*k2inv;
      }
    }
  }

  void FinalOutput() {
    Real E,Z,P,P2;
    ComputeInvariants(w,E,Z,P,P2);
    cout << endl;
    cout << "Energy = " << E << newl;
    cout << "Enstrophy = " << Z << newl;
    cout << "Palinstrophy = " << P << newl;
    cout << "Palinstrophy2 = " << P2 << newl;
  }

  void OutFrame(int it) {
    w0=0.0;
    u=0.0;
    v=0.0;
    ux=0.0;
    uy=0.0;
    vx=0.0;
    vy=0.0;
    A2u=0.0;
    A2v=0.0;

    for(int i=-mx+1; i < mx; ++i)
      for(int j=0; j < my; ++j)
        w0[i][j]=w(i,j);

    for(int i=-mx+1; i < mx; ++i) {
      for(int j=0; j < my; ++j) {
        Complex wij=w(i,j);
        Real k4=i*i+j*j;
        k4 *= k4;
        Real k2=i*i+j*j;
        Real k2inv=k2 > 0.0 ? 1.0/k2 : 0.0;
        Complex psi=wij*k2inv;
        Complex uij=I*j*psi;
        Complex vij=-I*i*psi;
        u[i][j]=uij;
        v[i][j]=vij;
        ux[i][j]=I*i*uij;
        uy[i][j]=I*j*uij;
        vx[i][j]=I*i*vij;
        vy[i][j]=I*j*vij;
        A2u[i][j]=k4*uij;
        A2v[i][j]=k4*vij;
      }
    }

//    unsigned nyp=ny0/2+1;

    fftwpp::HermitianSymmetrizeX(mx,my,mx,w0);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,u);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,v);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,ux);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,uy);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,vx);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,vy);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,A2u);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,A2v);

//    Pad2->pad(w0);
    /*
   // Zero Nyquist modes.
    for(int j=0; j < my; ++j) {
      w0(j)=0.0;
      u(j)=0.0;
      v(j)=0.0;
      ux(j)=0.0;
      uy(j)=0.0;
      vx(j)=0.0;
      vy(j)=0.0;
      A2u(j)=0.0;
      A2v(j)=0.0;
    }
    */

    F[0]=u;
    F[1]=v;
    F[2]=ux;
    F[3]=uy;
    F[4]=vx;
    F[5]=vy;
    F[6]=A2u;
    F[7]=A2v;

    Triplet0=0.0;
    Triplet=0.0;
    Triplet2=0.0;
    Norm1=0.0;
    Norm2=0.0;

    Convolution2->convolveRaw(F,multTriplet);

    f0.Set(F[0]);
    double E=f0[0][0].re;

    f0.Set(F[1]);
    double Z=f0[0][0].re;

    /*
    Backward2->fft0(w0);
    Backward2->fft0(u);
    Backward2->fft0(v);
    Backward2->fft0(ux);
    Backward2->fft0(uy);
    Backward2->fft0(vx);
    Backward2->fft0(vy);
    Backward2->fft0(A2u);
    Backward2->fft0(A2v);
    */

/*
    fw << 1 << ny0 << nx;
    for(int j=ny0-1; j >= 0; j--)
      for(unsigned i=0; i < nx; i++)
        fw << (float) wr(i,j);
    fw.flush();
*/

    /*
    cout << "H=" << H << endl;
    double h=H/sqrt(norm1*norm2);
    cout << h << endl;
    if(h > 1.0) h=1.0;
    if(h < -1.0) h=-1.0;
    cout << "Angle=" << acos(h)*180.0/PI << endl;
    */

    F[1]=f1;

    /*
    for(unsigned i=0; i < nx; i++) {
      for(unsigned j=0; j < ny0; j++) {
        Real bx=ur[i][j]*uxr[i][j]+vr[i][j]*uyr[i][j];
        Real by=ur[i][j]*vxr[i][j]+vr[i][j]*vyr[i][j];
        Real ax=A2ur[i][j];
        Real ay=A2vr[i][j];
        sum += bx*ax+by*ay;
        Z += wr[i][j]*wr[i][j];
        E += ur(i,j)*ur(i,j)+vr(i,j)*vr(i,j);
        norm1 += bx*bx+by*by;
        norm2 += ax*ax+ay*ay;
      }
    }
    */

    double scale=Convolution2->scale;
    cout << "energy = " << 0.5*E*scale << endl;
    cout << "energy = " << 0.5*Triplet0*scale << endl;
    cout << "enstrophy = " << 0.5*Z*scale << endl;
    cout << "enstrophy = " << 0.5*Triplet*scale << endl;
    cout << "palinstrophy = " << 0.5*Triplet2*scale << endl;
//    cout << "enstrophy = " << 0.5*Triplet/(nx*ny0) << endl;
//    cout << "Inner product=" << sum/((Nx+1)*(2*my-1)) << endl;
//    cout << "Angle=" << acos(sum/sqrt(norm1*norm2))*180.0/PI << endl;

/*
        f1[i][j]=w(i,j);

    fftwpp::HermitianSymmetrizeX(mx,my,mx,f1);

// Zero Nyquist modes.
    for(int j=0; j < my; ++j)
      f1(j)=0.0;

    Backward->fft0(f1);

    fw << 1 << 2*my << Nx+1;
    for(int j=2*my-1; j >= 0; j--)
      for(unsigned i=0; i <= Nx; i++)
        fw << (float) wr(i,j);
    fw.flush();
*/

    /*
// Zero Nyquist modes.
    for(int j=0; j < my; ++j)
      f1(j)=0.0;
    */
  }

  class FETL {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ,E;

  public:
    FETL(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                       Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                       DE(b->DE), DZ(b->DZ), E(b->E) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv=1.0/k;
      Real kinv2=kinv*kinv;
      Nu nuk2=b->nuk(k2);
      Real nuk2Z=nuk2*w2;
      TE[index] += kinv2*transfer;
      TZ[index] += transfer;
      Eps[index] += kinv2*eta;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
      DE[index] += kinv2*nuk2Z;
      DZ[index] += nuk2Z;
      E[index] += kinv*w2;
      Sij -= nuk2*wij;
    }
  };

  class FTL {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ;

  public:
    FTL(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                      Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                      DE(b->DE), DZ(b->DZ) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv2=1.0/k2;
      Nu nuk2=b->nuk(k2);
      Real nuk2Z=nuk2*abs2(wij);
      TE[index] += kinv2*transfer;
      TZ[index] += transfer;
      Eps[index] += kinv2*eta;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
      DE[index] += kinv2*nuk2Z;
      DZ[index] += nuk2Z;
      Sij -= nuk2*wij;
    }
  };

  class FET {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ,E;

  public:
    FET(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                      Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                      DE(b->DE), DZ(b->DZ), E(b->E) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv=1.0/k;
      Real kinv2=kinv*kinv;
      Real nuk2Z=b->nuk(k2)*w2;
      TE[index] += kinv2*transfer;
      TZ[index] += transfer;
      Eps[index] += kinv2*eta;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
      DE[index] += kinv2*nuk2Z;
      DZ[index] += nuk2Z;
      E[index] += kinv*w2;
    }
  };

  class FE {
    DNSBase *b;
    const vector& E;

  public:
    FE(DNSBase *b) : b(b), E(b->E) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      E[index] += abs2(wi[j])/k;
    }
  };

  class FL {
    DNSBase *b;

  public:
    FL(DNSBase *b) : b(b) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Complex wij=wi[j];
      Forcing->Force(wij,Si[j],i,j);
      Si[j] -= b->nuk(k2)*wij;
    }
  };

  class ForceStochastic {
    const vector& Eps,Eta,Zeta;
  public:
    ForceStochastic(DNSBase *b) : Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta) {}
    inline void operator()(const Vector& wi, const Vector&, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      double eta=Forcing->ForceStochastic(wi[j],i,j);
      Eps[index] += eta/k2;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
    }
  };

  class ForceStochasticNO {
  public:
    ForceStochasticNO(DNSBase *b) {}
    inline void operator()(const Vector& wi, const Vector&, int i, int j) {
      Forcing->ForceStochastic(wi[j],i,j);
    }
  };

  class InitializeValue {
  public:
    InitializeValue(DNSBase *b) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      wi[j]=InitialCondition->Value(i,j);
    }
  };

  class ForcingCount {
    DNSBase *b;
  public:
    ForcingCount(DNSBase *b) : b(b) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      if(Forcing->active(i,j)) {
        ++b->fcount;
      }
    }
  };

  class Invariants {
    DNSBase *b;
    Real &Energy,&Enstrophy,&Palinstrophy,&Palinstrophy2;
  public:
    Invariants(DNSBase *b) : b(b), Energy(b->Energy), Enstrophy(b->Enstrophy),
                             Palinstrophy(b->Palinstrophy),
                             Palinstrophy2(b->Palinstrophy2) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      Real w2=abs2(wi[j]);
      Enstrophy += w2;
      unsigned k2=i*i+j*j;
      Energy += w2/k2;
      Palinstrophy += k2*w2;
      Palinstrophy2 += k2*k2*w2;
    }
  };

  class InitwS {
    DNSBase *b;
  public:
    InitwS(DNSBase *b) : b(b) {}
    inline void operator()(Vector& wi, Vector& Si, int i) {
      Dimension(wi,b->w[i]);
      Dimension(Si,b->S[i]);
    }
  };

  class Initw {
    DNSBase *b;
  public:
    Initw(DNSBase *b) : b(b) {}
    inline void operator()(Vector& wi, Vector& Si, int i) {
      Dimension(wi,b->w[i]);
    }
  };

  class InitNone {
  public:
    InitNone(DNSBase *b) {}
    inline void operator()(Vector& wi, Vector& Si, int i) {
    }
  };

  class Count {
    DNSBase *b;
    const uvector& count;
  public:
    Count(DNSBase *b) : b(b), count(b->count) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      ++count[index];
    }
  };

  void NonLinearSource(const vector2& Src, const vector2& Y, double t) {
//    f0.Dimension(Nx+1,my,-mx,0);

    w.Set(Y[OMEGA]);
    f0.Set(Src[PAD]);

    f0[0][0]=0.0;
    f1[0][0]=0.0;

    // This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
    // requires only 4 FFTs per stage.
#pragma omp parallel for num_threads(threads)
    for(int i=-mx+1; i < mx; ++i) {
      Vector wi=w[i];
      Vector f0i=f0[i];
      Vector f1i=f1[i];
      rVector k2invi=k2inv[i];
      for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
        Complex wij=wi[j];
        Real k2invij=k2invi[j];
        Real jk2inv=j*k2invij;
        Real ik2inv=i*k2invij;
        f0i[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv); // u
        f1i[j]=Complex(wij.im*ik2inv,-wij.re*ik2inv); // v
      }
    }

    F[0]=f0;

    fftwpp::HermitianSymmetrizeX(mx,my,mx,f0);
    fftwpp::HermitianSymmetrizeX(mx,my,mx,f1);

    Convolution->convolveRaw(F,multadvection2);
    f0[0][0]=0.0;

    for(int i=-mx+1; i < mx; ++i) {
      Real i2=i*i;
      Vector f0i=f0[i];
      Vector f1i=f1[i];
      for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
        f0i[j]=i*j*f0i[j]+(i2-j*j)*f1i[j];
      }
    }
    fftwpp::HermitianSymmetrizeX(mx,my,mx,f0);

#if 0
    Real sum=0.0;
    for(int i=-mx+1; i < mx; ++i) {
      Vector wi=w[i];
      for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
        Complex wij=wi[j];
//      sum += (f0[i][j]*conj(wij)).re;
        sum += (f0[i][j]*conj(wij)/(i*i+j*j)).re;
      }
    }
    cout << sum << endl;
    cout << endl;
#endif
  }

  void Init(vector& T, const vector& Src) {
    Set(T,Src);
#pragma omp parallel for num_threads(threads)
    for(unsigned K=0; K < nshells; K++)
      T[K]=0.0;
  }

  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(TE,Src[TRANSFERE]);
      Init(TZ,Src[TRANSFERZ]);
      Init(Eps,Src[EPS]);
      Init(Eta,Src[ETA]);
      Init(Zeta,Src[ZETA]);
      Init(DE,Src[DISSIPATIONE]);
      Init(DZ,Src[DISSIPATIONZ]);
      Compute(FTL(this),Src,Y);
    }
    else
      Compute(FL(this),Src,Y);
  }

  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum) {
      Init(E,Src[EK]);
      Compute(FE(this),Src,Y);
    }
  }

  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(TE,Src[TRANSFERE]);
      Init(TZ,Src[TRANSFERZ]);
      Init(Eps,Src[EPS]);
      Init(Eta,Src[ETA]);
      Init(Zeta,Src[ZETA]);
      Init(DE,Src[DISSIPATIONE]);
      Init(DZ,Src[DISSIPATIONZ]);
      Init(E,Src[EK]);
      Compute(FET(this),Src,Y);
    }
  }

  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(TE,Src[TRANSFERE]);
      Init(TZ,Src[TRANSFERZ]);
      Init(Eps,Src[EPS]);
      Init(Eta,Src[ETA]);
      Init(Zeta,Src[ZETA]);
      Init(DE,Src[DISSIPATIONE]);
      Init(DZ,Src[DISSIPATIONZ]);
      Init(E,Src[EK]);
      Compute(FETL(this),Src,Y);
    } else
      Compute(FL(this),Src,Y);
  }

  template<class S, class T>
  void Loop(S init, T fcn)
  {
    Vector wi,Si;
    for(int i=-mx+1; i < mx; ++i) {
      init(wi,Si,i);
      for(int j=i <= 0 ? 1 : 0; j < my; ++j)
        fcn(wi,Si,i,j);
    }
  }

  template<class T>
  void Compute(T fcn, const vector2& Src, const vector2& Y)
  {
    S.Set(Src[OMEGA]);
    w.Set(Y[OMEGA]);

    Loop(InitwS(this),fcn);
  }

  void Stochastic(const vector2&Y, double, double dt)
  {
    if(!Forcing->Stochastic(dt)) return;
    w.Set(Y[OMEGA]);

    if(spectrum == 0) {
      Loop(Initw(this),ForceStochasticNO(this));
    } else {
      Set(Eps,Y[EPS]);
      Set(Eta,Y[ETA]);
      Set(Zeta,Y[ZETA]);
      Loop(Initw(this),ForceStochastic(this));
    }
  }

  Nu LinearCoeff(unsigned k) {
    unsigned i=k/my;
    unsigned j=k-my*i;
    return nuk(i*i+j*j);
  }

  Real nuk(double k2) {
    Real diss=0.0;
    if(k2 < kL2) diss += nuL*pow(k2,pL);
    if(k2 >= kH2) diss += nuH*pow(k2,pH);
    return diss;
  }

  virtual void ComputeInvariants(const array2<Complex> &w, Real& E, Real& Z,
                                 Real& P, Real& P2) {
    Energy=Enstrophy=Palinstrophy=Palinstrophy2=0.0;

    Loop(Initw(this),Invariants(this));

    E=Energy;
    Z=Enstrophy;
    P=Palinstrophy;
    P2=Palinstrophy2;
  }

  virtual Real getSpectrum(unsigned i) {
    double c=count[i];
    return c > 0 ? E[i].re*twopi/c : 0.0;
  }
  Real TE_(unsigned i) {return TE[i].re;}
  Real TZ_(unsigned i) {return TZ[i].re;}
  Real Eps_(unsigned i) {return Eps[i].re;}
  Real Eta_(unsigned i) {return Eta[i].re;}
  Real Zeta_(unsigned i) {return Zeta[i].re;}
  Real DE_(unsigned i) {return DE[i].re;}
  Real DZ_(unsigned i) {return DZ[i].re;}

  Real kb(unsigned i) {return i+0.5;}
  Real kc(unsigned i) {return i+1;}
};

extern InitialConditionBase *InitialCondition;

#endif
