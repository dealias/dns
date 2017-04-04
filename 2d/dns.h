#ifndef __dnsbase_h__
#define __dnsbase_h__ 1

#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "fftw++.h"
#include "convolution.h"
#include "explicit.h"
#include "Forcing.h"
#include "InitialCondition.h"
#include "Conservative.h"
#include "Exponential.h"
#include <sys/stat.h> // On Sun computers this must come after xstream.h

using namespace Array;
using namespace fftwpp;
using std::ostringstream;

typedef Array1<Var>::opt Vector;
typedef Array1<Real>::opt rVector;

extern unsigned spectrum;

extern int pH;
extern int pL;
 
class DNSBase {
protected:
  // Vocabulary:
  unsigned Nx;
  unsigned Ny;
  Real nuH,nuL;
  Real kH2,kL2;
  
  // TRANSFER, INJECTION, INJECTION2, and DISSIPATION must be contiguous!
  enum Field {PAD,OMEGA,TRANSFER,INJECTION,INJECTION2,DISSIPATION,EK};

  // derived variables:
  int mx,my; // size of data arrays

  Array2<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field;

  int tcount;
  unsigned fcount;

  unsigned nmode;
  unsigned nshells;  // Number of spectral shells

  Array2<Complex> f0,f1;
  Array2<Complex> S;
  array2<Complex> buffer;
  Complex *F[2];
  Complex *block;
  ImplicitHConvolution2 *Convolution;
  crfft2d *Backward;
  
  ifstream ftin;
  oxstream fwk,fw,fekvk,ftransfer;
  ofstream ft,fevt;
  
  uvector count;
  
  vector E; // Energy spectrum
  vector T; // Enstrophy and energy transfers (Pi_Z and Pi_E)
  vector I; // Enstrophy and energy injection rates (eps and eta)
  vector J; // Palenstrophy injection rate (zeta)
  vector D; // Enstrophy and energy dissipation rates
  
  Real Energy,Enstrophy,Palinstrophy;
  Array2<Real> k2inv;

public:
  void Initialize() {
    fevt << "# t\tE\tZ\tP" << endl;
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
  
    k2inv.Allocate(Nx,my,-mx+1,0);
    for(int i=-mx+1; i < mx; ++i) {
      int i2=i*i;
      rVector k2invi=k2inv[i];
      for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
        k2invi[j]=1.0/(i2+j*j);
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
  
  void FinalOutput() {
    Real E,Z,P;
    ComputeInvariants(w,E,Z,P);
    cout << endl;
    cout << "Energy = " << E << newl;
    cout << "Enstrophy = " << Z << newl;
    cout << "Palinstrophy = " << P << newl;
  }
  
  void OutFrame(int it) {
    for(int i=-mx+1; i < mx; ++i)
      for(int j=0; j < my; ++j)
        f1[i][j]=w(i,j);
    
    fftwpp::HermitianSymmetrizeX(mx,my,mx-1,f1);

    f1[-mx][0]=0; // Zero remaining Nyquist mode.
    
    Backward->fft0(f1);
     
    fw << 1 << 2*my << Nx+1;
    for(int j=2*my-1; j >= 0; j--)
      for(unsigned i=0; i <= Nx; i++)
        fw << (float) wr(i,j);
    fw.flush();
  }

  class FETL {
    DNSBase *b;
    const vector& E,T,I,J,D;
    
  public: 
    FETL(DNSBase *b) : b(b), E(b->E), T(b->T), I(b->I), J(b->J), D(b->D) {}
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
      E[index] += kinv*w2;
      T[index] += Complex(transfer,kinv2*transfer);
      I[index] += Complex(eta,kinv2*eta);
      J[index] += k2*eta;
      D[index] += Complex(nuk2Z,kinv2*nuk2Z);
      Sij -= nuk2*wij;
    }
  };
  
  class FTL {
    DNSBase *b;
    const vector& T,I,J,D;
    
  public: 
    FTL(DNSBase *b) : b(b), T(b->T), I(b->I), J(b->J), D(b->D) {}
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
      T[index] += Complex(transfer,kinv2*transfer);
      I[index] += Complex(eta,kinv2*eta);
      J[index] += k2*eta;
      D[index] += Complex(nuk2Z,kinv2*nuk2Z);
      Sij -= nuk2*wij;
    }
  };
  
  class FET {
    DNSBase *b;
    const vector& E,T,I,J,D;
    
  public: 
    FET(DNSBase *b) : b(b), E(b->E), T(b->T), I(b->I), J(b->J), D(b->D) {}
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
      E[index] += kinv*w2;
      T[index] += Complex(transfer,kinv2*transfer);
      I[index] += Complex(eta,kinv2*eta);
      J[index] += k2*eta;
      D[index] += Complex(nuk2Z,kinv2*nuk2Z);
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
    const vector& I,J;
  public: 
    ForceStochastic(DNSBase *b) : I(b->I), J(b->J) {}
    inline void operator()(const Vector& wi, const Vector&, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      double eta=Forcing->ForceStochastic(wi[j],i,j);
      Real kinv2=1.0/k2;
      I[index] += Complex(eta,kinv2*eta);
      J[index] += k2*eta;
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
    Real &Energy,&Enstrophy,&Palinstrophy;
  public: 
    Invariants(DNSBase *b) : b(b), Energy(b->Energy), Enstrophy(b->Enstrophy),
                             Palinstrophy(b->Palinstrophy) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      Real w2=abs2(wi[j]);
      Enstrophy += w2;
      unsigned k2=i*i+j*j;
      Energy += w2/k2;
      Palinstrophy += k2*w2;
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
    f0.Dimension(Nx+1,my,-mx,0);
  
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
    Convolution->convolve(F,multadvection2);
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
      Init(T,Src[TRANSFER]);
      Init(I,Src[INJECTION]);
      Init(J,Src[INJECTION2]);
      Init(D,Src[DISSIPATION]);
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
      Init(E,Src[EK]);
      Init(T,Src[TRANSFER]);
      Init(I,Src[INJECTION]);
      Init(J,Src[INJECTION2]);
      Init(D,Src[DISSIPATION]);
      Compute(FET(this),Src,Y);
    }
  }

  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(E,Src[EK]);
      Init(T,Src[TRANSFER]);
      Init(I,Src[INJECTION]);
      Init(J,Src[INJECTION2]);
      Init(D,Src[DISSIPATION]);
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
      Set(I,Y[INJECTION]);
      Set(J,Y[INJECTION2]);
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
                                 Real& P) {
    Energy=Enstrophy=Palinstrophy=0.0;
  
    Loop(Initw(this),Invariants(this));
  
    E=Energy;
    Z=Enstrophy;
    P=Palinstrophy;
  }
  
  virtual Real getSpectrum(unsigned i) {
    double c=count[i];
    return c > 0 ? E[i].re*twopi/c : 0.0;
  }
  Real TZ(unsigned i) {return T[i].re;}
  Real TE(unsigned i) {return T[i].im;}
  
  Real Eta(unsigned i) {return I[i].re;}
  Real Eps(unsigned i) {return I[i].im;}
  
  Real Zeta(unsigned i) {return J[i].re;}
  
  Real DZ(unsigned i) {return D[i].re;}
  Real DE(unsigned i) {return D[i].im;}
  
  Real kb(unsigned i) {return i+0.5;}
  Real kc(unsigned i) {return i+1;}
};

extern InitialConditionBase *InitialCondition;

#endif
