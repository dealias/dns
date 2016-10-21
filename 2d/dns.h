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
  
  enum Field {PAD,OMEGA,TRANSFER,EK};

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
  ofstream ft,fevt,fforce;
  
   uvector count;
  vector E; // Spectrum
  vector T; // Transfer
  Real Energy,Enstrophy,Palenstrophy;
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
    cout << "Palenstrophy = " << P << newl;
  }
  
  void OutFrame(int it) {
    for(int i=-mx+1; i < mx; ++i)
      for(int j=0; j < my; ++j)
        f1[i][j]=w(i,j);

    for(int j=0; j < my; ++j) // Fill in Nyquist mode assuming continuity.
      f1[-mx][j]=f1[-mx+1][j];
    
    Backward->fft0(f1);
     
    fw << 1 << 2*my << Nx+1;
    for(int j=2*my-1; j >= 0; j--)
      for(unsigned i=0; i <= Nx; i++)
        fw << (float) wr(i,j);
    fw.flush();
  }

  class FETL {
    DNSBase *b;
    const vector& E;
    const vector& T;
    
  public: 
    FETL(DNSBase *b) : b(b), E(b->E), T(b->T) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2);
      E[index] += Complex(w2/k,nu*w2);
      Complex& Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Sij,Tindex.im,i,j);
      Sij -= nu*wij;
    }
  };
  
  class FTL {
    DNSBase *b;
    const vector& T;
    
  public: 
    FTL(DNSBase *b) : b(b), T(b->T) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Nu nu=b->nuk(k2);
      Complex& Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Sij,Tindex.im,i,j);
      Sij -= nu*wij;
    }
  };
  
  class FET {
    DNSBase *b;
    const vector& E;
    const vector& T;
    
  public: 
    FET(DNSBase *b) : b(b), E(b->E), T(b->T) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2);
      E[index] += Complex(w2/k,nu*w2);
      Complex& Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Sij,Tindex.im,i,j);
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
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2);
      E[index] += Complex(w2/k,nu*w2);
    }
  };
  
  class FL {
    DNSBase *b;
    
  public: 
    FL(DNSBase *b) : b(b) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Complex wij=wi[j];
      Nu nu=b->nuk(k2);
      double T;
      Forcing->Force(wij,Si[j],T,i,j);
      Si[j] -= nu*wij;
    }
  };
  
  class ForceStochastic {
    const vector& T;
  public: 
    ForceStochastic(DNSBase *b) : T(b->T) {}
    inline void operator()(const Vector& wi, const Vector&, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Forcing->ForceStochastic(wi[j],T[index].im,i,j);
    }
  };
  
  class ForceStochasticNO {
  public: 
    ForceStochasticNO(DNSBase *b) {}
    inline void operator()(const Vector& wi, const Vector&, int i, int j) {
      double T;
      Forcing->ForceStochastic(wi[j],T,i,j);
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

        Complex w=0.0;
        double T=0.0;
        Complex S=0.0;
        Forcing->Force(w,S,T,i,j);
        if(S != 0.0)
          b->fforce << i << " " << j << " " << S << endl;
      }
    }
  };
  
  class Invariants {
    DNSBase *b;
    Real &Energy,&Enstrophy,&Palenstrophy;
  public: 
    Invariants(DNSBase *b) : b(b), Energy(b->Energy), Enstrophy(b->Enstrophy),
                             Palenstrophy(b->Palenstrophy) {}
    inline void operator()(const Vector& wi, const Vector& Si, int i, int j) {
      Real w2=abs2(wi[j]);
      Enstrophy += w2;
      unsigned k2=i*i+j*j;
      Energy += w2/k2;
      Palenstrophy += k2*w2;
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
      Init(T,Src[TRANSFER]);
      Init(E,Src[EK]);
      Compute(FET(this),Src,Y);
    }
  }

  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(T,Src[TRANSFER]);
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
      Set(T,Y[TRANSFER]);
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
    Energy=Enstrophy=Palenstrophy=0.0;
  
    Loop(Initw(this),Invariants(this));
  
    E=Energy;
    Z=Enstrophy;
    P=Palenstrophy;
  }
  
  virtual Real getSpectrum(unsigned i) {
    double c=count[i];
    return c > 0 ? T[i].re*twopi/c : 0.0;
  }
  Real Dissipation(unsigned i) {return T[i].im;}
  Real Pi(unsigned i) {return T[i].re;}
  Real Eta(unsigned i) {return T[i].im;}
  Real kb(unsigned i) {return i+0.5;}
  Real kc(unsigned i) {return i+1;}
};

extern InitialConditionBase *InitialCondition;

#endif
