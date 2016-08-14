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

extern unsigned spectrum;

extern int pH;
extern int pL;
 
class DNSBase {
protected:
  // Vocabulary:
  unsigned Nx;
  unsigned Ny;
  Real nuH, nuL;
  Real kH2, kL2;
  static const int xpad,ypad;
  
  enum Field {PAD,OMEGA,TRANSFER,EK};

  // derived variables:
  unsigned mx, my; // size of data arrays
  int imx; // (int) mx

  Real k0; // grid spacing factor
  Real k02; // k0^2
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
  
  Array2<Complex> f,g,h;
  uvector count;
  vector E; // Spectrum
  vector T; // Transfer
  Real Energy,Enstrophy,Palenstrophy;
  Array2<Real> k2inv;

public:
  void Initialize();
  void InitialConditions();
  void SetParameters();
  virtual void setcount();
  void FinalOutput();
  void OutFrame(int it);

  class FETL {
    DNSBase *b;
    const vector& E;
    const vector& T;
    Real k0;
    
  public: 
    FETL(DNSBase *b) : b(b), E(b->E), T(b->T), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2int);
      E[index] += Complex(w2/k,nu*w2);
      Complex& Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Sij,Tindex.im,k,i,j);
      Sij -= nu*wij;
    }
  };
  
  class FTL {
    DNSBase *b;
    const vector& T;
    Real k0;
    
  public: 
    FTL(DNSBase *b) : b(b), T(b->T), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Nu nu=b->nuk(k2int);
      Complex& Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Sij,Tindex.im,k,i,j);
      Sij -= nu*wij;
    }
  };
  
  class FET {
    DNSBase *b;
    const vector& E;
    const vector& T;
    Real k0;
    
  public: 
    FET(DNSBase *b) : b(b), E(b->E), T(b->T), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2int);
      E[index] += Complex(w2/k,nu*w2);
      Complex& Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Sij,Tindex.im,k,i,j);
    }
  };
  
  class FE {
    DNSBase *b;
    const vector& E;
    Real k0;
    
  public: 
    FE(DNSBase *b) : b(b), E(b->E), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2int);
      E[index] += Complex(w2/k,nu*w2);
    }
  };
  
  class FL {
    DNSBase *b;
    Real k0;
    
  public: 
    FL(DNSBase *b) : b(b), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      Complex wij=wi[j];
      Nu nu=b->nuk(k2int);
      double T;
      Forcing->Force(wij,Si[j],T,k,i,j);
      Si[j] -= nu*wij;
    }
  };
  
  class ForceStochastic {
    const vector& T;
    Real k0;
  public: 
    ForceStochastic(DNSBase *b) : T(b->T), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector&, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Forcing->ForceStochastic(wi[j],T[index].im,k);
    }
  };
  
  class ForceStochasticNO {
    Real k0;
  public: 
    ForceStochasticNO(DNSBase *b) : k0(b->k0) {}
    inline void operator()(const vector& wi, const vector&, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      double T;
      Forcing->ForceStochastic(wi[j],T,k);
    }
  };
  
  class InitializeValue {
    Real k0;
  public: 
    InitializeValue(DNSBase *b) : k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      wi[j]=InitialCondition->Value(k0*i,k0*j);
    }
  };
  
  class ForcingCount {
    DNSBase *b;
    Real k0;
  public: 
    ForcingCount(DNSBase *b) : b(b), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      if(Forcing->active(k)) ++b->fcount;
    }
  };
  
  class Invariants {
    DNSBase *b;
    Real &Energy,&Enstrophy,&Palenstrophy;
    Real k0;
  public: 
    Invariants(DNSBase *b) : b(b), Energy(b->Energy), Enstrophy(b->Enstrophy),
                             Palenstrophy(b->Palenstrophy), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      Real w2=abs2(wi[j]);
      Enstrophy += w2;
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      Real k2=k*k;
      Energy += w2/k2;
      Palenstrophy += k2*w2;
    }
  };
  
  class InitwS {
    DNSBase *b;
  public: 
    InitwS(DNSBase *b) : b(b) {}
    inline void operator()(vector& wi, vector& Si, int i) {
      Dimension(wi,b->w[i]);
      Dimension(Si,b->S[i]);
    }
  };
  
  class Initw {
    DNSBase *b;
  public: 
    Initw(DNSBase *b) : b(b) {}
    inline void operator()(vector& wi, vector& Si, int i) {
      Dimension(wi,b->w[i]);
    }
  };
  
  class InitNone {
  public: 
    InitNone(DNSBase *b) {}
    inline void operator()(vector& wi, vector& Si, int i) {
    }
  };
  
  class Count {
    DNSBase *b;
    const uvector& count;
  public: 
    Count(DNSBase *b) : b(b), count(b->count) {}
    inline void operator()(const vector& wi, const vector& Si, int i,
                           unsigned j) {
      unsigned k2int=i*i+j*j;
      Real kint=sqrt(k2int);
      unsigned index=(unsigned)(kint-0.5);
      ++count[index];
    }
  };
  
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);

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
    vector wi,Si;
    for(int i=-imx+1; i < imx; ++i) {
      init(wi,Si,i);
      for(unsigned j=i <= 0 ? 1 : 0; j < my; ++j)
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

  Real nuk(unsigned i2) {
    double k2=i2*k02;
    Real diss=0.0;
    if(k2 < kL2) diss += nuL*pow(k2,pL);
    if(k2 >= kH2) diss += nuH*pow(k2,pH);
    return diss;
  }

  virtual void ComputeInvariants(const array2<Complex>&, Real&, Real&, Real&);
  virtual Real getSpectrum(unsigned i) {
    double c=count[i];
    return c > 0 ? T[i].re*twopi/c : 0.0;
  }
  Real Dissipation(unsigned i) {return T[i].im;}
  Real Pi(unsigned i) {return T[i].re;}
  Real Eta(unsigned i) {return T[i].im;}
  Real kb(unsigned i) {return k0*(i+0.5);}
  Real kc(unsigned i) {return k0*(i+1);}
};

extern InitialConditionBase *InitialCondition;

#endif
