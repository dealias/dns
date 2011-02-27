#ifndef __dnsbase_h__
#define __dnsbase_h__ 1

#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "fftw++.h"
#include "convolution.h"
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
  static const int xpad,ypad;
  
  enum Field {OMEGA,TRANSFER,EK};

  // derived variables:
  unsigned mx, my; // size of data arrays
  unsigned origin; // linear index of Fourier origin.
  unsigned xorigin; // horizontal index of Fourier origin.

  Real k0; // grid spacing factor
  Real k02; // k0^2
  array2<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field;

  int tcount;

protected:  
  unsigned nmode;
  unsigned nshells;  // Number of spectral shells

  array2<Complex> f0,f1,g0,g1;
  array2<Complex> buffer;
  Complex *F[2];
  Complex *G[2];
  Complex *block;
  ImplicitHConvolution2 *Convolution;
  ExplicitHConvolution2 *Padded;

  ifstream ftin;
  oxstream fwk,fw,fekvk,ftransfer;
  ofstream ft,fevt;
  
  Array2<Complex> f,g,h;
  array1<unsigned>::opt count;
  vector E; // Spectrum
  vector T; // Transfer

public:
  DNSBase() {}
  virtual ~DNSBase() {}

  unsigned getNx() {return Nx;}
  unsigned getmx() {return mx;}
  unsigned getmy() {return my;}
  Real getk0() {return k0;}
  Real getk02() {return k02;}
  unsigned getxorigin() {return xorigin;}

  void InitialConditions();
  void Initialize();
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
    inline void operator()(const vector& wi, const vector& Si, unsigned I2,
                           unsigned j) {
      unsigned k2int=I2+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2int);
      E[index] += Complex(w2/k,nu*w2);
      Complex Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Tindex.im,k);
      Si[j]=Sij-nu*wij;
    }
  };
  
  class FTL {
    DNSBase *b;
    const vector& T;
    Real k0;
    
  public: 
    FTL(DNSBase *b) : b(b), T(b->T), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, unsigned I2,
                           unsigned j) {
      unsigned k2int=I2+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Nu nu=b->nuk(k2int);
      Complex Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Tindex.im,k);
      Si[j]=Sij-nu*wij;
    }
  };
  
  class FET {
    DNSBase *b;
    const vector& E;
    const vector& T;
    Real k0;
    
  public: 
    FET(DNSBase *b) : b(b), E(b->E), T(b->T), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, unsigned I2,
                           unsigned j) {
      unsigned k2int=I2+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      unsigned index=(unsigned)(kint-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Nu nu=b->nuk(k2int);
      E[index] += Complex(w2/k,nu*w2);
      Complex Sij=Si[j];
      Complex& Tindex=T[index];
      Tindex.re += realproduct(Sij,wij);
      Forcing->Force(wij,Tindex.im,k);
    }
  };
  
  class FE {
    DNSBase *b;
    const vector& E;
    Real k0;
    
  public: 
    FE(DNSBase *b) : b(b), E(b->E), k0(b->k0) {}
    inline void operator()(const vector& wi, const vector& Si, unsigned I2,
                           unsigned j) {
      unsigned k2int=I2+j*j;
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
    inline void operator()(const vector& wi, const vector& Si, unsigned I2,
                           unsigned j) {
      unsigned k2int=I2+j*j;
      Real kint=sqrt(k2int);
      Real k=k0*kint;
      Complex wij=wi[j];
      Nu nu=b->nuk(k2int);
      Complex Sij=Si[j];
      double T;
      Forcing->Force(wij,T,k);
      Si[j]=Sij-nu*wij;
    }
  };
  
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);

  void ZeroT(const vector2& Src) {
    Set(T,Src[TRANSFER]);
    for(unsigned K=0; K < nshells; K++)
      T[K]=0.0;  
  }
  
  void ZeroE(const vector2& Src) {
    Set(E,Src[EK]);
    for(unsigned K=0; K < nshells; K++)
      E[K]=0.0;
  }
    
  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    ZeroT(Src);
    ZeroE(Src);
    Compute(FETL(this),Src,Y);
  }

  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum) {
      ZeroT(Src);
      Compute(FTL(this),Src,Y);
    }
  }

  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      ZeroT(Src);
      ZeroE(Src);
      Compute(FET(this),Src,Y);
    }
  }

  template<class T>
  void Compute(T fcn, const vector2& Src, const vector2& Y)
  {
    f0.Set(Src[OMEGA]);
    w.Set(Y[OMEGA]);

    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      vector Si=f0[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j)
        fcn(wi,Si,I2,j);
    }
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      ZeroT(Src);
      ZeroE(Src);
      Compute(FETL(this),Src,Y);
    } else
      Compute(FL(this),Src,Y);
  }

  Nu LinearCoeff(unsigned k) {
    unsigned i=k/my;
    unsigned j=k-my*i;
    return nuk(i*i+j*j);
  }

  // TODO: use a 1D lookup table on i^2+j^2.
  Real nuk(unsigned i2) {
    double k2=i2*k02;
    return nuL*pow(k2,pL)+nuH*pow(k2,pH);
  }

  virtual void ComputeInvariants(const array2<Complex>&, Real&, Real&, Real&);
  void Stochastic(const vector2& Y, double, double);

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

//***** initial conditions *****//

extern InitialConditionBase *InitialCondition;
class Zero : public InitialConditionBase {
public:
  const char *Name() {return "Zero";}
  void Set(Complex *w, unsigned n) {
    for(unsigned i=0; i < n; i++)
      w[i]=0.0;
  }
};

#endif
