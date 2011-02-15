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
  static Real nuH, nuL;
  static const int xpad,ypad;
  
  enum Field {OMEGA,TRANSFER,EK};

  // derived variables:
  unsigned mx, my; // size of data arrays
  unsigned origin; // linear index of Fourier origin.
  unsigned xorigin; // horizontal index of Fourier origin.

  static Real k0; // grid spacing factor
  static Real k02; // k0^2
  array2<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field;

  int tcount;
public:  
  static Real etanorm;

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
  
  ImplicitHTConvolution2 *TConvolution;
  oxstream ftransferN;
  Array2<Complex> f,g,h;
  vector Tn;
  array1<unsigned>::opt count;
  static vector E; // Spectrum
  static vector T; // Transfer

public:
  DNSBase() {}
  virtual ~DNSBase() {}

  unsigned getNx() {return Nx;}
  unsigned getmx() {return mx;}
  unsigned getmy() {return my;}
  Real getk0() {return k0;}
  Real getk02() {return k02;}
  unsigned getxorigin() {return xorigin;}
  Real getetanorm() {return etanorm;}

  void InitialConditions();
  void Initialize();
  virtual void setcount();
  //  virtual void Output(int it)=0;
  void FinalOutput();
  void OutFrame(int it);

  typedef void (*SourceFcn)(const vector& wi, const vector& Si, unsigned I2,
                            unsigned j);
  
  static void FETL(const vector& wi, const vector& Si, unsigned I2, unsigned j);
  static void FTL(const vector& wi, const vector& Si, unsigned I2, unsigned j);
  static void FET(const vector& wi, const vector& Si, unsigned I2, unsigned j);
  static void FE(const vector& wi, const vector& Si, unsigned I2, unsigned j);
  static void FL(const vector& wi, const vector& Si, unsigned I2, unsigned j);

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
    Compute(FETL,Src,Y);
  }

  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum) {
      ZeroT(Src);
      Compute(FTL,Src,Y);
    }
  }

  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      ZeroT(Src);
      ZeroE(Src);
      Compute(FET,Src,Y);
    }
  }

  void Compute(SourceFcn fcn, const vector2& Src, const vector2& Y);
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      ZeroT(Src);
      ZeroE(Src);
      Compute(FETL,Src,Y);
    } else
      Compute(FL,Src,Y);
  }

  Nu LinearCoeff(unsigned k) {
    unsigned i=k/my;
    unsigned j=k-my*i;
    return nuk(i*i+j*j);
  }

  // TODO: use a 1D lookup table on i^2+j^2.
  static Real nuk(unsigned i2) {
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
