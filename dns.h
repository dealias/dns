#ifndef __dns_h__
#define __dns_h__ 1

#include "Array.h"
#include "kernel.h"
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

//class DNS : public ProblemBase {
class DNSBase {
 protected:
  // Vocabulary:
  int xpad, ypad; // these are always one.
  Real nuH, nuL;
  int pH;
  int pL;
  unsigned Nx;
  unsigned Ny;
  unsigned spectrum;
 
  // derived variables:
  enum Field {OMEGA,TRANSFER,EK};
  unsigned mx, my; // size of data arrays
  unsigned origin; // linear index of Fourier origin.
  unsigned xorigin; // horizontal index of Fourier origin.

  Real k0; // grid spacing factor
  Real k02; // k0^2
  array2<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field;
  vector T;

  int tcount;
  array1<unsigned>::opt count;

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
  
 public:
  DNSBase() {}
  DNSBase(unsigned Nx, unsigned my, Real k0): Nx(Nx), my(my), k0(k0) {
    block=ComplexAlign(3*Nx*my);
    mx=(Nx+1)/2;
    xorigin=mx-1;
    origin=xorigin*my;
    w.Dimension(Nx,my);
    f0.Dimension(Nx,my);
    f1.Dimension(Nx,my,block);
    g0.Dimension(Nx,my,block+Nx*my);
    g1.Dimension(Nx,my,block+2*Nx*my);
    F[0]=f0;
    F[1]=f1;
    G[0]=g0;
    G[1]=g1;
    Convolution=new ImplicitHConvolution2(mx,my,2);
  }
  ~DNSBase() {}

  unsigned getNx() {return Nx;}
  unsigned getmx() {return mx;}
  unsigned getmy() {return my;}
  Real getk0() {return k0;}
  Real getk02() {return k02;}
  unsigned getxorigin() {return xorigin;}

  void InitialConditions();
  void Initialize();
  //  virtual void Output(int it)=0;
  void FinalOutput();
  void OutFrame(int it);

  virtual void Spectrum(vector& S, const vector& y);
  void Transfer(const vector2& Src, const vector2& Y);
  void NonLinearSource(const vector& Src, const vector& Y, double t);
  void LinearSource(const vector& Src, const vector& Y, double t);

  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src[OMEGA],Y[OMEGA],t);
    if(spectrum) Transfer(Src,Y);
    LinearSource(Src[OMEGA],Y[OMEGA],t);
  }

  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum) Spectrum(Src[EK],Y[OMEGA]);
    HermitianSymmetrizeX(mx,my,xorigin,Src[OMEGA]);
  }

  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src[OMEGA],Y[OMEGA],t);
    if(spectrum) Transfer(Src,Y);
    NonConservativeSource(Src,Y,t);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    ConservativeSource(Src,Y,t);
    NonConservativeSource(Src,Y,t);
  }
  Nu LinearCoeff(unsigned k) {
    unsigned i=k/my;
    unsigned j=k-my*i;
    return nuk(k02*(i*i+j*j));
  }

  // TODO: use a 1D lookup table on i^2+j^2.
  Real nuk(Real k2) {
    return nuL*pow(k2,pL)+nuH*pow(k2,pH);
  }

  void ComputeInvariants(array2<Complex> &w,Real& E, Real& Z, Real& P);
  void Stochastic(const vector2& Y, double, double);

  Real Spectrum(unsigned int i) {
    double c=count[i];
    return c > 0 ? T[i].re*twopi/c : 0.0;
  }

  Real Dissipation(unsigned int i) {
    return T[i].im;
  }

  Real Pi(unsigned int i) {
    return T[i].re;
  }
  Real Eta(unsigned int i) {
    return T[i].im;
  }
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

//***** Source routines *****//

void DNSBase::LinearSource(const vector& wSrc, const vector& w0, double)
{
  w.Set(w0);
  f0.Set(wSrc);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector f0i=f0[i];
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j)
      f0i[j] -= nuk(k02*(I2+j*j))*wi[j];
  }
}

void DNSBase::NonLinearSource(const vector& wSrc, const vector& wY, double)
{
  w.Set(wY);
  f0.Set(wSrc);

  f0(origin)=0.0;
  f1(origin)=0.0;
  g0(origin)=0.0;
  g1(origin)=0.0;

  for(unsigned i=0; i < Nx; ++i) {
    Real kx=k0*((int) i-(int) xorigin);
    Real kx2=kx*kx;
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    vector g0i=g0[i];
    vector g1i=g1[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      Complex wij=wi[j];
      Complex kxw=Complex(-kx*wij.im,kx*wij.re);
      Complex kyw=Complex(-ky*wij.im,ky*wij.re);
      f0i[j]=kxw;
      f1i[j]=kyw;
      Real k2inv=1.0/(kx2+ky*ky);
      g0i[j]=k2inv*kyw;
      g1i[j]=-k2inv*kxw;
    }
  }
  F[0]=f0;
  //  Convolution->convolve(F,G);
  f0(origin)=0.0;

#if 0
  Real sum=0.0;
  for(unsigned i=0; i < Nx; ++i) {
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Complex wij=wi[j];
      sum += (f0[i][j]*conj(wij)).re;
    }
  }
  cout << sum << endl;
#endif

}

void DNSBase::Transfer(const vector2& Src, const vector2& Y)
{
  Set(T,Src[TRANSFER]);

  for(unsigned K=0; K < nshells; K++)
    T[K]=0.0;
  f0.Set(Src[OMEGA]);

  w.Set(Y[OMEGA]);
  Var factor=sqrt(2.0*dt)*crand_gauss();

  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    vector Si=f0[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real k=k0*sqrt(I2+j*j);
      T[(unsigned)(k-0.5)].re += realproduct(Si[j],wi[j]);
    }
  }

  Forcing->Force(f0,T);
}

void DNSBase::Spectrum(vector& S, const vector& y)
{
  w.Set(y);
  for(unsigned K=0; K < nshells; K++)
    S[K]=0.0;

  // Compute instantaneous angular sum over each circular shell.
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real k2=k02*(I2+j*j);
      Real k=sqrt(k2);
      Real w2=abs2(wi[j]);
      S[(unsigned)(k-0.5)] += Complex(w2/k,nuk(k2)*w2);
    }
  }
}

void DNSBase::Stochastic(const vector2&Y, double, double dt)
{
  w.Set(Y[OMEGA]);
  Forcing->Force(w,sqrt(2.0*dt)*crand_gauss());
}

//***** DNSBase Output routines *****//

void DNSBase::Initialize()
{
  fevt << "# t\tE\tZ\tP" << endl;
  if(spectrum) {
    for(unsigned i=0; i < nshells; i++)
      count[i]=0;
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
        count[(unsigned)(sqrt(I2+j*j)-0.5)]++;
      }
    }
  }
}

void DNSBase::OutFrame(int)
{
  //  w.Set(Y[OMEGA]); // FIXME
  unsigned int Nx0=Nx+xpad;
  unsigned int Ny0=Ny+ypad;
  unsigned int offset=Nx0/2-mx+1;
  for(unsigned int i=0; i < Nx; ++i) {
    unsigned int I=i+offset;
    for(unsigned int j=0; j < my; j++)
      buffer(I,j)=w(i,j);
  }

  Padded->pad(buffer);
  Padded->backwards(buffer,true);

  fw << 1 << Ny0 << Nx0;
  for(int j=Ny0-1; j >= 0; j--) {
    for(unsigned i=0; i < Nx0; i++) {
      fw << (float) wr(i,j);
    }
  }
  fw.flush();
}

void DNSBase::ComputeInvariants(array2<Complex> &w, Real& E, Real& Z, Real& P)
{
  E=Z=P=0.0;
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real w2=abs2(wi[j]);
      Z += w2;
      Real k2=k02*(I2+j*j);
      E += w2/k2;
      P += k2*w2;
    }
  }
}

void DNSBase::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(w,E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

#endif
