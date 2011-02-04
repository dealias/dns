#ifndef __dnsbase_h__
#define __dnsbase_h__ 1

#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "ArrayL.h"
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
extern unsigned casimir;

extern int pH;
extern int pL;

extern int circular;

class Hloop{
 private:
  unsigned Nx, my;
  unsigned mx, xorigin;
  Real k0;
  //unsigned (*Sindex)(unsigned, unsigned, Real);
  //void (*fptr)(Complex&,Complex&);
 public:
  Hloop() {};
  Hloop(unsigned Nx0, unsigned my0, Real k00) {
    Nx=Nx0;
    my=my0;
    k0=k00;
    mx=(Nx+1)/2;
    xorigin=mx-1;
  };
  ~Hloop() {};

  // loop over the Hermitian-symmetric array2 w, calculate
  // something, put it into the appropriate index of S.
  void Sloop(vector& S, const array2<Complex> w,
	     void (*fptr)(vector&,Complex)) {
    for(unsigned i=0; i < Nx; ++i) {
      vector wi=w[i];
      for(unsigned j=i > xorigin ? 1 : 0; j < my; ++j) {
	// TODO: add visible/invisible modes
	fptr(S,wi[j]);
      }
    }
  }

};

class DNSBase {
protected:
  // Vocabulary:
  unsigned Nx;
  unsigned Ny;
  Real nuH, nuL;
  static const int xpad,ypad;
  
  enum Field {OMEGA,TRANSFER,TRANSFERN,EK};
  enum SPEC {NOSPECTRUM, BINNED, INTERPOLATED, RAW}; 

  // derived variables:
  unsigned mx, my; // size of data arrays
  unsigned origin; // linear index of Fourier origin.
  unsigned xorigin; // horizontal index of Fourier origin.

  Real k0; // grid spacing factor
  Real k02; // k0^2
  array2<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field;

  int tcount;
  Real etanorm;

  unsigned (DNSBase::*Sindex)(unsigned, unsigned, Real); 
  unsigned SkBIN(unsigned I, unsigned j, Real k) {return (unsigned)(k-0.5);}
  unsigned SkRAW(unsigned I, unsigned j, Real k) {return  kval[I][j];}
  
  unsigned nmode;
  unsigned nshells;  // Number of spectral shells
  array1<unsigned> R2; //radii achieved for discrete spectrum
  array2L<unsigned> kval; // for discrete spectrum
  array2L<unsigned> area, areadown; // for interpolated spectrum


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
  
  void CasimirTransfer(const vector2& Src, const vector2& Y);
  ImplicitHTConvolution2 *TConvolution;
  oxstream ftransferN;
  Array2<Complex> f,g,h;
  vector Tn;

  void setSindex() {
    switch(spectrum) {
    case NOSPECTRUM:
      Sindex=NULL;
      break;
    case BINNED:
      Sindex=&DNSBase::SkBIN;
      break;
    case INTERPOLATED:
      msg(ERROR,"Interpolated spectrum not done yet.");
      break;
    case RAW:
      Sindex=&DNSBase::SkRAW;
      break;
    default:
      msg(ERROR,"Invalid choice of spectrum.");
    }
  }
  
  void check_rvn(DynVector<unsigned> & R2, const unsigned r2, 
		 const unsigned first)  {
    bool found=false;
    
    unsigned last=R2.Size();
    for(unsigned j=first; j < last; ++j) {
      if(r2 == R2[j]) {
	found=true;
	break;
      }
    }
    if(!found)
      R2.Push(r2);
  }
  

public:
  DNSBase() {setSindex();}
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
  
    setSindex();
  }
  virtual ~DNSBase() {}

  unsigned getNx() {return Nx;}
  unsigned getmx() {return mx;}
  unsigned getmy() {return my;}
  Real getk0() {return k0;}
  Real getk02() {return k02;}
  unsigned getxorigin() {return xorigin;}
  Real getetanorm() {return etanorm;}
  unsigned getkval(const unsigned i, const unsigned j) {return kval[i][j];}
  unsigned getR2(const unsigned i) {return R2[i];}

  void InitialConditions();
  void Initialize();
  virtual void setcount();
  void setcountBINNED();
  void setcountRAW(unsigned lambda2=1);
  //  virtual void Output(int it)=0;
  void FinalOutput();
  void OutFrame(int it);

  virtual void Spectrum(vector& S, const vector& y);
  void Transfer(const vector2& Src, const vector2& Y);
  void NonLinearSource(const vector& Src, const vector& Y, double t);
  void LinearSource(const vector& Src, const vector& Y, double t);

  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src[OMEGA],Y[OMEGA],t);
    if(spectrum != NOSPECTRUM) Transfer(Src,Y);
    LinearSource(Src[OMEGA],Y[OMEGA],t);
  }

  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum != NOSPECTRUM) Spectrum(Src[EK],Y[OMEGA]);
//    HermitianSymmetrizeX(mx,my,xorigin,Src[OMEGA]);
  }

  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src[OMEGA],Y[OMEGA],t);
    if(spectrum != NOSPECTRUM) Transfer(Src,Y);
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

  // TODO: consider using a lookup table on i2.
  Real nuk(unsigned i2) {
    double k2=i2*k02;
    return nuL*pow(k2,pL)+nuH*pow(k2,pH);
  }

  virtual void ComputeInvariants(const array2<Complex>&, Real&, Real&, Real&);
  void Stochastic(const vector2& Y, double, double);

  array1<unsigned>::opt count;
  vector T; // Transfer
  virtual Real getSpectrum(unsigned i) {
    double c=count[i];
    return c > 0 ? T[i].re*twopi/c : 0.0;
  }
  Real Dissipation(unsigned i) {return T[i].im;}
  Real Pi(unsigned i) {return T[i].re;}
  Real Eta(unsigned i) {return T[i].im;}
  Real kb(unsigned i) {
    if(spectrum == RAW)
      return i == 0 ? 0.5*k0 : k0*sqrt((Real) R2[i-1]);
    return k0*(i+0.5);
  }
  Real kc(unsigned i) {
    if(spectrum == RAW) 
      return k0*sqrt((Real) R2[i]);
    return k0*(i+1);
  }

  void findrads(DynVector<unsigned> &R2, array1<unsigned> nr, 
		unsigned m=0, Real lambda=0, unsigned Invisible=0)
  {
    for(unsigned i=1; i < my; ++i) {
      //unsigned start=nr[(unsigned) floor(sqrt((i-1)/2))];
      double nrstart=floor(sqrt((i-1)/2));
      unsigned start=nr[nrstart > 1 ? (unsigned) nrstart -1 : 0];
      start=0; // FIXME: restore and optimize.
      for(unsigned x=i-1; x <= i; ++x) {
	unsigned x2=x*x;
	unsigned ystopnow= i < x ? i : x;
	for(unsigned y= x == 0 ? 1 : 0; y <= ystopnow; ++y) {
	  if(isvisible(x,y,m,lambda,Invisible)) {
	    check_rvn(R2,x2+y*y,start);
	  }
	}
      }
      nr[i]=R2.Size();
    }
  }

  void killmodes(array2<Complex> &A) {
    if(circular) {
      unsigned m2=my*my;
      for(unsigned i=0; i < Nx; ++i) {
	vector Ai=A[i];
	unsigned I= i > xorigin ? i-xorigin : xorigin-i; 
	unsigned start=(unsigned) ceil(sqrt(m2-I*I));
	for(unsigned j=start; j < my; ++j) {
	  Ai[j]=0.0;
	}
      }
    }
  }

  virtual bool isvisible(unsigned I, unsigned j, 
			 unsigned m=0, Real lambda=0, unsigned Invsible=0) {
    if(circular) 
      return I*I + j*j <= my*my;
    return true;
  }

  virtual unsigned diagstart() {return 1;}
  virtual unsigned diagstop() {
    if(circular) 
      return (unsigned) ceil((my-1)/sqrt(2.0));
    return mx;
  }
  virtual unsigned mainjstart(unsigned I) {return 1;}
  virtual unsigned mainjstop(unsigned I) {
    if(circular) 
      return min(I,(unsigned) ceil(sqrt((my-1)*(my-1)-I*I)));
    return I;
  }
  virtual unsigned xoriginstart() {return 1;}
  virtual unsigned xoriginstop() {return my;}
  virtual unsigned bottomstart() {return 1;}
  virtual unsigned bottomstop() {return mx;}
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
