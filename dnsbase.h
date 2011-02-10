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

static const double sqrt2=sqrt(2.0);


using namespace Array;
using namespace fftwpp;
using std::ostringstream;

extern unsigned spectrum;
extern unsigned casimir;

extern int pH;
extern int pL;

extern int circular;

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

  void CountAxes(array1<unsigned>::opt &C, unsigned I)  {
    C[(this->*Sindex)(I,0,I)] += 2;
  }
  void CountDiag(array1<unsigned>::opt &C, unsigned I)  {
    C[(this->*Sindex)(I,I,sqrt2*I)] += 2;
  }
  void CountMain(array1<unsigned>::opt &C, unsigned I, unsigned j) {
    C[(this->*Sindex)(I,j,sqrt(I*I+j*j))] += 4;
  }

  void LinearAxes(unsigned I,
		  Complex wa0,Complex wa1,Complex& fa0,Complex& fa1)  {
    Real nukk=nuk(I*I);
    fa0 -= nukk*wa0;
    fa1 -= nukk*wa1;
  }
  void LinearDiag(unsigned I,
		  Complex wd0,Complex wd1,Complex& fd0,Complex& fd1)  {
    Real nukk=nuk(2*(I*I));
    fd0 -= nukk*wd0;
    fd1 -= nukk*wd1;
  }
  void LinearMain(unsigned I,unsigned j,
		  Complex w0,Complex w1,Complex w2,Complex w3,
		  Complex& f0,Complex& f1,Complex& f2,Complex& f3) {
    Real nukk=nuk(I*I+j*j);
    f0 -= nukk*w0;
    f1 -= nukk*w1;
    f2 -= nukk*w2;
    f3 -= nukk*w3;
  }

  void NonLinearAxes(unsigned I,
		     Complex wa,Complex wb,
		     Complex& f0a,Complex& f0b,
		     Complex& f1a,Complex& f1b,
		     Complex& g0a,Complex& g0b,
		     Complex& g1a,Complex& g1b)  {
    Real k2inv=1.0/(k02*(I*I));

    // wx[I]
    Real ky=k0*I;
    Complex kyw=Complex(-ky*wa.im,ky*wa.re);
    f0a=0;
    f1a=kyw;
    g0a=k2inv*kyw;
    g1a=0;

    // w[xorigin+I][0]
    Real kx=k0*I;
    Complex kxw=Complex(-kx*wb.im,kx*wb.re);
    f0b=kxw;
    f1b=0;
    g0b=0;
    g1b=-k2inv*kxw;
  }

  void NonLinearDiag(unsigned I,
		     Complex wa,Complex wb,
		     Complex& f0a,Complex& f0b,
		     Complex& f1a,Complex& f1b,
		     Complex& g0a,Complex& g0b,
		     Complex& g1a,Complex& g1b)  {
    Real k2inv=1.0/(k02*(2*I*I));

    // wi[I]
    Real kx=k0*I;
    Real ky=k0*I;
    Complex kxw=Complex(-kx*wa.im,kx*wa.re);
    Complex kyw=Complex(-ky*wa.im,ky*wa.re);
    f0a=kxw;
    f1a=kyw;
    g0a=k2inv*kyw;
    g1a=-k2inv*kxw;

    // wim[I]
    kx=-kx;
    kxw=Complex(-kx*wb.im,kx*wb.re);
    kyw=Complex(-ky*wb.im,ky*wb.re);
    f0b=kxw;
    f1b=kyw;
    g0b=k2inv*kyw;
    g1b=-k2inv*kxw;
  }

  void NonLinearMain(unsigned I,unsigned j,
		     Complex wa,Complex wb,Complex wc,Complex wd,
		     Complex& f0a,Complex& f0b,Complex& f0c,Complex& f0d,
		     Complex& f1a,Complex& f1b,Complex& f1c,Complex& f1d,
		     Complex& g0a,Complex& g0b,Complex& g0c,Complex& g0d,
		     Complex& g1a,Complex& g1b,Complex& g1c,Complex& g1d) {
    Real k2inv=1.0/(k02*(I*I+j*j));

    // a wi[j]
    Real kx=k0*I;
    Real ky=k0*j;
    Complex kxw=Complex(-kx*wa.im,kx*wa.re);
    Complex kyw=Complex(-ky*wa.im,ky*wa.re);
    f0a=kxw;
    f1a=kyw;
    g0a=k2inv*kyw;
    g1a=-k2inv*kxw;

    // b wim[j]
    kx=-k0*I;
    kxw=Complex(-kx*wb.im,kx*wb.re);
    kyw=Complex(-ky*wb.im,ky*wb.re);
    f0b=kxw;
    f1b=kyw;
    g0b=k2inv*kyw;
    g1b=-k2inv*kxw;

    // c w[xorigin-j][I]
    kx=-k0*j;
    ky=k0*I;
    kxw=Complex(-kx*wc.im,kx*wc.re);
    kyw=Complex(-ky*wc.im,ky*wc.re);
    f0c=kxw;;
    f1c=kyw;
    g0c=k2inv*kyw;
    g1c=-k2inv*kxw;

    // d w[xorigin+j][I]
    kx=-kx;
    kxw=Complex(-kx*wd.im,kx*wd.re);
    kyw=Complex(-ky*wd.im,ky*wd.re);
    f0d=kxw;
    f1d=kyw;
    g0d=k2inv*kyw;
    g1d=-k2inv*kxw;
  }


  void SpectrumAxes(unsigned I,vector& S,Complex wa0,Complex wa1)  {
    unsigned I2=I*I;
    Real Wall=abs2(wa0)+abs2(wa1);
    S[(this->*Sindex)(I,0,I)] += Complex(Wall/(k0*I),nuk(I2)*Wall);
    //S[kval[I][0]] += Complex(Wall/(k0*I),nuk(I2)*Wall);
  }
  void SpectrumDiag(unsigned I, vector& S,Complex wd0,Complex wd1)  {
    unsigned I2=I*I;
    Real Wall=abs2(wd0)+abs2(wd1);
    Real k=sqrt2*I;
    S[(this->*Sindex)(I,I,k)] += Complex(Wall/(k0*k),nuk(I2+I2)*Wall);
    //S[kval[I][I]] += Complex(Wall/(k0*k),nuk(2*I2)*Wall);
  }    
  void SpectrumMain(unsigned I,unsigned j, vector& S,
		    Complex w0,Complex w1,Complex w2,Complex w3) {
    Real Wall=abs2(w1)+abs2(w2)+abs2(w2)+abs2(w3);
    unsigned k2=I*I+j*j;
    Real k=sqrt(k2);
    S[(this->*Sindex)(I,j,k)] += Complex(Wall/(k0*k),nuk(k2)*Wall);
    //S[kval[I][j]] += Complex(Wall/(k0*k),nuk(k2)*Wall);
  }

  void TransferAxes(unsigned I, vector& T,
		    Complex wa0,Complex wa1,Complex fa0,Complex fa1)  {
    //T[kval[I][0]] += realproduct(wa0,fa0) + realproduct(wa1,fa1);
    T[(this->*Sindex)(I,0,I)] += realproduct(wa0,fa0) + realproduct(wa1,fa1);
  }
  void TransferDiag(unsigned I,vector& T,
		    Complex wd0,Complex wd1,Complex fd0,Complex fd1)  {
    //T[kval[I][I]] += realproduct(wd0,fd0) + realproduct(wd1,fd1);
    T[(this->*Sindex)(I,0,I)] += realproduct(wd0,fd0) + realproduct(wd1,fd1);
  }
  void TransferMain(unsigned I,unsigned j, vector& T,
		    Complex w0,Complex w1,Complex w2,Complex w3,
		    Complex f0,Complex f1,Complex f2,Complex f3) {
    T[(this->*Sindex)(I,j,sqrt(I*I+j*j))] 
      //T[kval[I][j]]
      += realproduct(w0,f0) + realproduct(w1,f1) 
      +  realproduct(w2,f2) + realproduct(w3,f3);
  }

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
  void setcountRAW();
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

//***** loop class *****//
typedef void (DNSBase::*Ca)(array1<unsigned>::opt&,unsigned);
typedef void (DNSBase::*Cm)(array1<unsigned>::opt&,unsigned,unsigned);

typedef void (DNSBase::*Sa)(unsigned,vector&,Complex,Complex);
typedef void (DNSBase::*Sm)(unsigned,unsigned,
			    vector&,Complex,Complex,Complex,Complex);

typedef void (DNSBase::*Ta)(unsigned,vector&,Complex,Complex,Complex,Complex);
typedef void (DNSBase::*Tm)(unsigned,unsigned,vector&,
			    Complex,Complex,Complex,Complex,
			    Complex,Complex,Complex,Complex);

typedef void (DNSBase::*La)(unsigned,Complex,Complex,Complex&,Complex&);
typedef void (DNSBase::*Lm)(unsigned,unsigned,
			    Complex,Complex,Complex,Complex,
			    Complex&,Complex&,Complex&,Complex&);

typedef void (DNSBase::*Na)(unsigned,Complex,Complex,
			    Complex&,Complex&,
			    Complex&,Complex&,
			    Complex&,Complex&,
			    Complex&,Complex&);
typedef void (DNSBase::*Nm)(unsigned,unsigned,
			    Complex,Complex,Complex,Complex,
			    Complex&,Complex&,Complex&,Complex&,
			    Complex&,Complex&,Complex&,Complex&,
			    Complex&,Complex&,Complex&,Complex&,
			    Complex&,Complex&,Complex&,Complex&);

class Hloop{
 private:
  unsigned N, m, xorigin;
  DNSBase *parent;
  unsigned radix, Invisible, m2;
  bool subgrid;
 public:
  Hloop() {radix=1; subgrid=false; Invisible=0;}
  Hloop(DNSBase *parent0) {
    setparams(parent0);
    radix=1;
    subgrid=false;
    Invisible=0;
  }
  Hloop(unsigned N0) {
    setparams(N);
    radix=1;
    subgrid=false;
    Invisible=0;
  }
  ~Hloop() {}

  void setparams(unsigned N0) {
    N=N0;
    m=(N+1)/2;
    xorigin=m-1;
    m2=(m-1)*(m-1);
  }
  void setparams(DNSBase *parent0) {
    parent=parent0;
    N=parent->getNx();
    m=parent->getmx();
    xorigin=parent->getxorigin();
    m2=(m-1)*(m-1);
  }
  void setradix(unsigned radix0) {radix=radix0;}
  void makesubgrid() {subgrid=true;}
  void setInvisible(unsigned inv) {Invisible=inv;}

  bool doaxes(unsigned I2) {return circular ? I2 <= m2: true;}
  bool dodiag(unsigned I2) {return circular ? 2*I2 <= m2: true;}
  unsigned jstop(unsigned I) {
    return  circular ? min((int)I,(int) ceil(sqrt(m2-I*I))) : I;
  }

  virtual bool isvisible(unsigned i, unsigned j) {return true;}
  
  void Rloop(DynVector<unsigned> &R2) {
    DynVector<unsigned> temp;
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      if(doaxes(I2)) temp.Push(I2);
      if(dodiag(I2)) temp.Push(2*I2);
      unsigned stop =jstop(I);
      for(unsigned j=1; j < stop; ++j) {
	temp.Push(I2+j*j);
      }
    }
    temp.sort();
    R2.Push(temp[0]);
    unsigned last=0;
    for(unsigned i=1; i < temp.Size(); ++i) {
      if(temp[i] != R2[last])
	R2[++last]=temp[i];
    }
  }

  void Cloop(array1<unsigned>::opt &C, Ca afp, Ca dfp, Cm mfp) {
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      if(doaxes(I2)) (parent->*afp)(C,I);
      if(dodiag(I2)) (parent->*dfp)(C,I);
      unsigned stop =jstop(I);
      for(unsigned j=1; j < stop; ++j) {
	(parent->*mfp)(C,I,j);
      }
    }
  }

  // loop over the Hermitian-symmetric array2 w, calculate
  // something, put it into the appropriate index of S, a spectrum-like array
  void Sloop(vector& S, const array2<Complex> w,Sa  afp,Sa  dfp, Sm mfp) {
    vector wx=w[xorigin];
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      vector wxi=w[xorigin+I];
      unsigned stop =jstop(I);
      for(unsigned j=1; j < stop; ++j)
	(parent->*mfp)(I,j,S,wi[j],wim[j],w[xorigin-j][I],w[xorigin+j][I]);
      if(doaxes(I2)) (parent->*afp)(I,S,wx[I],wxi[0]);
      if(dodiag(I2)) (parent->*dfp)(I,S,wi[I],wim[I]);
    }
  }
  
  void Tloop(vector& T, const array2<Complex> w, const array2<Complex> f,
	     Ta afp, Ta dfp, Tm mfp) {
    vector wx=w[xorigin];
    vector fx=f[xorigin];
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      vector fi=f[i];
      vector fim=f[im];
      unsigned stop =jstop(I);
      for(unsigned j=1; j < stop; ++j) {
	(parent->*mfp)(I,j,T,
		       wi[j],wim[j],w[xorigin-j][I],w[xorigin+j][I],
		       fi[j],fim[j],f[xorigin-j][I],f[xorigin+j][I]);
      }
      if(doaxes(I2))
	(parent->*afp)(I,T,wx[I],w[xorigin+I][0],fx[I],f[xorigin+I][0]);
      if(dodiag(I2))
	(parent->*dfp)(I,T,wi[I],wim[I],fi[I],fim[I]);
    }
  }

  void Nloop(const array2<Complex> w, 
	     array2<Complex> f0, array2<Complex> f1, 
	     array2<Complex> g0, array2<Complex> g1,
	     Na afp, Na dfp, Nm mfp) {
    vector wx=w[xorigin];
    vector f0x=f0[xorigin];
    vector f1x=f1[xorigin];
    vector g0x=g0[xorigin];
    vector g1x=g1[xorigin];
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      vector f0i=f0[i];
      vector f0im=f0[im];
      vector f1i=f1[i];
      vector f1im=f1[im];
      vector g0i=g0[i];
      vector g0im=g0[im];
      vector g1i=g1[i];
      vector g1im=g1[im];
      unsigned stop =jstop(I);
      for(unsigned j=1; j < stop; ++j) {
	(parent->*mfp)(I,j,
		       wi[j],wim[j],w[xorigin-j][I],w[xorigin+j][I],
		       f0i[j],f0im[j],f0[xorigin-j][I],f0[xorigin+j][I],
		       f1i[j],f1im[j],f1[xorigin-j][I],f1[xorigin+j][I],
		       g0i[j],g0im[j],g0[xorigin-j][I],g0[xorigin+j][I],
		       g1i[j],g1im[j],g1[xorigin-j][I],g1[xorigin+j][I]);
      }
      if(doaxes(I2))
	(parent->*afp)(I,wx[I],w[xorigin+I][0],
		       f0x[I],f0[xorigin+I][0],
		       f1x[I],f1[xorigin+I][0],
		       g0x[I],g0[xorigin+I][0],
		       g1x[I],g1[xorigin+I][0]);
      if(dodiag(I2))
	(parent->*dfp)(I,wi[I],wim[I],
		       f0i[I],f0im[I],
		       f1i[I],f1im[I],
		       g0i[I],g0im[I],
		       g1i[I],g1im[I]);
    }
  }
  
  void Lloop(const array2<Complex> w, array2<Complex> f,
	     La afp, La dfp, Lm mfp) {
    vector wx=w[xorigin];
    vector fx=f[xorigin];
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      vector fi=f[i];
      vector fim=f[im];
      unsigned stop =jstop(I);
      for(unsigned j=1; j < stop; ++j) {
	(parent->*mfp)(I,j,
		       wi[j],wim[j],w[xorigin-j][I],w[xorigin+j][I],
		       fi[j],fim[j],f[xorigin-j][I],f[xorigin+j][I]);
      }
      if(doaxes(I2))
	(parent->*afp)(I,wx[I],w[xorigin+I][0],fx[I],f[xorigin+I][0]);
      if(dodiag(I2))
	(parent->*dfp)(I,wi[I],wim[I],fi[I],fim[I]);
    }
  }

  
};


#endif
