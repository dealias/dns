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

  void SpectrumAxes(unsigned I,vector& S,Complex wa0,Complex wa1)  {
    unsigned I2=I*I;
    Real Wall=abs2(wa0)+abs2(wa1);
    S[(this->*Sindex)(I,0,I)] += Complex(Wall/(k0*I),nuk(I2)*Wall);
  }
  void SpectrumDiag(unsigned I, vector& S,Complex wd0,Complex wd1)  {
    unsigned I2=I*I;
    Real Wall=abs2(wd0)+abs2(wd1);
    Real k=sqrt2*I;
    S[(this->*Sindex)(I,I,k)] += Complex(Wall/(k0*k),nuk(I2+I2)*Wall);
  }    
  void SpectrumMain(unsigned I,unsigned j, vector& S,
		    Complex w0,Complex w1,Complex w2,Complex w3) {
    Real Wall=abs2(w1)+abs2(w2)+abs2(w2)+abs2(w3);
    unsigned k2=I*I+j*j;
    Real k=sqrt(k2);
    S[(this->*Sindex)(I,j,k)] += Complex(Wall/(k0*k),nuk(k2)*Wall);
  }

  void TransferAxes(unsigned I, vector& T,
		    Complex wa0,Complex wa1,Complex fa0,Complex fa1)  {
    T[(this->*Sindex)(I,0,I)] += realproduct(wa0,fa0) + realproduct(wa1,fa1);
  }
  void TransferDiag(unsigned I,vector& T,
		    Complex wd0,Complex wd1,Complex fd0,Complex fd1)  {
    T[(this->*Sindex)(I,0,I)] += realproduct(wd0,fd0) + realproduct(wd1,fd1);
  }
  void TransferMain(unsigned I,unsigned j, vector& T,
		    Complex w0,Complex w1,Complex w2,Complex w3,
		    Complex f0,Complex f1,Complex f2,Complex f3) {
    T[(this->*Sindex)(I,j,sqrt(I*I+j*j))] 
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
    if(spectrum != NOSPECTRUM) 
      Transfer(Src,Y); // FIXME: called more oft than it should be with C_RKs
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
    return nuk(i*i+j*j);
  }

  // TODO: consider using a lookup table on i2.
  Real nuk(unsigned i2) {
    double k2=i2*k02;
    return nuL*pow(k2,pL)+nuH*pow(k2,pH);
  }

  virtual void ComputeInvariants(const array2<Complex>&, Real&, Real&, Real&);
  void AddInvariants(unsigned k2int, Real w2, Real& E, Real& Z, Real& P) {
    Z += w2;
    Real k2=k2int*k02;
    E += w2/k2;
    P += k2*w2;
  }
  
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

// typedefs for pointers to functions in DNSBase
typedef void (DNSBase::*Ca)(array1<unsigned>::opt&,unsigned);
typedef void (DNSBase::*Cm)(array1<unsigned>::opt&,unsigned,unsigned);

typedef void (DNSBase::*Invarfp)(unsigned,Real,Real&,Real&,Real&);

typedef void (DNSBase::*Sa)(unsigned,vector&,Complex,Complex);
typedef void (DNSBase::*Sm)(unsigned,unsigned,
			    vector&,Complex,Complex,Complex,Complex);

typedef void (DNSBase::*Ta)(unsigned,vector&,Complex,Complex,Complex,Complex);
typedef void (DNSBase::*Tm)(unsigned,unsigned,vector&,
			    Complex,Complex,Complex,Complex,
			    Complex,Complex,Complex,Complex);
// TODO: go back to varargs again?

class Hloop{
 private:
  unsigned N, m, xorigin;
  DNSBase *parent;
  unsigned radix, Invisible, m2;
  bool subgrid;
  Real innerradius, innerradius2;
  unsigned Astart, Astop, Dstart, Dstop; 
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
    innerradius=0;
    subgrid=false;
    bounds();
  }
  void setparams(DNSBase *parent0) {
    parent=parent0;
    N=parent->getNx();
    m=parent->getmx();
    xorigin=parent->getxorigin();
    m2=(m-1)*(m-1);
    innerradius=0;
    subgrid=false;
    bounds();
  }
  void setinnerradius(unsigned nold) {
    Real a=max((nold-1.0),0.0);
    innerradius=sqrt(a*a/radix);
    bounds();
  }
  void setradix2() {radix=2;}
  void setradix4() {radix=4;}
  void makesubgrid() {subgrid=true;}
  void setInvisible(unsigned inv) {Invisible=inv;}

  void bounds() {
    innerradius2=innerradius*innerradius;
    Astop=m;
    Dstart=1;

    if(!subgrid) { // standard DNS
      Astart=1;
      if(circular) 
	Dstop=(unsigned) ceil(((Real) m)/sqrt2);
      else 
	Dstop=m;
      
    } else { // this is MDNS
      if(circular) { 
	if(Invisible == 0) // we are the first grid
	  Astart=1;
	else 
	  Astart=0; // FIXME: something to do with innerradius
      } else {
	if(radix == 2) {
	  Astart=Invisible/2;
	} else { // radix == 1 or radix == 4
	  Astart=Invisible;
	}
      }
    }
  }
 unsigned jstart(unsigned I) {
    if(subgrid) {
      if(circular) {
	msg(ERROR,"gotta enable circular in Hloop still");
	return 1+ (unsigned) floor(sqrt(innerradius2 - I*I));
      }
      if(radix == 2)
	return Invisible-I;
      if(radix == 4)
	return I >= Invisible ? 1 : Invisible;
    }
    return 1; // not a subgrid
  }
  unsigned jstop(unsigned I) {
    return  circular ? min((int)I,(int) ceil(sqrt(m2-I*I))) : I;
  }
  unsigned jkillstart(unsigned I) {
    if(!circular) return m; // no kill I
    if(I < Dstop) return m;
    return (int) ceil(sqrt(m2-I*I));
  }

  // FIXME: not sure if I'm going to use this or not.
  virtual bool isvisible(unsigned i, unsigned j) {return true;}
  
  void Rloop(DynVector<unsigned> &R2) {
    DynVector<unsigned> temp;
    for(unsigned I=Astart; I < Astop; ++I)
      temp.Push(I*I);
    for(unsigned I=Dstart; I < Dstop; ++I)
      temp.Push(2*I*I);
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      unsigned start =jstart(I);
      unsigned stop =jstop(I);
      for(unsigned j=start; j < stop; ++j) {
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
    //cout << "R2 " << R2  << endl; 
  }

  void Cloop(array1<unsigned>::opt &C, Ca afp, Ca dfp, Cm mfp) {
    for(unsigned I=Astart; I < Astop; ++I)
      (parent->*afp)(C,I);
    for(unsigned I=Dstart; I < Dstop; ++I)
      (parent->*dfp)(C,I);
    for(unsigned I=1; I < m; I++) {
      unsigned start =jstart(I);
      unsigned stop =jstop(I);
      for(unsigned j=start; j < stop; ++j) {
	(parent->*mfp)(C,I,j);
      }
    }
    //cout << "C:" << C << endl;
  }

  void Invariantsloop(const array2<Complex> w, Real &E, Real &Z, Real &P,
	       Invarfp fp) {
    vector wx=w[xorigin];
    for(unsigned I=Astart; I < Astop; ++I) {
      vector wxi=w[xorigin+I];
      (parent->*fp)(I*I,abs2(wx[I])+abs2(wxi[0]),E,Z,P);
    }
    for(unsigned I=Dstart; I < Dstop; ++I) {
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      (parent->*fp)(2*I*I,abs2(wi[I])+abs2(wim[I]),E,Z,P);
    }
    for(unsigned I=1; I < m; I++) {
      unsigned I2=I*I;
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      unsigned start =jstart(I);
      unsigned stop =jstop(I);
      for(unsigned j=start; j < stop; ++j) {
	(parent->*fp)(I2+j*j,
		      abs2(wi[j])+abs2(wim[j])+
		      abs2(w[xorigin-j][I])+abs2(w[xorigin+j][I]),
		      E,Z,P);
      }
    }
  }

  // loop over the Hermitian-symmetric array2 w, calculate
  // something, put it into the appropriate index of S, a spectrum-like array
  void Sloop(vector& S, const array2<Complex> w,Sa  afp,Sa  dfp, Sm mfp) {
    vector wx=w[xorigin];
    for(unsigned I=Astart; I < Astop; ++I) {
      vector wxi=w[xorigin+I];
      (parent->*afp)(I,S,wx[I],wxi[0]);
    }
    for(unsigned I=Dstart; I < Dstop; ++I) {
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      (parent->*dfp)(I,S,wi[I],wim[I]);
    }
    for(unsigned I=1; I < m; I++) {
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      unsigned start =jstart(I);
      unsigned stop =jstop(I);
      for(unsigned j=start; j < stop; ++j)
	(parent->*mfp)(I,j,S,wi[j],wim[j],w[xorigin-j][I],w[xorigin+j][I]);
    }
  }
  
  void Tloop(vector& T, const array2<Complex> w, const array2<Complex> f,
	     Ta afp, Ta dfp, Tm mfp) {
    vector wx=w[xorigin];
    vector fx=f[xorigin];
    for(unsigned I=Astart; I < Astop; ++I) {
      vector wxi=w[xorigin+I];
      vector fxi=f[xorigin+I];
      (parent->*afp)(I,T,wx[I],wxi[0],fx[I],fxi[0]);
    }
    for(unsigned I=Dstart; I < Dstop; ++I) {
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      vector fi=f[i];
      vector fim=f[im];
      (parent->*dfp)(I,T,wi[I],wim[I],fi[I],fim[I]);
    }

    for(unsigned I=1; I < m; I++) {
      unsigned i=xorigin+I;
      unsigned im=xorigin-I;
      vector wi=w[i];
      vector wim=w[im];
      vector fi=f[i];
      vector fim=f[im];
      unsigned start =jstart(I);
      unsigned stop =jstop(I);
      for(unsigned j=start; j < stop; ++j) {
	(parent->*mfp)(I,j,T,
		       wi[j],wim[j],w[xorigin-j][I],w[xorigin+j][I],
		       fi[j],fim[j],f[xorigin-j][I],f[xorigin+j][I]);
      }
    }
  }

  // FIXME: something is wrong here
  void killmodes(array2<Complex> A) {
    if(circular) {
      vector Ax=A[xorigin];
      for(unsigned I=Astop; I < m; ++I) 
	Ax[I]=A[xorigin+I][0]=0;
      for(unsigned I=Dstop; I < m; ++I) 
	A[xorigin+I][I]=A[xorigin-I][I]=0;
      for(unsigned I=Dstop; I < m; I++) {
	vector Ai=A[xorigin+I];
	vector Aim=A[xorigin-I];
	unsigned start =jkillstart(I);
	//cout << I << "->" << start << endl;
	for(unsigned j=start; j < m; ++j) {
	  //cout << " " << j ;
	  Ai[j]=Aim[j]=A[xorigin-j][I]=A[xorigin+j][I]=0;
	}
	//cout << endl;
      }
    }
    // TODO:
    // radix-2 and radix-4 kills as well?
    // via isvisible(unsigned I, unsigned j) ?
  }

};


#endif
