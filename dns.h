#ifndef __dns_h__
#define __dns_h__ 1

#include "Array.h"
#include "kernel.h"
#include "fftw++.h"
#include "convolution.h"
#include "Forcing.h"
#include "InitialCondition.h"

using namespace Array;
//using namespace fftwpp;
using std::ostringstream;

extern const char *ic;
//extern char *linearity="Power";
extern const char *forcing;

// Vocabulary
extern Real nuH, nuL;
extern int pH;
extern int pL;
extern unsigned Nx;
extern unsigned Ny;
extern Real eta;
extern Complex force;
extern Real kforce;
extern Real deltaf;
extern unsigned movie;
extern unsigned rezero;
extern unsigned spectrum;
extern Real icalpha;
extern Real icbeta;

extern int xpad;
extern int ypad;


// TODO: initial conditions go here

// TODO: forcing goes here
// TODO: once forcing, then move DNS::Transfer and DNS::Stochastic

// TODO: after ic and force done, vocab goes here

// TODO: After vocab is moved, then DNS:InitialConditions, etc.

//extern DNSVocabulary DNS_Vocabulary;

extern InitialConditionBase *InitialCondition;
extern ForcingBase *Forcing;



class DNS : public ProblemBase {
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
  Complex *block;
  Complex *F[2];
  Complex *G[2];
  
  fftwpp::ImplicitHConvolution2 *Convolution;
  fftwpp::ExplicitHConvolution2 *Padded;
  
  ifstream ftin;
  oxstream fwk,fw,fekvk,ftransfer;
  ofstream ft,fevt;

public:
  DNS();
  virtual ~DNS();
  
  void IndexLimits(unsigned& start, unsigned& stop,
		   unsigned& startT, unsigned& stopT,
		   unsigned& startM, unsigned& stopM) {
    start=Start(OMEGA);
    stop=Stop(OMEGA);
    startT=Start(TRANSFER);
    stopT=Stop(TRANSFER);
    startM=Start(EK);
    stopM=Stop(EK);
  }
  unsigned getNx() {return Nx;}
  unsigned getmx() {return mx;}
  unsigned getmy() {return my;}
  Real getk0() {return k0;}
  Real getk02() {return k02;}
  unsigned getxorigin() {return xorigin;}


  void InitialConditions();
  void Initialize();
  void Output(int it);
  void FinalOutput();
  void OutFrame(int it);
  
  void Spectrum(vector& S, const vector& y);
  void Transfer(const vector2& Src, const vector2& Y);
  
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  
  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) Transfer(Src,Y);
    LinearSource(Src,Y,t);
  }
  
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum) Spectrum(Src[EK],Y[OMEGA]);
    fftwpp::HermitianSymmetrizeX(mx,my,xorigin,Src[OMEGA]);
  }
  
  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
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

  void ComputeInvariants(Real& E, Real& Z, Real& P);
  void Stochastic(const vector2& Y, double, double);
  
  Real Spectrum(unsigned int i) {
    return T[i].re*twopi/count[i];
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
extern DNS *DNSProblem;



void DNS::LinearSource(const vector2& Src, const vector2& Y, double)
{
  w.Set(Y[OMEGA]);
  f0.Set(Src[OMEGA]);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector f0i=f0[i];
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j)
      f0i[j] -= nuk(k02*(I2+j*j))*wi[j];
  }
}


void DNS::Initialize()
{
  fevt << "#   t\t\t E\t\t\t Z" << endl;
  
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


void DNS::NonLinearSource(const vector2& Src, const vector2& Y, double)
{
  w.Set(Y[OMEGA]);
  f0.Set(Src[OMEGA]);
 
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
  Convolution->convolve(F,G);
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


void DNS::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

class cwrap{
public:
  static Real Spectrum(unsigned int i)
  {
    return DNSProblem->Spectrum(i);
  }

  static Real Dissipation(unsigned int i)
  {
    return DNSProblem->Dissipation(i);
  }
  
  static Real Pi(unsigned int i)
  {
    return DNSProblem->Pi(i);
  }
  
  static Real Eta(unsigned int i)
  {
    return DNSProblem->Eta(i);
  }
};

void DNS::Output(int it)
{
  Real E,Z,P;
	
  w.Set(y);
  ComputeInvariants(E,Z,P);
  fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;

  Complex *y=Y[0];
  if(output) out_curve(fw,y,"w",NY[0]);
  
  if(movie) OutFrame(it);
	
  if(spectrum) {
    ostringstream buf;
    Set(T,Y[EK]);
    buf << "ekvk" << dirsep << "t" << tcount;
    open_output(fekvk,dirsep,buf.str().c_str(),0);
    out_curve(fekvk,t,"t");
    out_curve(fekvk,cwrap::Spectrum,"Ek",nshells);
    out_curve(fekvk,cwrap::Dissipation,"nuk*Ek",nshells);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");

    Set(T,Y[TRANSFER]);
    buf.str("");
    buf << "transfer" << dirsep << "t" << tcount;
    open_output(ftransfer,dirsep,buf.str().c_str(),0);
    out_curve(ftransfer,t,"t");
    out_curve(ftransfer,cwrap::Pi,"Pi",nshells);
    out_curve(ftransfer,cwrap::Eta,"Eta",nshells);
    ftransfer.close();
    if(!ftransfer) msg(ERROR,"Cannot write to file transfer");
  }    

  tcount++;
  ft << t << endl;
  
  if(rezero && it % rezero == 0 && spectrum) {
    vector2 Y=Integrator->YVector();
    vector T=Y[TRANSFER];
    for(unsigned i=0; i < nshells; i++)
      T[i]=0.0;
    vector S=Y[EK];
    for(unsigned i=0; i < nshells; i++)
      S[i]=0.0;
  }
}




#endif
