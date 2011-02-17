#include "dnsbase.h"

const int DNSBase::xpad=1;
const int DNSBase::ypad=1;

//***** Source routines *****//

void DNSBase::NonLinearSource(const vector2& Src, const vector2& Y, double t)
{
  w.Set(Y[OMEGA]);
  f0.Set(Src[OMEGA]);
  
  f0(origin)=0.0;
  f1(origin)=0.0;
  g0(origin)=0.0;
  g1(origin)=0.0;

  for(unsigned i=0; i < Nx; ++i) {
    int I=(int) i-(int) xorigin;
    Real kx=k0*I;
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

vector DNSBase::E;
vector DNSBase::T;
Real DNSBase::k0;
Real DNSBase::k02;
Real DNSBase::nuH;
Real DNSBase::nuL;
Real DNSBase::etanorm;

void DNSBase::FETL(const vector& wi, const vector& Si, unsigned I2, unsigned j)
{
  unsigned k2int=I2+j*j;
  Real kint=sqrt(k2int);
  Real k=k0*kint;
  unsigned index=(unsigned)(kint-0.5);
  Complex wij=wi[j];
  Real w2=abs2(wij);
  Nu nu=nuk(k2int);
  E[index] += Complex(w2/k,nu*w2);
  Complex Sij=Si[j];
  Complex& Tindex=T[index];
  Tindex.re += realproduct(Sij,wij);
  Forcing->Force(wij,Tindex.im,k);
  Si[j]=Sij-nu*wij;
}

void DNSBase::FTL(const vector& wi, const vector& Si, unsigned I2, unsigned j)
{
  unsigned k2int=I2+j*j;
  Real kint=sqrt(k2int);
  Real k=k0*kint;
  unsigned index=(unsigned)(kint-0.5);
  Complex wij=wi[j];
  Nu nu=nuk(k2int);
  Complex Sij=Si[j];
  Complex& Tindex=T[index];
  Tindex.re += realproduct(Sij,wij);
  Forcing->Force(wij,Tindex.im,k);
  Si[j]=Sij-nu*wij;
}

void DNSBase::FET(const vector& wi, const vector& Si, unsigned I2, unsigned j)
{
  unsigned k2int=I2+j*j;
  Real kint=sqrt(k2int);
  Real k=k0*kint;
  unsigned index=(unsigned)(kint-0.5);
  Complex wij=wi[j];
  Real w2=abs2(wij);
  Nu nu=nuk(k2int);
  E[index] += Complex(w2/k,nu*w2);
  Complex Sij=Si[j];
  Complex& Tindex=T[index];
  Tindex.re += realproduct(Sij,wij);
  Forcing->Force(wij,Tindex.im,k);
}

void DNSBase::FE(const vector& wi, const vector& Si,unsigned I2, unsigned j)
{
  unsigned k2int=I2+j*j;
  Real kint=sqrt(k2int);
  Real k=k0*kint;
  unsigned index=(unsigned)(kint-0.5);
  Complex wij=wi[j];
  Real w2=abs2(wij);
  Nu nu=nuk(k2int);
  E[index] += Complex(w2/k,nu*w2);
}

void DNSBase::FL(const vector& wi, const vector& Si, unsigned I2, unsigned j)
{
  unsigned k2int=I2+j*j;
  Complex wij=wi[j];
  Nu nu=nuk(k2int);
  Complex Sij=Si[j];
  double T;
  Real k=k0*sqrt(k2int);
  Forcing->Force(wij,T,k);
  Si[j]=Sij-nu*wij;
}

void DNSBase::Compute(SourceFcn fcn, const vector2& Src, const vector2& Y)
{
  f0.Set(Src[OMEGA]);
  w.Set(Y[OMEGA]);

  Forcing->Set();

  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    vector Si=f0[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j)
      (*fcn)(wi,Si,I2,j);
  }
}

void DNSBase::Stochastic(const vector2&Y, double, double dt)
{
  w.Set(Y[OMEGA]);
  Set(T,Y[TRANSFER]);
  
  Forcing->SetStochastic(dt);
  
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real kint=sqrt(I2+j*j);
      unsigned index=(unsigned)(kint-0.5);
      Forcing->ForceStochastic(wi[j],T[index].im,k0*kint);
    }
  }
}

void DNSBase::Initialize()
{
  fevt << "# t\tE\tZ\tP" << endl;
  setcount();
}

void DNSBase::setcount()
{
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

void DNSBase::ComputeInvariants(const array2<Complex> &w, 
				Real& E, Real& Z, Real& P)
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

