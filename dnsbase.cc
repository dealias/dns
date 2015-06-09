#include "dnsbase.h"

const int DNSBase::xpad=1;
const int DNSBase::ypad=1;

void DNSBase::NonLinearSource(const vector2& Src, const vector2& Y, double t)
{
  w.Set(Y[OMEGA]);
  f0.Set(Src[OMEGA]);

  F0(origin0)=0.0;
  f1(origin0)=0.0;
  g0(origin0)=0.0;
  g1(origin0)=0.0;
  
  int imx=(int) mx;
#pragma omp parallel for num_threads(threads)
  for(int I=-imx+1; I < imx; ++I) {
    Real kx=k0*I;
    unsigned i=I+xorigin;
    unsigned i0=I+xorigin0;
    vector wi=w[i];
    vector f0i=F0[i0];
    vector f1i=f1[i0];
    vector g0i=g0[i0];
    vector g1i=g1[i0];
    rvector k2invi=k2inv[i];
    for(unsigned j=I <= 0 ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      Complex wij=wi[j];
      Complex kxw=Complex(-kx*wij.im,kx*wij.re);
      Complex kyw=Complex(-ky*wij.im,ky*wij.re);
      f0i[j]=kxw;
      f1i[j]=kyw;
      Real k2invij=k2invi[j];
      g0i[j]=k2invij*kyw;
      g1i[j]=-k2invij*kxw;
    }
  }

//  F[0]=f0;
//  Convolution->convolve(F,multbinary2);
  Convolution->convolve(F,F+2);
  
  for(int I=-imx+1; I < imx; ++I) {
    unsigned i=I+xorigin;
    unsigned i0=I+xorigin0;
    for(unsigned j=I <= 0 ? 1 : 0; j < my; ++j) {
      f0[i][j]=F0[i0][j];
    }
  }
  
  f0(origin)=0.0;
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,f0);
  
#if 0
  Real sum=0.0;
  for(int I=-imx+1; I < imx; ++I) {
    Real kx=k0*I;
    unsigned i=I+xorigin;
    vector wi=w[i];
    for(unsigned j=I <= 0 ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      Complex wij=wi[j];
//      sum += (f0[i][j]*conj(wij)).re;
      sum += (f0[i][j]*conj(wij)/(kx*kx+ky*ky)).re;
    }
  }
  cout << sum << endl;
  cout << endl;
#endif
}

void DNSBase::Initialize()
{
  fevt << "# t\tE\tZ\tP" << endl;
}

void DNSBase::SetParameters()
{
  setcount();
  fcount=0;
  
  Forcing->Init();
  
  Loop(InitNone(this),ForcingCount(this));

  fcount *= 2; // Account for Hermitian conjugate modes.

  Forcing->Init(fcount);
  
  k2inv.Allocate(Nx,my);
  int imx=(int) mx;
  for(int I=-imx+1; I < imx; ++I) {
    Real kx=k0*I;
    Real kx2=kx*kx;
    rvector k2invi=k2inv[I+xorigin];
    for(unsigned j=I <= 0 ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      k2invi[j]=1.0/(kx2+ky*ky);
    }
  }
}
  
void DNSBase::InitialConditions()
{
  w(xorigin,0)=0;
  
  Loop(Initw(this),InitializeValue(this));
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);
}

void DNSBase::setcount()
{
#pragma omp parallel for num_threads(threads)
  for(unsigned i=0; i < nshells; i++)
    count[i]=0;
  
  if(spectrum)
      Loop(InitNone(this),Count(this));
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
  Energy=Enstrophy=Palenstrophy=0.0;
  
  Loop(Initw(this),Invariants(this));
  
  E=Energy;
  Z=Enstrophy;
  P=Palenstrophy;
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

