#include "dnsbase.h"

const int DNSBase::xpad=0;
const int DNSBase::ypad=0;

void DNSBase::NonLinearSource(const vector2& Src, const vector2& Y, double t)
{
  f0.Dimension(Nx+1,my,-1,0);
  
  w.Set(Y[OMEGA]);
  f0.Set(Src[PAD]);

  for(unsigned int j=0; j < my; ++j) {
    f0(j)=0; // Required?
    f1(j)=0;
    g0(j)=0;
    g1(j)=0;
  }
  
  f0[xorigin][0]=0.0;
  f1[xorigin][0]=0.0;
  g0[xorigin][0]=0.0;
  g1[xorigin][0]=0.0;
  
  int imx=(int) mx;
#pragma omp parallel for num_threads(threads)
  for(int I=-imx+1; I < imx; ++I) {
    Real kx=k0*I;
    unsigned i=I+xorigin;
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    vector g0i=g0[i];
    vector g1i=g1[i];
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

  F[0]=f0;
  Convolution->convolve(F,multbinary2);
  f0[xorigin][0]=0.0;
  
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,&f0[0][0]);
  
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
  for(unsigned int j=0; j < my; ++j)
    f1(j)=0;
    
  for(unsigned int i=0; i < Nx; ++i)
    for(unsigned int j=0; j < my; ++j)
      f1[i][j]=w(i,j);
  
  Backward->fft0(f1,wr);

  fw << 1 << 2*my-1 << Nx+1;
  for(int j=2*my-2; j >= 0; j--) {
    for(unsigned i=0; i <= Nx; i++) {
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

