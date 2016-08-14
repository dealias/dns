#include "dnsbase.h"

const int DNSBase::xpad=0;
const int DNSBase::ypad=0;

void DNSBase::NonLinearSource(const vector2& Src, const vector2& Y, double t)
{
  f0.Dimension(Nx+1,my,-imx,0);
  
  w.Set(Y[OMEGA]);
  f0.Set(Src[PAD]);

  f0[0][0]=0.0;
  f1[0][0]=0.0;
  
  // This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
  // requires only 4 FFTs per stage.
#pragma omp parallel for num_threads(threads)
  for(int i=-imx+1; i < imx; ++i) {
    Real kx=k0*i;
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    rvector k2invi=k2inv[i];
    for(unsigned j=i <= 0 ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      Complex wij=wi[j];
      Real k2invij=k2invi[j];
      Real kyk2inv=ky*k2invij;
      Real kxk2inv=kx*k2invij;
      f0i[j]=Complex(-wij.im*kyk2inv,wij.re*kyk2inv); // u
      f1i[j]=Complex(wij.im*kxk2inv,-wij.re*kxk2inv); // v
    }
  }

  F[0]=f0;
  Convolution->convolve(F,multadvection2);
  f0[0][0]=0.0;
  
  for(int i=-imx+1; i < imx; ++i) {
    Real kx=k0*i;
    Real kx2=kx*kx;
    vector f0i=f0[i];
    vector f1i=f1[i];
    for(unsigned j=i <= 0 ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      Real ky2=ky*ky;
      f0i[j]=kx*ky*f0i[j]+(kx2-ky2)*f1i[j];
    }
  }
  
#if 0
  Real sum=0.0;
  for(int i=-imx+1; i < imx; ++i) {
    Real kx=k0*i;
    vector wi=w[i];
    for(unsigned j=i <= 0 ? 1 : 0; j < my; ++j) {
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
  
  k2inv.Allocate(Nx,my,-imx+1,0);
  for(int i=-imx+1; i < imx; ++i) {
    Real kx=k0*i;
    Real kx2=kx*kx;
    rvector k2invi=k2inv[i];
    for(unsigned j=i <= 0 ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      k2invi[j]=1.0/(kx2+ky*ky);
    }
  }
}
  
void DNSBase::InitialConditions()
{
  w(0,0)=0;
  
  Loop(Initw(this),InitializeValue(this));
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,w);
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
    
  for(int i=-imx+1; i < imx; ++i)
    for(unsigned int j=0; j < my; ++j)
      f1[i][j]=w(i,j);
  
  Backward->fft0(f1);

  fw << 1 << 2*my-1 << Nx+1;
  for(int j=2*my-2; j >= 0; j--) {
    for(unsigned int i=0; i <= Nx; i++) {
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

