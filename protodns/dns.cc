#include "Complex.h"
#include "convolution.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

// size of problem

int Nx=15;
int Ny=15;

int mx;
int my;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector2 w;
vector2 f0,f1,g0,g1;

ImplicitHConvolution2 *Convolution;

void init(vector2& w)
{
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      w[i][j]=Complex(i,j);
    }
  }
}
    
void Source(const vector2& w, vector2 &S)
{
  f0[0][0]=0.0;
  f1[0][0]=0.0;
  g0[0][0]=0.0;
  g1[0][0]=0.0;
  
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k2=i*i+j*j;
      Complex Ikxw=Complex(-i*w[i][j].im,i*w[i][j].re);
      Complex Ikyw=Complex(-j*w[i][j].im,j*w[i][j].re);
      f0[i][j]=Ikxw;
      f1[i][j]=Ikyw;
      double k2inv=1.0/k2;
      g0[i][j]=-Ikyw*k2inv;
      g1[i][j]=Ikxw*k2inv;
    }
  }

  Complex *F[]={f0,f1,g0,g1};
  Convolution->convolve(F,multbinary2);
  f0[0][0]=0.0; // Enforce no mean flow.
  
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,f0);
  
#if 1 // Check enstrophy and energy symmetries
  double sumE=0.0;
  double sumZ=0.0;
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      sumZ += (f0[i][j]*conj(w[i][j])).re;
      sumE += (f0[i][j]*conj(w[i][j])/(i*i+j*j)).re;
    }
  }
  cout << sumZ << endl;
  cout << sumE << endl;
  cout << endl;
#endif
}

int main(int argc, char* argv[])
{
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
    
  size_t align=sizeof(Complex);
  
  f0.Allocate(Nx,my,-mx+1,0,align);
  f1.Allocate(Nx,my,-mx+1,0,align);
  g0.Allocate(Nx,my,-mx+1,0,align);
  g1.Allocate(Nx,my,-mx+1,0,align);

  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,true,true,4,1);
  
  w.Allocate(Nx,my,-mx+1,0,align);
  
  init(w);
  
  w(0,0)=0;
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,w);

  Source(w,f0);
    
  return 0;
}
