#include <cmath>
#include <iomanip>
#include <fstream>
#include "Complex.h"
#include "convolution.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;


int Nx=255; // Number of modes in x direction
int Ny=255; // Number of modes in y direction

double Cs = 0.1;
double delta = 1/Nx;

double dt=1.0e-6;
double nu=0.003; // kinematic viscosity

int mx;
int my;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector2 w;
//added new vectors
vector2 f0,f1,f2,f3;

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution2 *Convolution;

void init(vector2& w)
{
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      w[i][j]=1.0/(i*i+j*j);
    }
  }
}

void Source(const vector2& w, vector2 &S)
{
  f0[0][0]=0.0; // Enforce no mean flow.
  f1[0][0]=0.0;
  fft2d Forward(Nx,my,-1);
  fft2d Backward(Nx,my,1);
  // This 2D version requires only 6 FFTs per stage (in the spirit
  // of Basdevant, J. Comp. Phys, 50, 1983).
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k2inv=1.0/(i*i+j*j);
      double jk2inv=j*k2inv;
      double ik2inv=i*k2inv;
      Complex wij=w[i][j];
      double u=Complex(-wij.im*jk2inv,wij.re*jk2inv);
      double v=Complex(wij.im*ik2inv,-wij.re*ik2inv);
      f0[i][j]=u;
      f1[i][j]=v;
      f2[i][j]=Complex(-i*u.im,i*u.re); // dudx
      f3[i][j]=Complex(-i*v.im,i*v.re)+Complex(-j*u.im,j*u.re); // dvdx + dudy
    }
  }
  Backward.fftNormalized(f0);
  Backward.fftNormalized(f1);
  Backward.fftNormalized(f2);
  Backward.fftNormalized(f3);
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double v=f1[i][j];
      double u=f0[i][j];
      double C=Cs*delta;
      double ux=f2[i][j];
      double s12=f3[i][j];
      f1[i][j]=v*v-u*u+4*C*C*sqrt(ux*ux+s12*s12)*ux; // B(x,t)
      f0[i][j]=u*v;
    }
  }

  Forward.fft(f0);
  Forward.fft(f1);

  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      f0[i][j]=(i*i-j*j)*f0[i][j]+i*j*f1[i][j]-nu*(i*i+j*j)*w[i][j];
    }
  }
}

void Spectrum()
{
  ofstream zkvk("zkvk",ios::out);

  int kmax=(int) hypot(mx-1,my-1);
  double Z[kmax+1];
  for(int k=0; k <= kmax; ++k) Z[k]=0.0;
     
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      int k=sqrt(i*i+j*j);
      Z[(int) (k+0.5)] += abs2(w[i][j]);
    }
  }
  zkvk << "# k\tZ(k)" << endl;
  
  for(int k=1; k <= kmax; ++k) {
    zkvk << k << "\t" << Z[k] << endl;
  }
}

void Output(int step, bool verbose=false)
{
  double E=0.0, Z=0.0, P=0.0;
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double w2=abs2(w[i][j]);
      double k2=i*i+j*j;
      P += k2*w2;
      Z += w2;
      E += w2/k2;
    }
  }
  if(verbose) {
    cout << "t=" << step*dt << endl;
    cout << "Energy=" << E << endl;
    cout << "Enstrophy=" << Z << endl;
    cout << "Palenstrophy=" << P << endl;
  }
  ezvt << E << "\t" << Z << "\t" << P << endl;
}

int main(int argc, char* argv[])
{
  int n;
  cout << "Number of time steps? " << endl;
  cin >> n;

  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  size_t align=sizeof(Complex);

  f0.Allocate(Nx,my,-mx+1,0,align);
  f1.Allocate(Nx,my,-mx+1,0,align);
  f2.Allocate(Nx,my,-mx+1,0,align);
  f3.Allocate(Nx,my,-mx+1,0,align);

  Convolution=new ImplicitHConvolution2(mx,my,true,true,4,2);

  w.Allocate(Nx,my,-mx+1,0,align);

  init(w);
  w(0,0)=0.0; // Enforce no mean flow.

  cout.precision(15);

  for(int step=0; step < n; ++step) {
    Output(step,step == 0);
     Source(w,f0);
     for(int i=-mx+1; i < mx; ++i) {
       for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
	 w[i][j] += f0[i][j]*dt;
       }
     }
//     cout << "[" << step << "] ";
  }
  cout << endl;
  Output(n,true);
  Spectrum();

  return 0;
}
