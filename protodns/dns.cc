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
  // This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
  // requires only 4 FFTs per stage.
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k2inv=1.0/(i*i+j*j);
      double jk2inv=j*k2inv;
      double ik2inv=i*k2inv;
      f0[i][j]=Complex(-w[i][j].im*jk2inv,w[i][j].re*jk2inv); // u
      f1[i][j]=Complex(w[i][j].im*ik2inv,-w[i][j].re*ik2inv); // v
      f2[i][j] = Complex(i*f0[i][j].im,i*f0[i][j].re); //dvdy
      f3[i][j] = Complex(i*f1[i][j].im,i*f1[i][j].re)+Complex(j*f2[i][j].im,j*f2[i][j].re); //dvdx + dudy
    }
  }
  Backward.fftNormalized(f0);
  Backward.fftNormalized(f1);
  Backward.fftNormalized(f2);
  Backward.fftNormalized(f3);
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
        f2[i][j]= f1[i][j]*f1[i][j] - f0[i][j]*f0[i][j] - 4*(Cs*delta)*(Cs*delta)*sqrt((f2[i][j]*f2[i][j])+(f3[i][j])*(f3[i][j]))*f2[i][j]; // B(x,t)
        f0[i][j]= f0[i][j]*f1[i][j]; // u1*u2
    }
  }

  Forward.fft(f0);
  Forward.fft(f2);

  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      f0[i][j]=i*j*f0[i][j]+(i*i-j*j)*f2[i][j]-nu*(i*i+j*j)*w[i][j];
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

  Convolution=new ImplicitHConvolution2(mx,my,true,true,2,2);

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
