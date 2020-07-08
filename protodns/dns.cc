#include <cmath>
#include <iomanip>
#include <fstream>
#include "Complex.h"
#include "convolution.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

int Nx=1023; // Number of modes in x direction
int Ny=1023; // Number of modes in y direction

double dt=1.0e-6;
double nu=0.0; // kinematic viscosity

int mx;
int my;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector2 w;
vector2 f0,f1;

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution2 *Convolution;

void init(vector2& w)
{
  for(int i=-mx+1; i < mx; ++i) {
    vector wi=w[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      wi[j]=1.0/(i2+j*j);
    }
  }
}
    
void Source(const vector2& w, vector2 &S)
{
  f0[0][0]=0.0; // Enforce no mean flow.
  f1[0][0]=0.0;
  
  // This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
  // requires only 4 FFTs per stage.
  for(int i=-mx+1; i < mx; ++i) {
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k2inv=1.0/(i2+j*j);
      double jk2inv=j*k2inv;
      double ik2inv=i*k2inv;
      Complex wij=wi[j];
      f0i[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv); // u_k
      f1i[j]=Complex(wij.im*ik2inv,-wij.re*ik2inv); // v_k
    }
  }

  Complex *F[]={f0,f1};

  // u_k, v_k -> F{v^2-u^2}_k, F{u*v}_k
  Convolution->convolve(F,multadvection2);
  
  for(int i=-mx+1; i < mx; ++i) {
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      int j2=j*j;
      f0i[j]=i*j*f0i[j]+(i2-j2)*f1i[j]-nu*(i2+j2)*wi[j];
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
    vector wi=w[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      int k=sqrt(i2+j*j);
      Z[(int) (k+0.5)] += abs2(wi[j]);
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
    vector wi=w[i];
    double i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double w2=abs2(wi[j]);
      double k2=i2+j*j;
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

  Convolution=new ImplicitHConvolution2(mx,my,true,true,2,2);
  
  w.Allocate(Nx,my,-mx+1,0,align);
  
  init(w);
  w[0][0]=0.0; // Enforce no mean flow.

  cout.precision(15);
  
  for(int step=0; step < n; ++step) {
    Output(step,step == 0);
     Source(w,f0);
     for(int i=-mx+1; i < mx; ++i) {
       vector wi=w[i];
       vector f0i=f0[i];
       for(int j=(i <= 0 ? 1 : 0); j < my; ++j)
	 wi[j] += f0i[j]*dt;
     }
     cout << "[" << step << "] " << flush;
  }
  cout << endl;
  Output(n,true);
  Spectrum();
     
  return 0;
}
