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

double dt=1.0e-8;
double nu=0.0; // kinematic viscosity

int mx;
int my;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector2 w;
vector2 f0,f1,g0,g1;

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
  
  for(int i=0; i < mx ;++i){
    for(int j=0; j < my ;++j){
       double k2=i*i+j*j;     
       f0[i][j] += -nu*k2*w[i][j];
    }
  }
  f0(0,0)=0.0; // Enforce no mean flow.
  
  HermitianSymmetrizeX(mx,my,mx-1,f0);
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
  double E=0.0, Z=0.0;
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double w2=abs2(w[i][j]);
      double k2=i*i+j*j;
      Z += w2;
      E += w2/k2;
    }
  }
  if(verbose) {
    cout << "t=" << step*dt << endl;
    cout << "Energy=" << E << endl;
    cout << "Enstrophy=" << Z << endl;
  }
  ezvt << E << "\t" << Z << endl;
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
  g0.Allocate(Nx,my,-mx+1,0,align);
  g1.Allocate(Nx,my,-mx+1,0,align);

  Convolution=new ImplicitHConvolution2(mx,my,true,true,4,1);
  
  w.Allocate(Nx,my,-mx+1,0,align);
  
  init(w);
  w(0,0)=0.0; // Enforce no mean flow.
  HermitianSymmetrizeX(mx,my,mx-1,w);

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
