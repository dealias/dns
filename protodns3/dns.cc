#include <cmath>
#include <iomanip>
#include <fstream>
#include "Complex.h"
#include "convolution.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

int Nx=15; // Number of modes in x direction
int Ny=15; // Number of modes in y direction
int Ny=15; // Number of modes in z direction

double dt=1.0e-3;
double nu=0.0; // kinematic viscosity

int mx;
int my;
int mz;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector3 u,v,w;
vector3 f0,f1,g0,g1;

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution2 *Convolution;

void init(vector3& u, vector3& v, vector3& w)
{
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int l=-my+1; l < my; ++l) {
	u[i][j][l]=1.0/(i*i+j*j);
	v[i][j][l]=1.0/(i*i+j*j);
        w[i][j][l]=1.0/(i*i+j*j);
    }
  }
}
    
  void Source_component(const vector3& u, const vector3& v, const vector3& w, char component, vector3 &S){
    f0[0][0]=0.0;
    f1[0][0]=0.0;
    g0[0][0]=0.0;
    g1[0][0]=0.0;
  
    switch (char component) {
    case 'u':
      for(int i=-mx+1; i < mx; ++i) {
	for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
	  for(int l=((i || j) <=0?1:0); l < mz; ++l{
	      Complex Ikxu=Complex(-i*u[i][j][l].im,i*u[i][j][l].re);
	      Complex Ikyu=Complex(-j*u[i][j][l].im,j*u[i][j][l].re);
	      Complex Ikzu=Complex(-l*u[i][j][l].im,l*u[i][j][l].re);
	      f0[i][j][l]=Ikxu;
	      f1[i][j][l]=Ikyu;
	      f2[i][j][l]=Ikzu;
	      
	      g0[i][j][l]=u[i][j][l];
	      g1[i][j][l]=v[i][j][l];
	      g3[i][j][l]=w[i][j][l];
	  }
	}
      }
      Complex *F[]={f0,f1,f2,g0,g1,g3};
      Convolution->convolve(F,multbinary2);
      for(int i=0; i<mx ;++i){
	for(int j=0; j<my ;++j){
	  for(int l=0; l<mz ;++l){
	     double kx2=i*i;
	     f0[i][j][l] += -nu*i*i*u[i][j][l];
	  }
	}
      }
    break;
  case 'v':
    ....
    break;
  case 'w':
    ....
    break;
  default:
    cout << "Character does not match any of the components!" < <endl;
  }


      
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
  
  for(int i=0; i<mx ;++i){
    for(int j=0; j<my ;++j){
       double k2=i*i+j*j;     
       f0[i][j] += -nu*k2*w[i][j];
    }
  }
  f0(0,0)=0.0; // Enforce no mean flow.
  
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,f0);
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
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  mz=(Nz+1)/2;
  size_t align=sizeof(Complex);
  
  f0.Allocate(Nx,my,-mx+1,0,align);
  f1.Allocate(Nx,my,-mx+1,0,align);
  
  g0.Allocate(Nx,my,-mx+1,0,align);
  g1.Allocate(Nx,my,-mx+1,0,align);

  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,true,true,4,1);
  
  w.Allocate(Nx,my,-mx+1,0,align);
  
  init(w);
  w(0,0)=0.0; // Enforce no mean flow.
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,w);

  int n=10000;
  
  cout.precision(15);
  
  for(int step=0; step < n; ++step) {
    Output(step,step == 0);
     Source(w,f0);
     for(int i=0; i < mx; ++i) {
       for(int j=0; j < my; ++j) {
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
