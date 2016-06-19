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
typedef Array3<Complex> vector3;

vector3 u,v,w;
vector3 f0,f1,g0,g1;

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution2 *Convolution;

void init(vector3& u, vector3& v, vector3& w)
{
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=-my+1; k < my; ++k) {
        Complex val=1.0/(i*i+j*j+k*k);
	u[i][j][k]=val;
	v[i][j][k]=val;
        w[i][j][k]=val;
      }
    }
  }
}  

void multadvection3(double **F, unsigned int m,
                    const unsigned int indexsize,
                    const unsigned int *index,
                    unsigned int r, unsigned int threads)
{
  double* F0=F[0];
  double* F1=F[1];
  double* F2=F[2];
  double* F3=F[3];
  double* F4=F[4];
  double* F5=F[5];
  double* F6=F[6];
  double* F7=F[7];
  double* F8=F[8];
  double* F9=F[9];
  double* F10=F[10];
  
  for(unsigned int j=0; j < m; ++j) {
    // I*(u*kx+v*ky+w*kz)*u
    // I*(u*kx+v*ky+w*kz)*v
    // I*(u*kx+v*ky+w*kz)*w
    F0[j]=F0[j]*F3[j]+F1[j]*F4[j]+F2[j]*F5[j];
    F1[j]=F0[j]*F6[j]+F1[j]*F7[j]+F2[j]*F8[j];
    F2[j]=F0[j]*F9[j]+F1[j]*F10[j]-F2[j]*(F3[j]+F7[j]);
  }
}

void Source(const vector4& u, vector4 &S)
{
  f0[0][0][0]=0.0;
  f1[0][0][0]=0.0;
  f2[0][0][0]=0.0;
  f3[0][0][0]=0.0;
  f4[0][0][0]=0.0;
  f5[0][0][0]=0.0;
  f6[0][0][0]=0.0;
  f7[0][0][0]=0.0;
  f8[0][0][0]=0.0;
  f9[0][0][0]=0.0;
  f10[0][0][0]=0.0;
  
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        Complex U=u[0][i][j][k];
        Complex V=u[1][i][j][k];
        Complex W=u[2][i][j][k];
        
        f0[i][j][k]=U;
        f1[i][j][k]=V;
        f2[i][j][k]=W;

        f3[i][j][k]=Complex(-i*U.im,i*U.re);  // I*kx*U
        f4[i][j][k]=Complex(-j*U.im,j*U.re);  // I*ky*U
        f5[i][j][k]=Complex(-k*U.im,k*U.re);  // I*kz*U

        f6[i][j][k]=Complex(-i*V.im,i*V.re);  // I*kx*V
        f7[i][j][k]=Complex(-j*V.im,j*V.re);  // I*ky*V
        f8[i][j][k]=Complex(-k*V.im,k*V.re);  // I*kz*V

        f9[i][j][k]=Complex(-i*W.im,i*W.re);  // I*kx*W
        f10[i][j][k]=Complex(-j*W.im,j*W.re); // I*ky*W
      }
    } 
  }
  
  Complex *F[]={f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10};
  Convolution->convolve(F,multadvection3);
  
  vector3 S0=S[0];
  vector3 S1=S[1];
  vector3 S2=S[2];
  
  // The purpose of pressure is to enforce incompressibility!
  // Apply projection operator.
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        Complex s0=f0[i][j][k];
        Complex s1=f1[i][j][k];
        Complex s2=f2[i][j][k];
        // Calculate -i*P
        Complex miP=(i*s0+j*s1+k*s2)/(i*i+j*j+k*k);
        S0[i][j][k]=i*miP-s0;
        S1[i][j][k]=j*miP-s1;
        S2[i][j][k]=k*miP-s2;
      }
    }
  }
  
/*  
    for(int i=0; i < mx; ++i) {
    for(int j=0; j < my; ++j) {
    for(int k=0; k < mz; ++k) {
    double k2=i*i+j*j+k*k;
    f0[i][j][k] += -nu*k2*u[i][j][k];
    }
    }
    }
*/

  S0(0,0,0)=0.0; // Enforce no mean flow.
  S1(0,0,0)=0.0;
  S2(0,0,0)=0.0;
  
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,S0);
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,S1);
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,S2);
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
