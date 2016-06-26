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
int Nz=15; // Number of modes in z direction

double dt=1.0e-3;
double nu=0.0; // kinematic viscosity

int mx;
int my;
int mz;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;
typedef Array3<Complex> vector3;
typedef Array4<Complex> vector4;

vector4 u;
vector3 f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10;

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution3 *Convolution;

void init(vector4& u)
{
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        Complex val=1.0/(i*i+j*j+k*k);
        for(int c=0; c < 3; ++c) {
          u[c][i][j][k]=val;
        }
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
    for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
    for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
    double k2=i*i+j*j+k*k;
    S0[i][j][k] += -nu*k2*S0;
    S1[i][j][k] += -nu*k2*S1;
    S2[i][j][k] += -nu*k2*S2;
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

inline double hypot(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
}

inline double abs2(Complex x, Complex y, Complex z)
{
  return abs2(x)+abs2(y)+abs2(z);
}

void Spectrum()
{
  ofstream ekvk("ekvk",ios::out);
  
  int kmax=(int) hypot(mx-1,my-1,mz-1);
  double E[kmax+1];
  double Z[kmax+1];
  for(int K=0; K <= kmax; ++K) Z[K]=0.0;
     
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        int K=sqrt(i*i+j*j+k*k);
        int index=(int) (K+0.5);
        double e=abs2(u[0][i][j][k],u[1][i][j][k],u[2][i][j][k]);
        E[index] += e;
        Z[index] += K*K*e;
      }
    }
  }
    
  ekvk << "# k\tE(k)" << endl;
  
  for(int k=1; k <= kmax; ++k)
    ekvk << k << "\t" << E[k] << "\t" << Z[k] << endl;
}

void Output(int step, bool verbose=false)
{
  double E=0.0, Z=0.0;
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
	double k2=i*i+j*j+k*k;
        double e=abs2(u[0][i][j][k],u[1][i][j][k],u[2][i][j][k]);
	E += e;
	Z += k2*e;
      }
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
  vector4 S;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  mz=(Nz+1)/2;
  size_t align=sizeof(Complex);

  f0.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f1.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f2.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f3.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f4.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f5.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f6.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f7.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f8.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f9.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);
  f10.Allocate(Nx,Ny,mz,-mx+1,-my+1,0,align);

  Convolution=new fftwpp::ImplicitHConvolution3(mx,my,mz,true,true,true,11,3);
  
  u.Allocate(3,Nx,Ny,mz,0,-mx+1,-my+1,0,align);
  S.Allocate(3,Nx,Ny,mz,0,-mx+1,-my+1,0,align);
  
  init(u);
  u(0,0,0,0)=0.0; // Enforce no mean flow.
  u(1,0,0,0)=0.0;
  u(2,0,0,0)=0.0;
  
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,u[0]);
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,u[1]);
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,u[2]);

  int n=10;
  
  cout.precision(15);
  
  for(int step=0; step < n; ++step) {
    Output(step,step == 0);
    Source(u,S);
    for(int i=-mx+1; i < mx; ++i) {
      for(int j=-my+1; j < my; ++j) {
	for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
          for(int c=0; c < 3; ++c) {
            u[c][i][j][k] += S[c][i][j][k]*dt;
          }
        }
      }
    }
     cout << "[" << step << "] ";
  }
  cout << endl;
  Output(n,true);
  Spectrum();
     
  return 0;
}
