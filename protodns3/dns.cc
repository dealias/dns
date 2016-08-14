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

double dt=1.0e-8;
double nu=0.0; // kinematic viscosity

int mx;
int my;
int mz;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;
typedef Array3<Complex> vector3;
typedef Array4<Complex> vector4;

vector4 u;
vector3 f0,f1,f2,f3,f4,f5;

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution3 *Convolution;

void init(vector4& u)
{
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        Complex val=1.0/(i*i+j*j+k*k);
        Complex U=val;
        Complex V=val;
        Complex W=val;
        if(k != 0)
          W=-(i*U+j*V)/k;
        else if(j != 0)
          V=-(i*U)/j;
        else
          U=0.0;
        u[0][i][j][k]=U;
        u[1][i][j][k]=V;
        u[2][i][j][k]=W;
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
  
  for(unsigned int j=0; j < m; ++j) {
    double u=F0[j];
    double v=F1[j];
    double w=F2[j];
    F0[j]=u*u;
    F1[j]=u*v;
    F2[j]=u*w;
    F3[j]=v*v;
    F4[j]=v*w;
    F5[j]=w*w;
  }
}

void Source(const vector4& u, vector4 &S)
{
  f0[0][0][0]=0.0;
  f1[0][0][0]=0.0;
  f2[0][0][0]=0.0;
  
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        f0[i][j][k]=u[0][i][j][k];
        f1[i][j][k]=u[1][i][j][k];
        f2[i][j][k]=u[2][i][j][k];
      }
    } 
  }

  Complex *F[]={f0,f1,f2,f3,f4,f5};
  Convolution->convolve(F,multadvection3);
  
  vector3 S0=S[0];
  vector3 S1=S[1];
  vector3 S2=S[2];
  
  S0(0,0,0)=0.0; // Enforce no mean flow.
  S1(0,0,0)=0.0;
  S2(0,0,0)=0.0;
  
  // The purpose of pressure is to enforce incompressibility!
  // Apply projection operator.
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        Complex S00=f0[i][j][k];
        Complex S01=f1[i][j][k];
        Complex S02=f2[i][j][k];
        Complex S11=f3[i][j][k];
        Complex S12=f4[i][j][k];
        Complex S22=f5[i][j][k];
        
        Complex s0=Complex(0.0,i*S00+j*S01+k*S02);
        Complex s1=Complex(0.0,i*S01+j*S11+k*S12);
        Complex s2=Complex(0.0,i*S02+j*S12+k*S22);
        
        // Calculate -i*P
        Complex miP=(i*s0+j*s1+k*s2)/(i*i+j*j+k*k);
        S0[i][j][k]=i*miP-s0;
        S1[i][j][k]=j*miP-s1;
        S2[i][j][k]=k*miP-s2;
      }
    }
  }
    
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,S0);
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,S1);
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,S2);
  
#if 0
  Complex sum=0.0;
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
//        sum += S0[i][j][k]*i+S1[i][j][k]*j+S2[i][j][k]*k;
//        sum += u[0][i][j][k]*i+u[1][i][j][k]*j+u[2][i][j][k]*k;
        sum += (S0[i][j][k]*conj(u[0][i][j][k])).re
          +(S1[i][j][k]*conj(u[1][i][j][k])).re
          +(S2[i][j][k]*conj(u[2][i][j][k])).re;
      }
    }
  }
  cout << "sum=" << sum << endl;
  cout << endl;
#endif  
  
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
        double nuk2=nu*(i*i+j*j+k*k);
        S0[i][j][k] -= nuk2*S0;
        S1[i][j][k] -= nuk2*S1;
        S2[i][j][k] -= nuk2*S2;
      }
    }
  }
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
  for(int K=0; K <= kmax; ++K) {
    E[K]=0.0;
    Z[K]=0.0;
  }
     
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=-my+1; j < my; ++j) {
      for(int k=(j < 0 || (j == 0 && i <= 0)) ? 1 : 0; k < mz; ++k) {
	int k2=i*i+j*j+k*k;
        int K=sqrt(k2);	
        int index=(int) (K+0.5);
        double e=abs2(u[0][i][j][k],u[1][i][j][k],u[2][i][j][k]);
        E[index] += e;
        Z[index] += k2*e;
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
    cout << endl;
  }
  ezvt << E << "\t" << Z << endl;
}

int main(int argc, char* argv[])
{
  int n;
  cout << "Number of time steps? " << endl;
  cin >> n;
  cout << endl;

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

  Convolution=new ImplicitHConvolution3(mx,my,mz,true,true,true,3,6);
  
  u.Allocate(3,Nx,Ny,mz,0,-mx+1,-my+1,0,align);
  S.Allocate(3,Nx,Ny,mz,0,-mx+1,-my+1,0,align);
  
  init(u);
  u(0,0,0,0)=0.0; // Enforce no mean flow.
  u(1,0,0,0)=0.0;
  u(2,0,0,0)=0.0;
  
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,u[0]);
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,u[1]);
  HermitianSymmetrizeXY(mx,my,mz,mx-1,my-1,u[2]);
  
  cout.precision(15);
  
  for(int step=0; step < n; ++step) {
//    Output(step,step == 0);
    Output(step,true);
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
