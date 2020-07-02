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

#ifdef __SSE2__
static inline Vec sqrt(const Vec& z)
{
  return _mm_sqrt_pd(z);
}
#else
static inline Vec sqrt(const Vec& z)
{
  return Vec(sqrt(z.x),sqrt(z.y));
}
#endif

// 2D Navier-Stokes advection a la Basdevant with Smagorinsky subgrid model
// requiring only 4 inputs and 2 outputs.
double Cd;

void multSmagorinsky2(double **F, unsigned int m,
                    const unsigned int indexsize,
                    const unsigned int *index,
                    unsigned int r, unsigned int threads)
{
  double* F0=F[0];
  double* F1=F[1];
  double* F2=F[2];
  double* F3=F[3];

#ifdef __SSE2__
  unsigned int m1=m-1;
  PARALLEL(
    for(unsigned int j=0; j < m1; j += 2) {
      double *F0j=F0+j;
      double *F1j=F1+j;
      double *F2j=F2+j;
      double *F3j=F3+j;
      Vec u=LOAD(F0j);
      Vec v=LOAD(F1j);
      Vec ux=LOAD(F2j);
      Vec s12=LOAD(F3j);
      STORE(F0j,v*v-u*u+4*Cd*Cd*sqrt(ux*ux+s12*s12)*ux); // B(x,t)
      STORE(F1j,u*v);
    }
    );
  if(m % 2) {
    double u=F0[m1];
    double v=F1[m1];
    double ux=F2[m1];
    double s12=F3[m1];
    F0[m1]=v*v-u*u+4*Cd*Cd*sqrt(ux*ux+s12*s12)*ux; // B(x,t)
    F1[m1]=u*v;
  }
#else
  for(unsigned int j=0; j < m; ++j) {
    double u=F0[j];
    double v=F1[j];
    double ux=F2[j];
    double s12=F3[j];
    F0[j]=v*v-u*u+4*Cd*Cd*sqrt(ux*ux+s12*s12)*ux; // B(x,t)
    F1[j]=u*v;
  }
#endif
}

void Source(const vector2& w, vector2 &S)
{
  f0[0][0]=0.0; // Enforce no mean flow.
  f1[0][0]=0.0;
  f2[0][0]=0.0;
  f3[0][0]=0.0;

  // This 2D version requires only 6 FFTs per stage (in the spirit
  // of Basdevant, J. Comp. Phys, 50, 1983).
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k2inv=1.0/(i*i+j*j);
      double jk2inv=j*k2inv;
      double ik2inv=i*k2inv;
      Complex wij=w[i][j];
      Complex u=Complex(-wij.im*jk2inv,wij.re*jk2inv);
      Complex v=Complex(wij.im*ik2inv,-wij.re*ik2inv);
      f0[i][j]=u;
      f1[i][j]=v;
      f2[i][j]=Complex(-i*u.im,i*u.re); // dudx
      f3[i][j]=Complex(-i*v.im,i*v.re)+Complex(-j*u.im,j*u.re); // dvdx + dudy
    }
  }
q
  Complex *F[]={f0,f1,f2,f3};
  Convolution->convolve(F,multSmagorinsky2);

  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      f0[i][j]=i*j*f0[i][j]+(i*i-j*j)*f1[i][j]-nu*(i*i+j*j)*w[i][j];
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

  Cd=Cs*delta;

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
