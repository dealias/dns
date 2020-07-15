#include <cmath>
#include <iomanip>
#include <fstream>
#include "Complex.h"
#include "convolution.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;


int Nx=127; // Number of modes in x direction
int Ny=127; // Number of modes in y direction

double Cs=0.1;
double delta=twopi/Nx;

double dt=1.0e-4;
double nu=0.002; // kinematic viscosity
double nuL=0.2;  // damping (uniform friction across all scales)

double kforce=10.0;
double deltaf=1.0;
Complex force(1.0,1.0);

int mx;
int my;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector2 w; // vorticity
vector2 f0,f1,f2,f3; // inputs and outputs to advection routine
vector2 F; // constant vorticity forccing

ofstream ezvt("ezvt",ios::out);

ImplicitHConvolution2 *Convolution;

void init(vector2& w, vector2 &F)
{
  for(int i=-mx+1; i < mx; ++i) {
    vector wi=w[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      wi[j]=1.0/(i2+j*j);
    }
  }

  for(int i=-mx+1; i < mx; ++i) {
    vector Fi=F[i];
    int i2=i*i;
    double halfdeltaf=0.5*deltaf;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k=sqrt(i2+j*j);
      Fi[j]=abs(k-kforce) < halfdeltaf ? force : 0.0;
    }
  }
}

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
      //double nut = Cd*Cd*SQRT(ux*ux+s12*s12);
      STORE(F0j,v*v-u*u+4*Cd*Cd*SQRT(ux*ux+s12*s12)*ux); // B(x,t)
      STORE(F1j,u*v-Cd*Cd*SQRT(ux*ux+s12*s12)*s12); // A(x,t)
    }
    );
  if(m % 2) {
    double u=F0[m1];
    double v=F1[m1];
    double ux=F2[m1];
    double s12=F3[m1];
    double nut = Cd*Cd*sqrt(ux*ux+s12*s12);
    F0[m1]=v*v-u*u+4*nut*ux; // B(x,t)
    F1[m1]=u*v-nut*s12;// A(x,t)
  }
#else
  for(unsigned int j=0; j < m; ++j) {
    double u=F0[j];
    double v=F1[j];
    double ux=F2[j];
    double s12=F3[j];
    double nut = Cd*Cd*sqrt(ux*ux+s12*s12);
    F0[j]=v*v-u*u+4*nut*ux; // B(x,t)
    F1[j]=u*v-nut*s12;// A(x,t)
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
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    vector f2i=f2[i];
    vector f3i=f3[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double k2inv=1.0/(i2+j*j);
      double jk2inv=j*k2inv;
      double ik2inv=i*k2inv;
      Complex wij=wi[j];
      Complex u=Complex(-wij.im*jk2inv,wij.re*jk2inv);
      Complex v=Complex(wij.im*ik2inv,-wij.re*ik2inv);
      f0i[j]=u;
      f1i[j]=v;
      f2i[j]=Complex(-i*u.im,i*u.re); // F{dudx}
      f3i[j]=Complex(-i*v.im,i*v.re)+Complex(-j*u.im,j*u.re); // F{dvdx + dudy}
    }
  }

  Complex *fAll[]={f0,f1,f2,f3};
  Convolution->convolve(fAll,multSmagorinsky2);

  for(int i=-mx+1; i < mx; ++i) {
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    vector Fi=F[i];
    int i2=i*i;
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      int j2=j*j;
        f0i[j]=i*j*f0i[j]+(i2-j2)*f1i[j]-nu*(i2+j2)*wi[j]-nuL*wi[j]+Fi[j];
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
  double t=step*dt;
  if(verbose) {
    cout << "t=" << t << endl;
    cout << "Energy=" << E << endl;
    cout << "Enstrophy=" << Z << endl;
    cout << "Palenstrophy=" << P << endl;
  }
  ezvt << t << "\t" << E << "\t" << Z << "\t" << P << endl;
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

  F.Allocate(Nx,my,-mx+1,0,align);

  Convolution=new ImplicitHConvolution2(mx,my,true,true,4,2);

  Cd=Cs*delta;

  w.Allocate(Nx,my,-mx+1,0,align);

  init(w,F);
  w[0][0]=0.0; // Enforce no mean flow.
  F[0][0]=0.0;

  cout.precision(15);
  int kmax=(int) hypot(mx-1,my-1);
  double aZ[kmax+1];
  for(int k=0; k <= kmax; ++k) aZ[k]=0.0;
  double counter = 0.0;
    
  for(int step=0; step < n; ++step) {
    //Output(step,step == 0);
     Source(w,f0);
     for(int i=-mx+1; i < mx; ++i) {
       vector wi=w[i];
       vector f0i=f0[i];
       for(int j=(i <= 0 ? 1 : 0); j < my; ++j)
	 wi[j] += f0i[j]*dt;
         
     }
     cout << "[" << step << "] " << flush;
     if (step >= (n*9)/10){
         counter += 1.0;
         for(int i=-mx+1; i < mx; ++i) {
           vector wi=w[i];
           int i2=i*i;
           for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
             int k=sqrt(i2+j*j);
             aZ[(int) (k+0.5)] += abs2(wi[j]);
           }
         }
     }
  }
  cout << endl;
  Output(n,true);
    
  ofstream av_zkvk("av_zkvk",ios::out);
    counter += 1.0;
  for(int i=-mx+1; i < mx; ++i) {
  vector wi=w[i];
  int i2=i*i;
  for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
    int k=sqrt(i2+j*j);
    aZ[(int) (k+0.5)] += abs2(wi[j]);
    }
  }
  av_zkvk << "# k\tZ(k)" << endl;
  cout << counter << "\n";
  for(int k=1; k <= kmax; ++k) {
      av_zkvk << k << "\t" << aZ[k]/counter << endl;
  }
  Spectrum();
  

  return 0;
}
