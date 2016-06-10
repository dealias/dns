#include <cmath>
#include <iomanip>
#include <fstream>
#include "Complex.h"
#include "convolution.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;



// size of problem

int Nx=15;
int Ny=15;

int mx;
int my;
const double nu=0.0;

typedef Array1<Complex>::opt vector;
typedef Array2<Complex> vector2;

vector2 w;
vector2 f0,f1,g0,g1;

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
  
  for(int i=0; i<mx ;++i){
    for(int j=0; j<my ;++j){
       double k2=i*i+j*j;     
       f0[i][j] += -nu*k2*w[i][j];
    }
  }
  f0(0,0)=0.0; // Enforce no mean flow.
  
  fftwpp::HermitianSymmetrizeX(mx,my,mx-1,f0);
  
#if 1 // Check enstrophy and energy symmetries
  double sumE=0.0;
  double sumZ=0.0;
  for(int i=-mx+1; i < mx; ++i) {
    for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
      double product=real(f0[i][j]*conj(w[i][j]));
      sumE += product/(i*i+j*j);
      sumZ += product;
    }
  }
  cout << "sum wS=" << sumZ << endl;
  cout << "sum wS/k^2=" << sumE << endl;
  cout << endl;
#endif
}

int main(int argc, char* argv[])
{
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  const int kmax=(int) hypot(mx-1,my-1);
  double Z[kmax+1]={};
  double E[kmax+1]={};
  double Enstrophy=0,Energy=0;
  ofstream zout("Zk.dat",ios::out);
  ofstream eout("Ek.dat",ios::out);
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

  double dt=1.0e-6;
  for(int step=0; step < 100; ++step){
     Source(w,f0);
     for(int i=0; i < mx; ++i){
       for(int j=0; j < my; ++j){
	 w[i][j] += f0[i][j]*dt;
       }
     }
     cout << "step:" << step << endl;
     int k=0;
     for(int i=0; i < mx; ++i){
       for(int j=0; j < my; ++j){
	 k=sqrt(i*i+j*j);
	 Z[(int) (k+0.5)] += 0.5*abs(w[i][j])*abs(w[i][j]);
       }
     }
  }
  cout << endl << endl;
  cout << "k" << setw(15) << "Z[k]" << setw(15) << "E[k]" << endl;
  cout << "--------------------------------" << endl;  
  for(int i=1; i <= kmax; ++i){
    E[i]=Z[i]/(i*i);
    cout << i << setw(15) << Z[i] << setw(15) << E[i] << endl;
    zout << i << "\t" << Z[i] <<endl;
    eout << i << "\t" << E[i] <<endl;
    Energy +=2*E[i];
    Enstrophy +=2*Z[i];
  }
  cout << "--------------------------------" << endl;
  cout << "Energy=" << "\t\t" << Energy << endl;
  cout << "Enstrophy=" <<  "\t" << Enstrophy << endl;
  cout << "--------------------------------" << endl;  
  zout.close();
  eout.close();  
  return 0;
}
