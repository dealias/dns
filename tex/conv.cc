using namespace std;
#include "Complex.h"
#include "convolution.h"

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// Number of iterations.
unsigned int N0=10000000;
unsigned int N=0;
unsigned int m=12;
unsigned int M=1;
  
const Complex I(0.0,1.0);
const double E=exp(1.0);
const double F=sqrt(3.0);
const double G=sqrt(5.0);

bool Direct=false, Implicit=true, Explicit=false, Test=false;

using namespace std;

#include <sys/time.h>

inline double seconds()
{
  static timeval lasttime;
  timeval tv;
  gettimeofday(&tv,NULL);
  double seconds=tv.tv_sec-lasttime.tv_sec+
    ((double) tv.tv_usec-lasttime.tv_usec)/1000000.0;
  lasttime=tv;
  return seconds;
}

inline void init(Complex *f, Complex *g, int M=1) 
{
  unsigned int Mm=M*m;
  double factor=1.0/sqrt(M);
  for(unsigned int i=0; i < Mm; i += m) {
    Complex *fi=f+i;
    Complex *gi=g+i;
    if(Test) {
      for(unsigned int k=0; k < m; k++) fi[k]=factor*F*pow(E,k*I);
      for(unsigned int k=0; k < m; k++) gi[k]=factor*G*pow(E,k*I);
//    for(unsigned int k=0; k < m; k++) fi[k]=factor*F*k;
//    for(unsigned int k=0; k < m; k++) gi[k]=factor*G*k;
    } else {
      fi[0]=1.0*factor;
      for(unsigned int k=1; k < m; k++) fi[k]=factor*Complex(3.0,2.0);
      gi[0]=2.0*factor;
      for(unsigned int k=1; k < m; k++) gi[k]=factor*Complex(5.0,3.0);
    }
  }
}

int main(int argc, char* argv[])
{
#ifndef __SSE2__
  fftw::effort |= FFTW_NO_SIMD;
#endif  
  
#ifdef __GNUC__	
  optind=0;
#endif	
  for (;;) {
    int c = getopt(argc,argv,"deiptM:N:m:");
    if (c == -1) break;
		
    switch (c) {
      case 0:
        break;
      case 'd':
        Direct=true;
        break;
      case 'e':
        Explicit=true;
        Implicit=false;
        break;
      case 'i':
        Implicit=true;
        Explicit=false;
        break;
      case 'p':
        break;
      case 'M':
        M=atoi(optarg);
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 't':
        Test=true;
        break;
      case 'm':
        m=atoi(optarg);
        break;
    }
  }

  unsigned int n=3*m-2;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  if(N == 0) {
    N=N0/n;
    if(N < 10) N=10;
  }
  cout << "N=" << N << endl;
  
  unsigned int np=Explicit ? n/2+1 : m;
  if(Implicit) np *= M;
    
  Complex *f=ComplexAlign(np);
  Complex *g=ComplexAlign(np);

  Complex *h0=NULL;
  if(Test) h0=ComplexAlign(m);

  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(Implicit) {
    ImplicitHConvolution C(m,M);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g,M);
      seconds();
      C.convolve(f,g);
      sum += seconds();
    }
     
    cout << endl;
    cout << "Implicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << f[i] << endl;
    else cout << f[0] << endl;
    if(Test)
      for(unsigned int i=0; i < m; i++) h0[i]=f[i];
  }
  
  if(Explicit) {
    sum=0.0;
    ExplicitHConvolution C(n,m,f);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      C.convolve(f,g);
      sum += seconds();
    }
    
    cout << endl;
    cout << "Explicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << f[i] << endl;
    else cout << f[0] << endl;
    if(Test) 
      for(unsigned int i=0; i < m; i++) h0[i]=f[i];
  }
  
  if(Direct) {
    DirectHConvolution C(m);
    init(f,g);
    Complex *h=ComplexAlign(m);
    seconds();
    C.convolve(h,f,g);
    sum=seconds();
    
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
    if(Test) 
      for(unsigned int i=0; i < m; i++) h0[i]=h[i];
    deleteAlign(h);
  }

  if(Test) {
    Complex *h=ComplexAlign(m);
    double error=0.0;
    cout << endl;
    double norm=0.0;
    long long M=m;
    for(long long k=0; k < M; k++) {
      h[k]=F*G*(2*M-1-k)*pow(E,k*I);
//      h[k]=F*G*(4*m*m*m-6*(k+1)*m*m+(6*k+2)*m+3*k*k*k-3*k)/6.0;
      error += abs2(h0[k]-h[k]);
      norm += abs2(h[k]);
    }
    error=sqrt(error/norm);
    cout << "error=" << error << endl;
    if (error > 1e-12)
      cerr << "Caution! error=" << error << endl;
    deleteAlign(h);
  }
  
  deleteAlign(f);
  deleteAlign(g);
}
