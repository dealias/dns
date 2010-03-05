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
  
const double E=exp(1.0);

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

inline void init(Complex *f, Complex *g) 
{
  if(Test) {
    for(unsigned int i=0; i < m; i++) f[i]=g[i]=i*E;
  } else {
    f[0]=1.0;
    for(unsigned int i=1; i < m; i++) f[i]=Complex(3.0,2.0);
    g[0]=2.0;
    for(unsigned int i=1; i < m; i++) g[i]=Complex(5.0,3.0);
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
    int c = getopt(argc,argv,"deiptN:m:");
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
    
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);

  Complex *h0=NULL;
  if(Test) h0=FFTWComplex(m);

  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(Implicit) {
    unsigned int c=m/2;
    Complex *u=FFTWComplex(c+1);
    Complex *v=FFTWComplex(c+1);
    ImplicitHConvolution C(m,u,v);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      C.convolve(f,g,u,v);
      sum += seconds();
    }
    FFTWdelete(v);
    FFTWdelete(u);
    
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
    Complex *h=FFTWComplex(m);
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
    FFTWdelete(h);
    if(Test) 
      for(unsigned int i=0; i < m; i++) h0[i]=h[i];
  }

  if(Test) {
    Complex *h=FFTWComplex(m);
    double error=0.0;
    cout << endl;
    cout << "Exact:" << endl;
    double norm=0.0;
    for(unsigned int k=0; k < m; k++) {
       h[k]=E*E*(4.0*m*m*m-6.0*(k+1)*m*m+(6.0*k+2.0)*m+3.0*k*k*k-3.0*k)/6.0;
      norm += abs2(h[k]);
      error += abs2(h[k]-h0[k]);
    }
    error=sqrt(error/norm);
    cout << "error=" << error << endl;
    if (error > 1e-12)
      cerr << "Caution! error=" << error << endl;
    FFTWdelete(h);
  }
  
  FFTWdelete(f);
  FFTWdelete(g);
}
